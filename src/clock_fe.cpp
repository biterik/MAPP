#include "clock_fe.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
#define SMALL 1.0e-18
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_FE::Clock_FE(MAPP* mapp,int narg
    ,char** arg):Clock(mapp)
{
    e_tol=0.0;
    a_tol=1.0e-6;
    min_del_t=1.0e-12;
    max_del_t=1.0e4;
    
    if(narg<3)
        error->abort("invalid arguments");
    no_steps=atoi(arg[2]);
    if(no_steps<=0)
        error->abort("no of steps should "
                     "be greater than zero");
    
    if(narg>3)
    {
        if((narg-3)%2!=0)
            error->abort("wrong number of inputs");
        int iarg=3;
        while(iarg<narg)
        {

            if(strcmp(arg[iarg],"a_tol")==0)
            {
                iarg++;
                a_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"e_tol")==0)
            {
                iarg++;
                e_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"min_del_t")==0)
            {
                iarg++;
                min_del_t=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_del_t")==0)
            {
                iarg++;
                max_del_t=atof(arg[iarg]);
                iarg++;
            }
            else
                error->abort("unknown keyword: %s",arg[iarg]);
        }
    }
    

    if(a_tol<=0.0)
        error->abort("a_tol should be larger than zero");
    if(e_tol<0.0)
        error->abort("e_tol should be larger than zero");
    if(min_del_t<=0.0)
        error->abort("min_del_t should be larger than zero");
    if(max_del_t<=0.0)
        error->abort("max_del_t should be larger than zero");
    if(min_del_t>=max_del_t)
        error->abort("max_del_t should be larger than min_del_t");
    
    
   
    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n].dim;
    dof_lcl=atoms->natms*c_dim;
    MPI_Allreduce(&dof_lcl,&dof_tot,1,MPI_INT,MPI_SUM,world);
    
    CREATE1D(y_0,dof_lcl);
    CREATE1D(dy_0,dof_lcl);
    CREATE1D(dy_1,dof_lcl);
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_FE::~Clock_FE()
{
    delete thermo;
    
    if(dof_lcl)
    {
        delete [] y_0;
        delete [] dy_0;
        delete [] dy_1;
    }
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_FE::init()
{
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<TYPE0>(0,atoms->vectors[0].dim,"f");
    
    c_n=atoms->find_exist("c");
    c_d_n=atoms->find("c_d");
    
    if(cdof_n>-1)
    {
        int dof_n=atoms->find("dof");
        vecs_comm=new VecLst(mapp,5,0,c_n,c_d_n,cdof_n,dof_n);

    }
    else
    {
        vecs_comm=new VecLst(mapp,3,0,c_n,c_d_n);
    }
    
    vecs_comm->add_update(0);
    atoms->reset_comm(vecs_comm);
    
    
    forcefield->init();
    atoms->ph_setup(1,vecs_comm);
    neighbor->init();
    neighbor->create_list(0,1);
    atoms->store_0();
    
    
    forcefield->create_2nd_neigh_lst();
    
    
    
    
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    forcefield->force_calc(1,energy_stress);
    
    thermo->update(fe_idx,energy_stress[0]);
    thermo->update(stress_idx,6,&energy_stress[1]);
    thermo->update(time_idx,0.0);
    
    
    if(write!=NULL)
        write->init();
    thermo->init();
    
    delete [] energy_stress;
    
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    memcpy(y_0,c,dof_lcl*sizeof(TYPE0));
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_FE::fin()
{
    if(write!=NULL)
        write->fin();
    
    forcefield->fin();
    neighbor->fin();
    
    delete vecs_comm;
    
    thermo->fin();
    atoms->x2s(atoms->natms);
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void Clock_FE::run()
{
    TYPE0 curr_t=0.0;
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    
    TYPE0 del_t,err_lcl,err_tot,tmp0,tmp1;
    TYPE0 tot_ratio,ratio;
    
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    del_t=min_del_t;
    eq_ratio=1.0;
    int istep=0;
    
    while (eq_ratio>=1.0 && istep <no_steps)
    {
        thermo->start_force_time();
        forcefield->c_d_calc();
        rectify(c_d);
        thermo->stop_force_time();
        
        memcpy(y_0,c,dof_lcl*sizeof(TYPE0));
        memcpy(dy_0,c_d,dof_lcl*sizeof(TYPE0));
        
        ratio=1.0;
        for(int i=0;i<dof_lcl;i++)
        {
            tmp0=y_0[i]+dy_0[i]*del_t;
            if(tmp0>1.0)
                ratio=MIN((1.0-y_0[i])/(tmp0-y_0[i]),ratio);
            else if(tmp0<0.0)
                ratio=MIN((y_0[i]-0.0)/(y_0[i]-tmp0),ratio);
        }
        
        MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
        del_t*=ratio;
        
        err_tot=1.0;
        while (err_tot>=1.0)
        {
            for(int i=0;i<dof_lcl;i++)
                c[i]+=0.5*del_t*c_d[i];
            
            thermo->start_comm_time();
            atoms->update(c_n);
            thermo->stop_comm_time();
            
            thermo->start_force_time();
            forcefield->c_d_calc();
            rectify(c_d);
            thermo->stop_force_time();
            
            memcpy(dy_1,c_d,dof_lcl*sizeof(TYPE0));
            err_lcl=0.0;
            for(int i=0;i<dof_lcl;i++)
            {
                tmp0=dy_1[i]-dy_0[i];
                err_lcl+=tmp0*tmp0;
            }
            err_tot=0.0;
            MPI_Allreduce(&err_lcl,&err_tot,1,MPI_TYPE0,MPI_SUM,world);
            err_tot=0.5*del_t*sqrt(err_tot/static_cast<TYPE0>(dof_tot))/a_tol;
            
            
            if(del_t==min_del_t && err_tot>=1.0)
                err_tot=0.0;
            
            if(err_tot>=1.0)
            {
                del_t*=0.9/err_tot;
                if (ratio<=0.5)
                    ratio=0.5;
                del_t*=ratio;
                del_t=MAX(del_t,min_del_t);
            }
        }
        
        for(int i=0;i<dof_lcl;i++)
            c[i]=y_0[i]+del_t*dy_0[i];
        
        thermo->start_comm_time();
        atoms->update(c_n);
        thermo->stop_comm_time();
        
        curr_t+=del_t;
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step())
        {
            thermo->start_force_time();
            forcefield->force_calc(1,energy_stress);
            thermo->stop_force_time();
            thermo->update(fe_idx,energy_stress[0]);
            thermo->update(stress_idx,6,&energy_stress[1]);
            thermo->update(time_idx,curr_t);
        }
        
        //cout<< err_tot <<endl;
        
        
        ratio=0.9/err_tot;
        
        if(ratio>=2.0)
            ratio=2.0;
        else if (ratio<=0.5)
            ratio=0.5;
        
        del_t*=ratio;
        del_t=MAX(del_t,min_del_t);
        del_t=MIN(del_t,max_del_t);

        
        
        tmp1=0.0;
        for(int i=0;i<dof_lcl;i++)
        {
            tmp1+=dy_0[i]*dy_0[i];
        }
        
        MPI_Allreduce(&tmp1,&eq_ratio,1,MPI_TYPE0,MPI_SUM,world);
        eq_ratio=sqrt(eq_ratio/static_cast<TYPE0>(dof_tot))/e_tol;
        
        istep++;
        step_no++;
    }
    
    thermo->start_force_time();
    forcefield->force_calc(1,energy_stress);
    thermo->stop_force_time();
    thermo->update(fe_idx,energy_stress[0]);
    thermo->update(stress_idx,6,&energy_stress[1]);
    delete [] energy_stress;
    
    
    
}

