#include "clock_cn.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_CN::Clock_CN(MAPP* mapp,int narg
                   ,char** arg):Clock(mapp)
{


    min_gamma=1.0e-10;
    gamma_red=0.8;
    slope=0.4;
    max_iter=50;
    m_tol=1.0e-5;
    a_tol=1.0e-10;
    min_del_t=1.0e-12;
    max_del_t=1.0e4;
    initial_del_t=-1.0;

    
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
            if(strcmp(arg[iarg],"min_gamma")==0)
            {
                iarg++;
                min_gamma=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"red_gamma")==0)
            {
                iarg++;
                gamma_red=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"slope")==0)
            {
                iarg++;
                slope=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_iter")==0)
            {
                iarg++;
                max_iter=atoi(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"m_tol")==0)
            {
                iarg++;
                m_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"a_tol")==0)
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
            else if(strcmp(arg[iarg],"initial_del_t")==0)
            {
                iarg++;
                initial_del_t=atof(arg[iarg]);
                if(initial_del_t<=0.0)
                    error->abort("initial_del_t should be larger than zero");
                iarg++;
            }
            else
                error->abort("unknown keyword: %s",arg[iarg]);
        }
    }
    
    if(min_gamma<=0.0)
        error->abort("min_gamma should be larger than zero");
    if(gamma_red<=0.0 || gamma_red>=1.0)
        error->abort("gamma_red should be between zero and one");
    if(slope<=0.0 || slope>=1.0)
        error->abort("slope should be between zero and one");
    if(max_iter<=0)
        error->abort("max_iter should be larger than zero");
    if(m_tol<=0.0)
        error->abort("m_tol should be larger than zero");
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
    
    CREATE1D(t,2);
    CREATE1D(dy,2);
    for(int i=0;i<2;i++)
        CREATE1D(dy[i],dof_lcl);
    
    CREATE1D(y0,dof_lcl);
    CREATE1D(dy0,dof_lcl);
    CREATE1D(y_0,dof_lcl);
    CREATE1D(a,dof_lcl);
    CREATE1D(g,dof_lcl);
    CREATE1D(g0,dof_lcl);
    CREATE1D(c0,dof_lcl);
    CREATE1D(h,dof_lcl);

    

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_CN::~Clock_CN()
{
    delete thermo;
    
    if(dof_lcl)
    {
        for(int i=0;i<2;i++)
            delete [] dy[i];
        
        delete [] y0;
        delete [] dy0;
        delete [] y_0;
        delete [] a;
        delete [] g0;
        delete [] c0;
        delete [] g;
        delete [] h;
    }
    
    delete [] t;
    delete [] dy;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_CN::init()
{
    TYPE0* c;
    TYPE0* c_d;
    
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<TYPE0>(0,atoms->vectors[0].dim,"f");
    
    vecs_comm=new VecLst(mapp,3,0,c_n,c_d_n);
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
    
    t[0]=0.0;
    thermo->start_force_time();
    forcefield->c_d_calc();
    thermo->stop_force_time();
    atoms->vectors[c_n].ret(c);
    atoms->vectors[c_d_n].ret(c_d);
    
    rectify(c_d);
    if(initial_del_t<0.0)
    {
        TYPE0 sum,sum_lcl;
        sum_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
            sum_lcl+=c_d[i]*c_d[i];
        sum=0.0;
        MPI_Allreduce(&sum_lcl,&sum,1,MPI_TYPE0,MPI_SUM,world);
        sum=sqrt(sum/static_cast<TYPE0>(dof_tot));
        initial_del_t=MAX(a_tol/sum,min_del_t);
        initial_del_t=MIN(initial_del_t,max_del_t);
    }
    
    t[1]=-initial_del_t;
    memcpy(dy[0],c_d,dof_lcl*sizeof(TYPE0));
    memcpy(dy[1],c_d,dof_lcl*sizeof(TYPE0));
    memcpy(y0,c,dof_lcl*sizeof(TYPE0));
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_CN::fin()
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
 init
 --------------------------------------------*/
void Clock_CN::run()
{
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    int chk;
    TYPE0 del_t,cost,err1,ratio,del_t_tmp;
    TYPE0* tmp_dy;
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    del_t=initial_del_t;
    
    for(int istep=0;istep<no_steps;istep++)
    {
        chk=interpolate(del_t);
        while(chk==-1)
        {
            del_t*=0.8;
            del_t=MAX(del_t,min_del_t);
            //cout << del_t <<endl;
            chk=interpolate(del_t);
        }
        cost=solve(del_t);
        while(err>=1.0 || cost>=1.0)
        {
            
            err1=MAX(err,cost);
            ratio=pow(0.5/err1,1.0/3.0);
            
            if(ratio<0.5)
                ratio=0.5;
            else if(ratio>0.9)
                ratio=0.9;
            
            del_t_tmp=del_t*ratio;
            if(del_t_tmp<min_del_t)
                del_t=min_del_t;
            else
                del_t=del_t_tmp;
            
            
            chk=interpolate(del_t);
            
            while(chk==-1)
            {
                del_t*=0.8;
                chk=interpolate(del_t);
            }
            
            cost=solve(del_t);
        }
        
        //printf("err %e\n",err);
        
        t[1]=t[0];
        t[0]+=del_t;
        tmp_dy=dy[1];
        dy[1]=dy[0];
        dy[0]=tmp_dy;
        memcpy(dy[0],c_d,dof_lcl*sizeof(TYPE0));
        memcpy(y0,c,dof_lcl*sizeof(TYPE0));
        
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step())
        {
            forcefield->force_calc(1,energy_stress);
            thermo->update(fe_idx,energy_stress[0]);
            thermo->update(stress_idx,6,&energy_stress[1]);
            thermo->update(time_idx,t[0]);
        }
        
        
        
        ratio=pow(0.5/err,1.0/3.0);
        
        if(ratio<0.5)
            ratio=0.5;
        else if(ratio>2.0)
            ratio=2.0;
        
        del_t_tmp=del_t*ratio;
        if(del_t_tmp<min_del_t)
            del_t=min_del_t;
        else
            del_t=del_t_tmp;
        
        
        step_no++;
    }
    delete [] energy_stress;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
int Clock_CN::interpolate(TYPE0 del_t)
{
    TYPE0 c0=del_t+0.5*del_t*del_t/(t[0]-t[1]);
    TYPE0 c1=-0.5*del_t*del_t/(t[0]-t[1]);

    
    int ret_val=1;
    int tot_ret_val;
    int idof=0;
    while(idof<dof_lcl && ret_val==1)
    {
        y_0[idof]=y0[idof]+c0*dy[0][idof]+c1*dy[1][idof];
        a[idof]=-y0[idof]-0.5*del_t*dy[0][idof];
        if(y_0[idof]<0.0 || y_0[idof]>1.0)
        {
            ret_val=-1;
        }
        idof++;
    }
    beta=del_t*0.5;
    
    MPI_Allreduce(&ret_val,&tot_ret_val,1,MPI_INT,MPI_MIN,world);
    
    return tot_ret_val;
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
TYPE0 Clock_CN::solve(TYPE0 del_t)
{
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    
    TYPE0 gamma,max_gamma=1.0;
    TYPE0 inner,tmp;
    TYPE0 ratio;
    TYPE0 g0_g0,g_g,g_g0,g_h;
    TYPE0 curr_cost,ideal_cost,cost;
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    int chk;

    
    memcpy(c,y_0,dof_lcl*sizeof(TYPE0));
    thermo->start_comm_time();
    atoms->update(c_n);
    thermo->stop_comm_time();
    
    curr_cost=forcefield->calc_g(0,beta,a,g);
    rectify(g);
    memcpy(h,g,dof_lcl*sizeof(TYPE0));
    
    inner=0.0;
    for(int i=0;i<dof_lcl;i++)
        inner+=g[i]*g[i];
    g0_g0=0.0;
    MPI_Allreduce(&inner,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    g_h=g0_g0;
    
    int iter=0;
    
    while(curr_cost>m_tol*static_cast<TYPE0>(dof_tot)
          && iter<max_iter && max_gamma>min_gamma)
    {
        memcpy(g0,g,dof_lcl*sizeof(TYPE0));
        memcpy(c0,c,dof_lcl*sizeof(TYPE0));
        
        gamma=0.99;
        
        for(int i=0;i<dof_lcl;i++)
        {
            tmp=c0[i]+h[i];
            
            if(tmp>1.0)
            {
                gamma=MIN((1.0-c0[i])/(tmp-c0[i]),gamma);
            }
            else if(tmp<0.0)
            {
                gamma=MIN((c0[i]-0.0)/(c0[i]-tmp),gamma);
            }
        }
        MPI_Allreduce(&gamma,&max_gamma,1,MPI_TYPE0,MPI_MIN,world);
        
        max_gamma*=0.9;
        chk=1;
        
        cost=curr_cost;
        while(chk && max_gamma>min_gamma)
        {
            for(int i=0;i<dof_lcl;i++)
                c[i]=c0[i]+max_gamma*h[i];
            
            thermo->start_comm_time();
            atoms->update(c_n);
            thermo->stop_comm_time();
            
            thermo->start_force_time();
            curr_cost=forcefield->calc_g(0,beta,a,g);
            rectify(g);
            thermo->stop_force_time();
            
            ideal_cost=cost-slope*max_gamma*g_h;
            if(curr_cost<ideal_cost)
                chk=0;
            max_gamma*=gamma_red;
            if(max_gamma<=min_gamma)
            {
                memcpy(c,c0,dof_lcl*sizeof(TYPE0));
                
                thermo->start_comm_time();
                atoms->update(c_n);
                thermo->stop_comm_time();
                
                thermo->start_force_time();
                curr_cost=forcefield->calc_g(1,beta,a,g);
                thermo->stop_force_time();
            }
        }
        
        if(chk==0)
        {
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                inner+=g[i]*g0[i];
            g_g0=0.0;
            MPI_Allreduce(&inner,&g_g0,1,MPI_TYPE0,MPI_SUM,world);
            
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                inner+=g[i]*g[i];
            g_g=0.0;
            MPI_Allreduce(&inner,&g_g,1,MPI_TYPE0,MPI_SUM,world);
            
            ratio=(g_g-g_g0)/g0_g0;
            
            g0_g0=g_g;
            
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
            {
                h[i]*=ratio;
                h[i]+=g[i];
                inner+=h[i]*g[i];
            }
            g_h=0.0;
            MPI_Allreduce(&inner,&g_h,1,MPI_TYPE0,MPI_SUM,world);
            
            if(g_h<0.0)
            {
                memcpy(h,g,dof_lcl*sizeof(TYPE0));
                g_h=g_g;
            }
        }
        
        
        iter++;
    }
    
    
    rectify(c_d);

    
    TYPE0 err_lcl=0.0;
    TYPE0 tmp0;
    err=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp0=(y_0[i]-c[i]);
        err_lcl+=tmp0*tmp0;
    }
    
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    err=(1.0/3.0)*(del_t/(del_t+t[0]-t[1]))*sqrt(err/static_cast<TYPE0>(dof_tot))/a_tol;
    
    return curr_cost/(m_tol*static_cast<TYPE0>(dof_tot));
}
