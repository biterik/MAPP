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
    if(narg!=4)
        error->abort("invalid arguments");
    
    no_steps=atoi(arg[2]);
    delta_t=atof(arg[3]);
    delta_c=1.0e-6;
    min_gamma=1.0e-30;
    gamma_red=0.8;
    slope=0.4;
    MAX_ITER=50;
    RTOL=1.0e-8;
    ATOL=1.0e-4;
    min_del_t=1.0e-12;
    

    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n].dim;
    tot_dim=atoms->natms*c_dim;
    MPI_Allreduce(&tot_dim,&tot_dof,1,MPI_INT,MPI_SUM,world);
    
    CREATE1D(y_0,tot_dim);
    CREATE1D(a,tot_dim);
    CREATE1D(g,tot_dim);
    CREATE1D(g0,tot_dim);
    CREATE1D(c0,tot_dim);
    CREATE1D(h,tot_dim);

    

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_CN::~Clock_CN()
{
    delete thermo;
    int natms=atoms->natms;
    
    if(natms)
    {
        delete [] y_0;
        delete [] a;
        delete [] g0;
        delete [] c0;
        delete [] g;
        delete [] h;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_CN::init()
{
    TYPE0* c;
    
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
    
    atoms->vectors[c_n].ret(c);
    
    forcefield->c_d_calc();
    
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
    TYPE0 err,new_del_t;
    TYPE0 curr_t=0.0;
    
    for(int istep=0;istep<no_steps;istep++)
    {

        //err=1.0;
        //while (err>=1.0) {
            prepare();
            err=solve();
            
          //  if(err>=1.0)
            //    delta_t/=err;
        //}
        
        curr_t+=delta_t;
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step())
        {
            forcefield->force_calc(1,energy_stress);
            thermo->update(fe_idx,energy_stress[0]);
            thermo->update(stress_idx,6,&energy_stress[1]);
            thermo->update(time_idx,curr_t);
        }
        
        new_del_t=step_size();
        if(new_del_t>2.0*delta_t)
            delta_t*=2.0;
        
        step_no++;
    }
    delete [] energy_stress;
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
TYPE0 Clock_CN::step_size()
{
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    TYPE0 sug_del_t,sum0,tot_sum0;
    sum0=0.0;
    for(int i=0;i<tot_dim;i++)
        sum0+=c_d[i]*c_d[i];
    
    tot_sum0=0.0;
    MPI_Allreduce(&sum0,&tot_sum0,1,MPI_TYPE0,MPI_MIN,world);
    tot_sum0=sqrt(tot_sum0)/static_cast<TYPE0>(tot_dof);
    sug_del_t=delta_c/tot_sum0;

    return sug_del_t;
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
void Clock_CN::prepare()
{
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    TYPE0 ratio,tot_ratio,tmp;
    
    ratio=1.0;
    for(int i=0;i<tot_dim;i++)
    {
        tmp=c_d[i]*delta_t+c[i];
        if(tmp>1.0)
            ratio=MIN((1.0-c[i])/(tmp-c[i]),ratio);
        else if(tmp<0.0)
            ratio=MIN((c[i]-0.0)/(c[i]-tmp),ratio);
    }
    MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
    delta_t*=ratio;
    
    alpha=delta_t*0.5;
    for(int i=0;i<tot_dim;i++)
    {
        a[i]=-c[i]-0.5*delta_t*c_d[i];
        c[i]+=delta_t*c_d[i];
    }
    atoms->update(c_n);
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
TYPE0 Clock_CN::solve()
{
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    
    TYPE0 gamma,max_gamma=1.0;
    TYPE0 inner,tmp;
    TYPE0 ratio;
    TYPE0 g0_g0,g_g,g_g0,g_h;
    TYPE0 curr_cost,ideal_cost,cost;
    
    
    int chk;

    curr_cost=forcefield->calc_g(0,alpha,a,g);
    
    memcpy(h,g,tot_dim*sizeof(TYPE0));
    
    inner=0.0;
    for(int i=0;i<tot_dim;i++)
        inner+=g[i]*g[i];
    g0_g0=0.0;
    MPI_Allreduce(&inner,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    g_h=g0_g0;
    
    int iter=0;
    
    while(curr_cost>RTOL && iter<MAX_ITER && max_gamma>min_gamma)
    {
        //cout << "cost " << curr_cost<<endl;
        memcpy(g0,g,tot_dim*sizeof(TYPE0));
        memcpy(c0,c,tot_dim*sizeof(TYPE0));
        
        gamma=1.0;
        
        for(int i=0;i<tot_dim;i++)
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
        
        chk=1;
        
        cost=curr_cost;
        while(chk && max_gamma>min_gamma)
        {
            for(int i=0;i<tot_dim;i++)
                c[i]=c0[i]+max_gamma*h[i];
            
            atoms->update(c_n);
            
            curr_cost=forcefield->calc_g(1,alpha,a,g);
            ideal_cost=cost-slope*max_gamma*g_h;
            if(curr_cost<ideal_cost)
                chk=0;
            max_gamma*=gamma_red;
            if(max_gamma<=min_gamma)
            {
                memcpy(c,c0,tot_dim*sizeof(TYPE0));
                atoms->update(c_n);
                curr_cost=forcefield->calc_g(0,alpha,a,g);
            }
        }
        
        if(chk==0)
        {
            inner=0.0;
            for(int i=0;i<tot_dim;i++)
                inner+=g[i]*g0[i];
            g_g0=0.0;
            MPI_Allreduce(&inner,&g_g0,1,MPI_TYPE0,MPI_SUM,world);
            
            inner=0.0;
            for(int i=0;i<tot_dim;i++)
                inner+=g[i]*g[i];
            g_g=0.0;
            MPI_Allreduce(&inner,&g_g,1,MPI_TYPE0,MPI_SUM,world);
            
            ratio=(g_g-g_g0)/g0_g0;
            
            g0_g0=g_g;
            
            inner=0.0;
            for(int i=0;i<tot_dim;i++)
            {
                h[i]*=ratio;
                h[i]+=g[i];
                inner+=h[i]*g[i];
            }
            g_h=0.0;
            MPI_Allreduce(&inner,&g_h,1,MPI_TYPE0,MPI_SUM,world);
            
            if(g_h<0.0)
            {
                memcpy(h,g,tot_dim*sizeof(TYPE0));
                g_h=g_g;
            }
        }
        
        
        iter++;
    }
    
    return curr_cost/(RTOL*static_cast<TYPE0>(tot_dof));
}
