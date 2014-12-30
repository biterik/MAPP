
#include "clock_adams.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
#define SMALL 1.0e-18
#define MAX_ITER 50
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_Adams::Clock_Adams(MAPP* mapp,int narg
                     ,char** arg):Clock(mapp)
{
    if(narg!=5)
        error->abort("invalid arguments");
    
    order=atoi(arg[2]);
    no_steps=atoi(arg[3]);
    delta_t=atof(arg[4]);
    
    
    gamma_red=0.8;
    slope=0.4;
    
    
    if(order>6 || order<1)
        error->abort("BDF cannot operate with order: %d",order);
    
    CREATE1D(beta,order);
    CREATE1D(t,order);
    CREATE1D(beta_mod,order);
    
    if(order==1)
    {
        beta[0]=1.0;
    }
    else if(order==2)
    {
        beta[0]=1.0/2.0;
        beta[1]=1.0/2.0;
    }
    else if(order==3)
    {
        beta[0]=5.0/12.0;
        beta[1]=8.0/12.0;
        beta[2]=-1.0/12.0;
    }
    else if(order==4)
    {
        beta[0]=9.0/24.0;
        beta[1]=19.0/24.0;
        beta[2]=-5.0/24.0;
        beta[3]=1.0/24.0;
    }
    else if(order==5)
    {
        beta[0]=251.0/720.0;
        beta[1]=646.0/720.0;
        beta[2]=-264.0/720.0;
        beta[3]=106.0/720.0;
        beta[4]=-19.0/720.0;
    }
    else if(order==6)
    {
        beta[0]=475.0/1440.0;
        beta[1]=1427.0/1440.0;
        beta[2]=-798.0/1440.0;
        beta[3]=482.0/1440.0;
        beta[4]=-173.0/1440.0;
        beta[5]=27.0/1440.0;
    }
    
    
    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n].dim;
    tot_dim=atoms->natms*c_dim;
    
    CREATE1D(a,tot_dim);
    CREATE1D(g,tot_dim);
    CREATE1D(g0,tot_dim);
    CREATE1D(c0,tot_dim);
    CREATE1D(h,tot_dim);
    CREATE1D(y,order-1);
    for(int i=0;i<order-1;i++)
        CREATE1D(y[i],tot_dim);
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_Adams::~Clock_Adams()
{
    delete [] t;
    delete [] beta;
    delete [] beta_mod;
    
    
    delete thermo;
    
    int natms=atoms->natms;
    for(int i=0;i<order-1;i++)
        if(natms)
            delete [] y[i];
    if(order>1)
        delete [] y;
    
    if(natms)
    {
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
void Clock_Adams::init()
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
    
    ave_err=0.0;
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_Adams::fin()
{
    ave_err=ave_err/static_cast<TYPE0>(no_steps);
    
    if(write!=NULL)
        write->fin();
    
    forcefield->fin();
    neighbor->fin();
    
    delete vecs_comm;
    thermo->fin();
    if(atoms->my_p_no==0)
        printf("ave. err.: %5.4e\n",ave_err);
    
    atoms->x2s(atoms->natms);
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_Adams::run()
{
    TYPE0 curr_err;
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    TYPE0* tmp_y;
    
    for(int istep=0;istep<no_steps;istep++)
    {
        if(istep<order)
        {
            for(int i=0;i<tot_dim;i++)
                a[i]=-c[i];
            curr_err=solve(1.0);
        }
        else
        {
            for(int i=0;i<tot_dim;i++)
            {
                a[i]=-c[i];
                for(int j=0;j<order-1;j++)
                    a[i]-=delta_t*beta[j+1]*y[j][i];
            }
            curr_err=solve(beta[0]);
        }
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step())
        {
            forcefield->force_calc(1,energy_stress);
            thermo->update(fe_idx,energy_stress[0]);
            thermo->update(stress_idx,6,&energy_stress[1]);
            thermo->update(time_idx,curr_err);
        }
        
        if(order>1)
        {
            tmp_y=y[order-2];
            for(int i=order-2;i>0;i--)
                y[i]=y[i-1];
            y[0]=tmp_y;
            memcpy(y[0],c_d,tot_dim*sizeof(TYPE0));
        }
        
        step_no++;
    }
    
    delete [] energy_stress;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
TYPE0 Clock_Adams::solve(TYPE0 bet)
{
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    TYPE0 gamma,max_gamma=1.0,min_gamma;
    TYPE0 inner,tmp;
    TYPE0 ratio;
    TYPE0 g0_g0,g_g,g_g0,g_h;
    TYPE0 curr_cost,ideal_cost,cost;
    
    min_gamma=1.0e-16;
    int chk;
    
    curr_cost=forcefield->calc_g(0,delta_t*bet,a,g);
    
    TYPE0 tot_ratio;
    
    forcefield->c_d_calc();
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
    for(int i=0;i<tot_dim;i++)
        c[i]+=c_d[i]*delta_t*tot_ratio;
    atoms->update(c_n);
    
    memcpy(h,g,tot_dim*sizeof(TYPE0));
    
    inner=0.0;
    for(int i=0;i<tot_dim;i++)
        inner+=g[i]*g[i];
    g0_g0=0.0;
    MPI_Allreduce(&inner,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    g_h=g0_g0;
    
    int iter=0;
    
    while(curr_cost>1.0e-8 && iter<MAX_ITER && max_gamma>min_gamma)
    {
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
            
            curr_cost=forcefield->calc_g(1,delta_t*bet,a,g);
            ideal_cost=cost-slope*max_gamma*g_h;
            if(curr_cost<ideal_cost)
                chk=0;
            max_gamma*=gamma_red;
            if(max_gamma<=0.0)
            {
                memcpy(c,c0,tot_dim*sizeof(TYPE0));
                atoms->update(c_n);
            }
        }
        
        curr_cost=forcefield->calc_g(0,delta_t*bet,a,g);
        
        
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
        iter++;
        
    }
    
    ave_err+=curr_cost;
    return curr_cost;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_Adams::interpolate(TYPE0 del_t)
{
    TYPE0 tmp0=1.0,tmp1;
    
    for(int i=2;i<order;i++)
    {
        beta_mod[i]=0.0;
        for(int j=2;j<order;j++)
        {
            tmp0=beta[j];
            tmp1=tot_time-(j-1.0)*del_t;
            for(int k=2;k<order;k++)
                if(k!=i)
                    tmp0*=(tmp1-t[k]);
            
            beta_mod[i]+=tmp0;
        }
        
        tmp0=1.0;
        
        for(int j=2;j<order;j++)
            if(i!=j)
                tmp0*=t[i]-t[j];
        
        beta_mod[i]/=tmp0;
    }
    
    beta_mod[1]=beta[1];
    beta_mod[0]=beta[0];
    
    for(int i=0;i<tot_dim;i++)
    {
        for(int j=0;j<order-1;j++)
            a[i]-=del_t*beta_mod[j+1]*y[j][i];
    }
}

