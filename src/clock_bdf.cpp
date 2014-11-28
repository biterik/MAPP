#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
#include "clock_bdf.h"
#define SMALL 1.0e-18
#define MAX_ITER 50
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_BDF::Clock_BDF(MAPP* mapp,int narg
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
    
    CREATE1D(alpha,order);
    if(order==1)
    {
        beta=1.0;
        alpha[0]=-1.0;
    }
    else if(order==2)
    {
        beta=2.0/3.0;
        alpha[0]=-4.0/3.0;
        alpha[1]=1.0/3.0;
    }
    else if(order==3)
    {
        beta=6.0/11.0;
        alpha[0]=-18.0/11.0;
        alpha[1]=9.0/11.0;
        alpha[2]=-2.0/11.0;
    }
    else if(order==4)
    {
        beta=12.0/25.0;
        alpha[0]=-48.0/25.0;
        alpha[1]=36.0/25.0;
        alpha[2]=-16.0/25.0;
        alpha[3]=3.0/25.0;
    }
    else if(order==5)
    {
        beta=60.0/137.0;
        alpha[0]=-300.0/137.0;
        alpha[1]=300.0/137.0;
        alpha[2]=-200.0/137.0;
        alpha[3]=75.0/137.0;
        alpha[4]=-12.0/137.0;
    }
    else if(order==6)
    {
        beta=60.0/147.0;
        alpha[0]=-360.0/147.0;
        alpha[1]=450.0/147.0;
        alpha[2]=-400.0/147.0;
        alpha[3]=225.0/147.0;
        alpha[4]=-72.0/147.0;
        alpha[5]=10.0/147.0;
    }

    
    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n].dim;
    tot_dim=atoms->natms*c_dim;

    CREATE1D(y_list,order);
    
    CREATE1D(a,tot_dim);
    CREATE1D(g,tot_dim);
    CREATE1D(g0,tot_dim);
    CREATE1D(c0,tot_dim);
    CREATE1D(h,tot_dim);
    CREATE1D(y,order);
    for(int i=0;i<order;i++)
        CREATE1D(y[i],tot_dim);
    
    
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_BDF::~Clock_BDF()
{
    delete [] alpha;
    delete [] y_list;
    
    delete thermo;
    
    int natms=atoms->natms;
    for(int i=0;i<order;i++)
        if(natms)
            delete [] y[i];
    
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
void Clock_BDF::init()
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
    
    thermo->update(pe_idx,energy_stress[0]);
    thermo->update(stress_idx,6,&energy_stress[1]);
    
    
    if(write!=NULL)
        write->init();
    thermo->init();
    
    delete [] energy_stress;
    
    for(int i=0;i<order;i++)
        y_list[i]=i;

    atoms->vectors[c_n].ret(c);
    
    for(int i=0;i<tot_dim;i++)
        y[y_list[0]][i]=c[i];
    
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_BDF::fin()
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
void Clock_BDF::run()
{
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
  
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    int no;
    
    for(int istep=0;istep<no_steps;istep++)
    {

        
        if(istep<order)
        {
            for(int i=0;i<tot_dim;i++)
                a[i]=-y[y_list[0]][i];
            
            solve(1.0);
        }
        else
        {
            for(int i=0;i<tot_dim;i++)
            {
                a[i]=0;
                for(int j=0;j<order;j++)
                    a[i]+=alpha[j]*y[y_list[j]][i];
            }
            solve(beta);
        }
        
        
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step())
        {
            forcefield->force_calc(1,energy_stress);
            thermo->update(pe_idx,energy_stress[0]);
            thermo->update(stress_idx,6,&energy_stress[1]);
        }
        

        no=y_list[order-1];
        for(int i=order-1;i>0;i--)
            y_list[i]=y_list[i-1];
        y_list[0]=no;
        memcpy(y[y_list[0]],c,tot_dim*sizeof(TYPE0));
        
        step_no++;
    }
    
    forcefield->force_calc(1,energy_stress);
    thermo->update(pe_idx,energy_stress[0]);
    thermo->update(stress_idx,6,&energy_stress[1]);
    delete [] energy_stress;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_BDF::solve(TYPE0 bet)
{
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    
    TYPE0 gamma,max_gamma=1.0,min_gamma;
    TYPE0 inner,tmp;
    TYPE0 ratio;
    TYPE0 g0_g0,g_g,g_g0,g_h;
    TYPE0 curr_cost,ideal_cost,cost;
//  TYPE0 upper=1.0-SMALL;
//  TYPE0 lower=SMALL;
    
    
    
    
    
    
    

    min_gamma=1.0e-16;
    int chk;

    /*
    for(int i=0;i<tot_dim;i++)
    {
        if(c[i]==0.0)
            c[i]+=0.001;
        else if (c[i]==1.0)
            c[i]-=0.001;
    }
     */
    atoms->update(c_n);

    
    curr_cost=forcefield->calc_g(0,delta_t*bet,a,g);
    
     /*
    forcefield->calc_y();
   

    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);

    
    for(int i=0;i<tot_dim;i++)
        printf("y[%d]: %e0\n",i,c_d[i]);
    for(int i=0;i<tot_dim;i++)
        printf("g[%d]: %e\n",i,g[i]);
    
    */
    
    memcpy(h,g,tot_dim*sizeof(TYPE0));
    
    inner=0.0;
    for(int i=0;i<tot_dim;i++)
        inner+=g[i]*g[i];
    g0_g0=0.0;
    MPI_Allreduce(&inner,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    g_h=g0_g0;
    
    int iter=0;
    
    while(curr_cost>SMALL && iter<MAX_ITER && max_gamma>min_gamma)
    {
        

        memcpy(g0,g,tot_dim*sizeof(TYPE0));
        memcpy(c0,c,tot_dim*sizeof(TYPE0));
        

        /*
        for(int i=0;i<tot_dim;i++)
        {
            tmp=c0[i]+h[i];
            if(c[i]==1.0)
            {
                if(h[i]>0.0)
                    h[i]=0.0;
            }
            else if(c[i]==0.0)
            {
                if(h[i]<0.0)
                    h[i]=0.0;
            }
        }
         */

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
    
        //max_gamma*=0.1;
        /*
        inner=0.0;
        for(int i=0;i<tot_dim;i++)
            inner+=h[i];
        printf("innnnnn  %e \n",inner);
          */  
        
        cost=curr_cost;
        //printf("beg  %e \n",g_h);
        while(chk && max_gamma>min_gamma)
        {
            for(int i=0;i<tot_dim;i++)
                c[i]=c0[i]+max_gamma*h[i];
            
            atoms->update(c_n);
            
            curr_cost=forcefield->calc_g(1,delta_t*bet,a,g);
            ideal_cost=cost-slope*max_gamma*g_h;
            //ideal_cost=cost;
            if(curr_cost<ideal_cost)
                chk=0;
            //printf("dddd gamma %e %e\n",max_gamma,curr_cost);
            max_gamma*=gamma_red;
            if(max_gamma<=0.0)
            {
                memcpy(c,c0,tot_dim*sizeof(TYPE0));
                atoms->update(c_n);
            }
        }
        
        curr_cost=forcefield->calc_g(0,delta_t*bet,a,g);
        //printf("end %e %e\n",curr_cost,max_gamma);
        
        
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
        
        //printf("g_h %lf \n",g_h);
        
        if(g_h<0.0)
        {
            
            memcpy(h,g,tot_dim*sizeof(TYPE0));
            g_h=g_g;
        }
        iter++;
        
    }
    //printf("cost: %e\n",curr_cost);
    
    
}

