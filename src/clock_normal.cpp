#include "clock_normal.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
#define SMALL 1.0e-18
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_NORMAL::Clock_NORMAL(MAPP* mapp,int narg
,char** arg):Clock(mapp)
{
    if(narg!=4)
        error->abort("invalid arguments");
    
    
    no_steps=atoi(arg[2]);
    delta_t=atof(arg[3]);
    
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_NORMAL::~Clock_NORMAL()
{
    delete thermo;

}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_NORMAL::init()
{
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<TYPE0>(0,atoms->vectors[0].dim,"f");
    
    c_n=atoms->find_exist("c");
    c_d_n=atoms->find("c_d");
    
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

}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_NORMAL::fin()
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
void Clock_NORMAL::run()
{
    
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    
    
    TYPE0 tmp,curr_t;
    TYPE0 ratio,tot_ratio;
    TYPE0 upper=1.0-SMALL;
    TYPE0 lower=SMALL;
    int natms=atoms->natms;
    int tot_dim=atoms->vectors[c_n].dim*natms;
    
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    TYPE0 tot_t;
    
    for(int istep=0;istep<no_steps;istep++)
    {
        tot_t=delta_t;
        while(tot_t>0.0)
        {
            forcefield->calc_y();
            ratio=1.0;
            
            for(int i=0;i<tot_dim;i++)
            {
                tmp=c_d[i]*tot_t+c[i];
                if(tmp>upper)
                    ratio=MIN((upper-c[i])/(tmp-c[i]),ratio);
                else if(tmp<lower)
                    ratio=MIN((c[i]-lower)/(c[i]-tmp),ratio);
            }
            MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
            
            curr_t=tot_ratio*tot_t;
            
            for(int i=0;i<tot_dim;i++)
                c[i]+=c_d[i]*curr_t;
            atoms->update(c_n);
            tot_t-=curr_t;
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
        

        step_no++;

    }
    
    forcefield->force_calc(1,energy_stress);
    thermo->update(pe_idx,energy_stress[0]);
    thermo->update(stress_idx,6,&energy_stress[1]);
    delete [] energy_stress;

}
