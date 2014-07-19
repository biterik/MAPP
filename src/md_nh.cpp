/*--------------------------------------------
 Created by Sina on 06/20/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <cmath>
#include <stdlib.h>
#include "md_nh.h"
#include "memory.h"
#include "ff.h"
#include "neighbor.h"
#include "error.h"
#include "rand_engine.h"
#include "atoms.h"
#include "atom_types.h"
#include "thermo_dynamics.h"
#include "write.h"
#define TAU 1
#define XYZ 2
#define YZ 3
#define ZX 4
#define XY 5
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MD_NH::MD_NH(MAPP* mapp,int narg,char** arg)
: MD(mapp)
{
    if(forcefield==NULL)
        error->abort("force field is not initiated");
    if(atoms->dimension!=3)
        error->abort("to use md, dimension of the box should be 3");
    
    
    //the defaults
    no_it_eta=1;
    no_ch_eta=3;
    t_tar=0;
    chk_create_vel=0;
    dt=0.0;

    
    int iarg=2;
    if(!strcmp(arg[iarg],"ntv"))
    {
        chk_stress=0;
        iarg++;
    }
    else if(!strcmp(arg[iarg],"nt\tau"))
    {
        chk_stress=TAU;
        iarg++;
    }
    else if(!strcmp(arg[iarg],"ntp")&&((iarg+2)<=narg))
    {
        
        iarg++;
        if(!strcmp(arg[iarg],"iso"))
        {
            chk_stress=XYZ;
            iarg++;
        }
        else if(!strcmp(arg[iarg],"yz"))
        {
            chk_stress=YZ;
            iarg++;
        }
        else if(!strcmp(arg[iarg],"zx"))
        {
            chk_stress=ZX;
            iarg++;
        }
        else if(!strcmp(arg[iarg],"xy"))
        {
            chk_stress=XY;
            iarg++;
        }
        else
            error->abort("wrong run command"
                " coupling style: %s",arg[iarg]);
        
    }
    else error->abort("wrong run command style %s",arg[iarg]);
  
    

    // create and allocate the memory for
    // necessary

    if(chk_stress)
    {
        CREATE2D(M1,3,3);
        CREATE2D(M2,3,3);
        CREATE1D(chk_tau,6);
        for(int i=0;i<6;i++)
            chk_tau[i]=0;
        CREATE1D(v_per_atm,6);
        CREATE1D(tau_freq,6);
        no_it_peta=1;
        no_ch_peta=3;
        
    }

    CREATE1D(ke_curr,6);

    
    
    while(iarg<narg)
    {
        if(!strcmp(arg[iarg],"temp")&&((iarg+2)<=narg))
        {
            iarg++;
            t_tar=atof(arg[iarg]);
            if(t_tar<=0)
                error->abort("temperature cannot be equal or less than zero");
            iarg++;
            if(atof(arg[iarg])<=0)
                error->abort("temperature period cannot be equal or less than zero");
            t_freq=1.0/atof(arg[iarg]);
            iarg++;
        }
        else if(!strcmp(arg[iarg],"stress")&&((iarg+7)<=narg))
        {
            if(chk_stress!=TAU)
                error->abort("wrong command");
            iarg++;
            for (int i=0;i<6;i++)
            {
                tau_tar[i]=atof(arg[iarg]);
                iarg++;
                
                if (atof(arg[iarg])<=0)
                    error->abort("tau period cannot be or smaller than zero");
                tau_freq[i]=1.0/atof(arg[iarg]);

                chk_tau[i]=1;
                iarg++;
            }
        }
        else if(!strcmp(arg[iarg],"ave")&&((iarg+2)<=narg))
        {
            if(chk_stress<=TAU)
                error->abort("wrong command");
            iarg++;
            if(chk_stress==XYZ)
            {
                tau_tar[0]=tau_tar[1]=tau_tar[2]=atof(arg[iarg]);
                chk_tau[0]=chk_tau[1]=chk_tau[2]=1;
                iarg++;
                if (atof(arg[iarg])<=0)
                    error->abort("tau period cannot be or smaller than zero");
                tau_freq[0]=tau_freq[1]=tau_freq[2]=1.0/atof(arg[iarg]);
                iarg++;
            }
            else if(chk_stress==YZ)
            {
                tau_tar[1]=tau_tar[2]=atof(arg[iarg]);
                chk_tau[1]=chk_tau[2]=1;
                iarg++;
                if (atof(arg[iarg])<=0)
                    error->abort("tau period cannot be or smaller than zero");
                tau_freq[1]=tau_freq[2]=1.0/atof(arg[iarg]);
                iarg++;
            }
            else if(chk_stress==ZX)
            {
                tau_tar[0]=tau_tar[2]=atof(arg[iarg]);
                chk_tau[0]=chk_tau[2]=1;
                iarg++;
                if (atof(arg[iarg])<=0)
                    error->abort("tau period cannot be or smaller than zero");
                tau_freq[0]=tau_freq[2]=1.0/atof(arg[iarg]);
                iarg++;
            }
            else if(chk_stress==XY)
            {
                tau_tar[0]=tau_tar[1]=atof(arg[iarg]);
                chk_tau[0]=chk_tau[1]=1;
                iarg++;
                if (atof(arg[iarg])<=0)
                    error->abort("tau period cannot be or smaller than zero");
                tau_freq[0]=tau_freq[1]=1.0/atof(arg[iarg]);
                iarg++;
            }
            else
                error->abort("wrong ave style in run command");
                
        }
        else if(!strcmp(arg[iarg],"eta_iteration")&&((iarg+2)<=narg))
        {
            iarg++;
            no_it_eta=atoi(arg[iarg]);
            if (no_it_eta<1)
                error->abort("no_it_eta cannot be smaller than 1");
            iarg++;
        }
        else if(!strcmp(arg[iarg],"peta_iteration")&&((iarg+2)<=narg))
        {
            if(chk_stress==0)
                error->abort("wrong command");
            iarg++;
            no_it_peta=atoi(arg[iarg]);
            if (no_it_peta<1)
                error->abort("no_it_peta cannot be smaller than 1");
            iarg++;
        }
        else if(!strcmp(arg[iarg],"eta_chains")&&((iarg+2)<=narg))
        {
            iarg++;
            no_ch_eta=atoi(arg[iarg]);
            if (no_ch_eta<3)
                error->abort("no_ch_eta cannot be smaller than 3");
            iarg++;
        }
        else if(!strcmp(arg[iarg],"peta_chains")&&((iarg+2)<=narg))
        {
            if(chk_stress==0)
                error->abort("wrong command");
            iarg++;
            no_ch_peta=atoi(arg[iarg]);
            if (no_ch_peta<3)
                error->abort("no_ch_peta cannot be smaller than 3");
            iarg++;
        }
        else if(!strcmp(arg[iarg],"create_vel")&&((iarg+2)<=narg))
        {
            chk_create_vel=1;
            iarg++;
            seed=atoi(arg[iarg]);
            iarg++;
        }
        else error->abort("wrong command: %s",arg[iarg]);
    }
    
    // check wether all the relavant quantities are assigned or not
    // check wether the temperature is set
    if(t_tar==0)
        error->abort("temperature is not set yet");
    // check wether the stresses are set
    if(chk_stress)
    {
        int temp=1;
        if(chk_stress==TAU)
        {
            for(int i=0;i<6;i++)
                temp*=chk_tau[i];
        }
        else if(chk_stress==XYZ)
        {
            temp=chk_tau[0]*chk_tau[1]*chk_tau[2];
            for(int i=0;i<6;i++)
                if(i!=0&&i!=1&&i!=2)
                    temp*=1-chk_tau[i];
            
        }
        else if(chk_stress==YZ)
        {
            temp=chk_tau[1]*chk_tau[2];
            for(int i=0;i<6;i++)
                if(i!=1&&i!=2)
                    temp*=1-chk_tau[i];
            
        }
        else if(chk_stress==ZX)
        {
            temp=chk_tau[0]*chk_tau[2];
            for(int i=0;i<6;i++)
                if(i!=0&&i!=2)
                    temp*=1-chk_tau[i];
            
        }
        else if(chk_stress==XY)
        {
            temp=chk_tau[0]*chk_tau[1];
            for(int i=0;i<6;i++)
                if(i!=0&&i!=1)
                    temp*=1-chk_tau[i];
            
        }
        else
            error->abort("wrong stress style");
        
        if(temp==0)
            error->abort("not all the stresses are set");
    }
        
    // allocate the Nose Hoover chains stuff
    
    if(chk_stress)
    {
        CREATE1D(peta_d,no_ch_peta);
        CREATE1D(peta_dd,no_ch_peta);
        CREATE1D(peta_m,no_ch_peta);
        for(int i=0;i<no_ch_peta;i++)
            peta_d[i]=0.0;
    }
    
    CREATE1D(eta_d,no_ch_eta);
    CREATE1D(eta_dd,no_ch_eta);
    CREATE1D(eta_m,no_ch_eta);
    for(int i=0;i<no_ch_eta;i++)
        eta_d[i]=0.0;

    for (int i=0;i<6;i++)
        ke_curr[i]=0.0;
    
    tau_freq_m=0.0;
    if(chk_stress)
        for(int i=0;i<6;i++)
            if(chk_tau[i])
                tau_freq_m=MAX(tau_freq_m,tau_freq[i]);
    
    CREATE1D(x_ave,3);
    CREATE1D(x_ave_tot,3);
    
    boltz=0.0;
    
    
    char** argss;
    int nargs=mapp->parse_line((char*)"KE Temp. PE S_xx S_yy S_zz S_yz S_zx S_xy",argss);
    ke_idx=0;
    temp_idx=1;
    pe_idx=2;
    stress_idx=3;
    
    thermo=new ThermoDynamics(mapp,nargs,argss);
    for(int i=0;i<nargs;i++)
        delete [] argss[i];
    delete [] argss;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
MD_NH::~MD_NH()
{
    delete thermo;
    
    if(chk_stress)
    {
        for(int i=0;i<3;i++)
        {
            delete [] M1[i];
            delete [] M2[i];
        }
        delete [] M1;
        delete [] M2;
        delete [] ke_curr;
        delete [] chk_tau;
        delete [] v_per_atm;
        delete [] tau_tar;
        delete [] tau_freq;
    }
    
    delete [] eta_d;
    delete [] eta_dd;
    delete [] eta_m;
    
    if(chk_stress)
    {
        delete [] peta_d;
        delete [] peta_dd;
        delete [] peta_m;
    }
    
    delete [] x_ave_tot;
    delete [] x_ave;
    
}
/*--------------------------------------------
 setup before start of run
 --------------------------------------------*/
void MD_NH::init()
{
    if(boltz==0.0)
        error->abort("boltzmann constant is not set");
    if(dt==0.0)
        error->abort("time step is not set");
    
    no_dof=atoms->tot_natms*3-3;
    ke_tar=t_tar*boltz*no_dof;
    dt2=0.5*dt;
    dt4=0.25*dt;
    dt8=0.125*dt;
    x_n=atoms->find("x");
    type_n=atoms->find("type");
    
    x_d_n=atoms->find_exist("x_d");
    if(x_d_n<0)
    {
        x_d_n=atoms->add<TYPE0>(0, 3,"x_d");
        if(chk_create_vel==0)
            error->abort("since no velocity exists you should include create_vel command");
    }
    
    f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<TYPE0>(0, 3,"f");
    
    /*
     0. set the atomic vectors communication
     we need force and type in addition to
     x vector; however, we do not need to
     include force since their pervious
     values are not important and are set to
     0; when a communication happens atoms
     allocates the f vector with the
     appropriate size
     */
    vecs_comm=new VecLst(mapp,3,x_n,x_d_n,type_n);
    vecs_comm->add_update(x_n);
    atoms->reset_comm(vecs_comm);
    /*--------------------------------------*/
    
    /*
     1. setup the cutoffs for neighbor lists
     and also the maximum for communication
     in atoms class (this is needed for
     finding the phantom atoms and updating
     them);
     */
    forcefield->init();
    /*--------------------------------------*/
    
    /*
     2. find the phantom atoms and
     communicate them between the processors
     using the values assigned by the
     previous step;
     */
    atoms->ph_setup(1,vecs_comm);
    /*--------------------------------------*/
    
    /*
     3. creates the neighbor list for the
     bins;
     */
    neighbor->init();
    neighbor->create_list(0,1);

    /*--------------------------------------*/


    if(chk_create_vel)
        create_vel(seed,t_tar);
    else
        init_vel(t_tar);
    
    zero_f();
    
    
    
    
    TYPE0* enst;
    CREATE1D(enst,7);
    forcefield->force_calc(1,enst);

    
    thermo->update(stress_idx,6,&enst[1]);
    thermo->update(pe_idx,enst[0]);
    thermo->update(ke_idx,ke_cur);
    thermo->update(temp_idx,t_cur);
    delete [] enst;
    

    
    for (int i=0;i<no_ch_eta;i++)
        eta_m[i] = boltz * t_tar/(t_freq*t_freq);
    eta_m[0]*=no_dof;

    for (int i=0;i<no_ch_eta;i++)
        eta_d[i]=0.0;
    
    for (int i=1;i<no_ch_eta;i++)
        eta_dd[i] = (eta_m[i-1]*eta_d[i-1]*eta_d[i-1]-
                    boltz*t_tar)/eta_m[i];
    
    if(chk_stress)
    {
       
        for (int i=0;i<no_ch_peta;i++)
            peta_m[i] = boltz * t_tar/(tau_freq_m*tau_freq_m);
        
        
        for (int i=0;i<no_ch_peta;i++)
            peta_d[i]=0.0;
        
        for (int i=1;i<no_ch_peta;i++)
            peta_dd[i] = (peta_m[i-1]*peta_d[i-1]*peta_d[i-1]-
                         boltz*t_tar)/peta_m[i];
         
    }
    
    
    thermo->init();
    if(write!=NULL)
        write->init();
}
/*--------------------------------------------
 finalize after the run is complete
 --------------------------------------------*/
void MD_NH::fin()
{
    thermo->fin();
    if(write!=NULL)
        write->fin();
    if(atoms->my_p_no==0)
        fprintf(output,"\n");
    
    forcefield->fin();
    neighbor->fin();
    delete vecs_comm;
    atoms->x2s(atoms->natms);
}
/*--------------------------------------------
 MDrun
 --------------------------------------------*/
void MD_NH::run(int no_stps)
{
    TYPE0* enst;
    CREATE1D(enst,7);
    if(chk_stress)
    {
        for(int i=0;i<no_stps;i++)
        {
            update_NH_tau(dt2);
            update_NH_T(dt2);
            update_omega_d(dt2);
            update_x_d_xpnd(dt2);
            update_x_d(dt2);
            update_x(dt);
            
            zero_f();
            thermo->thermo_print();
            if(write!=NULL)
                write->write();
            
            forcefield->force_calc(1,enst);
            
            if(thermo->test_prev_step()|| i==no_stps-1)
            {
                thermo->update(stress_idx,6,&enst[1]);
                thermo->update(pe_idx,enst[0]);
                thermo->update(ke_idx,ke_cur);
                thermo->update(temp_idx,t_cur);
                
            }
            
            for (int j=0;j<6;j++)
                v_per_atm[j]=enst[1+j];
            
            update_x_d(dt2);
            update_x_d_xpnd(dt2);
            update_omega_d(dt2);
            update_NH_T(dt2);
            update_NH_tau(dt2);
            step_no++;
        }
    }
    else
    {

        for(int i=0;i<no_stps;i++)
        {
            update_NH_T(dt2);
            update_x_d(dt2);
            update_x(dt);
            zero_f();
            thermo->thermo_print();
            if(write!=NULL)
                write->write();
            
            if(thermo->test_prev_step()|| i==no_stps-1)
            {
                forcefield->force_calc(1,&enst[0]);
                thermo->update(stress_idx,6,&enst[1]);
                thermo->update(pe_idx,enst[0]);
                thermo->update(ke_idx,ke_cur);
                thermo->update(temp_idx,t_cur);
                
            }
            else
            {
                forcefield->force_calc(0,&enst[0]);
            }

            update_x_d(dt2);
            update_NH_T(dt2);
            step_no++;
        }

    }
    
    delete [] enst;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_NH::update_H(TYPE0 dlt)
{
    TYPE0 dlt2,dlt4,dlt8,exfac;
    dlt2=0.5*dlt;
    dlt4=0.25*dlt;
    dlt8=0.125*dlt;
    
    TYPE0 H0[3][3];
    TYPE0 H0_inv[3][3];
    
    M3EQV(atoms->H,H0);
    
    for (int i=0;i<2;i++)
    {
        M3INV_TRI_LOWER(H0,H0_inv);
        
        if (chk_tau[4])
        {
            exfac=exp(dlt8*omega_d[0]);
            H0[2][0]*=exfac;
            H0[2][0]+=dlt4*(omega_d[5]*H0[2][1]
                            +omega_d[4]*H0[2][2]);
            H0[2][0]*=exfac;
        }
        
        if (chk_tau[3])
        {
            exfac=exp(dlt4*omega_d[1]);
            H0[2][1]*=exfac;
            H0[2][1]+=dlt2*omega_d[3]*H0[2][2];
            H0[2][1]*=exfac;
        }
        
        if (chk_tau[5])
        {
            exfac=exp(dlt4*omega_d[0]);
            H0[1][0]*=exfac;
            H0[1][0]+=dlt2*omega_d[5]*H0[1][1];
            H0[1][0]*=exfac;
        }
        
        if (chk_tau[4])
        {
            exfac=exp(dlt8*omega_d[0]);
            H0[2][0]*=exfac;
            H0[2][0]+=dlt4*(omega_d[5]*H0[2][1]+omega_d[4]*H0[2][2]);
            H0[2][0]*=exfac;
        }
        
        if (chk_tau[0])
        {
            exfac=exp(dlt*omega_d[0]);
            H0[0][0]*=exfac;
        }
        
        if (chk_tau[1])
        {
            exfac=exp(dlt*omega_d[1]);
            H0[1][1]*=exfac;
            H0[1][0]*=exfac;
        }
        
        if (chk_tau[2])
        {
            exfac=exp(dlt*omega_d[2]);
            H0[2][2]*=exfac;
            H0[2][1]*=exfac;
            H0[2][0]*=exfac;
        }
        
        if (chk_tau[4])
        {
            exfac=exp(dlt8*omega_d[0]);
            H0[2][0]*=exfac;
            H0[2][0]+=dlt4*(omega_d[5]*H0[2][1]+omega_d[4]*H0[2][2]);
            H0[2][0]*=exfac;
        }
        
        if (chk_tau[5])
        {
            exfac=exp(dlt4*omega_d[0]);
            H0[1][0]*=exfac;
            H0[1][0]+=dlt2*omega_d[5]*H0[1][1];
            H0[1][0]*=exfac;
        }
        
        if (chk_tau[3])
        {
            exfac=exp(dlt4*omega_d[1]);
            H0[2][1]*=exfac;
            H0[2][1]+=dlt2*omega_d[3]*H0[2][2];
            H0[2][1]*=exfac;
        }
        
        if (chk_tau[4])
        {
            exfac=exp(dlt8*omega_d[0]);
            H0[2][0]*=exfac;
            H0[2][0]+=dlt4*(omega_d[5]*H0[2][1]+omega_d[4]*H0[2][2]);
            H0[2][0]*=exfac;
        }
        
        if (i==0) M3MUL_TRI_LOWER(H0_inv,H0,M1);
        else if (i==1) M3MUL_TRI_LOWER(H0_inv,H0,M2);
    }
    
    M3EQV(H0,atoms->H);
    M3INV_TRI_LOWER(atoms->H,atoms->B);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_NH::update_x(TYPE0 dlt)
{
    //TYPE0* x=(TYPE0*)atoms->vectors[x_n].ret_vec();
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    
    //TYPE0* x_d=(TYPE0*)atoms->vectors[x_d_n].ret_vec();
    TYPE0* x_d;
    atoms->vectors[x_d_n].ret(x_d);
    
    int natms=atoms->natms;
    int icomp;
    x_ave[0]=x_ave[1]=x_ave[2]=0.0;
    x_ave_tot[0]=x_ave_tot[1]=x_ave_tot[2]=0.0;
    
    if (chk_stress)
    {
        update_H(0.5*dlt);
        TYPE0 xt[3];
        for(int i=0;i<natms;i++)
        {
            icomp=3*i;
            x_ave[0]-=x[icomp];
            x_ave[1]-=x[icomp+1];
            x_ave[2]-=x[icomp+2];
            
            xt[0]=x[icomp]*M1[0][0]
            +x[icomp+1]*M1[1][0]
            +x[icomp+2]*M1[2][0];
            xt[1]=x[icomp+1]*M1[1][1]
            +x[icomp+2]*M1[2][1];
            xt[2]=x[icomp+2]*M1[2][2];
            
            xt[0]+=x_d[icomp]*dlt;
            xt[1]+=x_d[icomp+1]*dlt;
            xt[2]+=x_d[icomp+2]*dlt;
            
            x[icomp]=xt[0]*M2[0][0]
            +xt[1]*M2[1][0]
            +xt[2]*M2[2][0];
            x[icomp+1]=xt[1]*M2[1][1]
            +xt[2]*M2[2][1];
            x[icomp+2]=xt[2]*M2[2][2];
            
            x_ave[0]+=x[icomp];
            x_ave[1]+=x[icomp+1];
            x_ave[2]+=x[icomp+2];
        }
    }
    else
    {
        for(int i=0;i<natms;i++)
        {
            icomp=3*i;
            x[icomp]+=x_d[icomp]*dlt;
            x[icomp+1]+=x_d[icomp+1]*dlt;
            x[icomp+2]+=x_d[icomp+2]*dlt;
            x_ave[0]+=x_d[icomp]*dlt;
            x_ave[1]+=x_d[icomp+1]*dlt;
            x_ave[2]+=x_d[icomp+2]*dlt;
        }
    }
    
    MPI_Allreduce(&x_ave[0],&x_ave_tot[0],3,MPI_TYPE0,MPI_SUM,world);
    x_ave_tot[0]*=1.0/(atoms->tot_natms);
    x_ave_tot[1]*=1.0/(atoms->tot_natms);
    x_ave_tot[2]*=1.0/(atoms->tot_natms);
    
    for(int i=0;i<natms;i++)
    {
        icomp=3*i;
        x[icomp]-=x_ave_tot[0];
        x[icomp+1]-=x_ave_tot[1];
        x[icomp+2]-=x_ave_tot[2];
    }

    if (chk_stress)
        atoms->update_0(1,1,vecs_comm);
    else
        atoms->update_0(0,1,vecs_comm);

}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_NH::update_x_d(TYPE0 dlt)
{
    

    /*
    TYPE0* x_d=(TYPE0*)atoms->vectors[x_d_n].ret_vec();
    TYPE0* f=(TYPE0*)atoms->vectors[f_n].ret_vec();
    int* type=(int*)atoms->vectors[type_n].ret_vec();
     */
    TYPE0* x_d;
    atoms->vectors[x_d_n].ret(x_d);
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    int* type;
    atoms->vectors[type_n].ret(type);
    
    
    TYPE0* mass=atom_types->mass;
    int natms=atoms->natms;
    TYPE0 temp[6];
    int icomp;
    
    for(int i=0;i<6;i++)
        temp[i]=0.0;
    
    for(int i=0;i<natms;i++)
    {
        icomp=3*i;
        
        x_d[icomp]+=f[icomp]*dlt/mass[type[i]];
        x_d[icomp+1]+=f[icomp+1]*dlt/mass[type[i]];
        x_d[icomp+2]+=f[icomp+2]*dlt/mass[type[i]];
        
        temp[0]+=mass[type[i]]*x_d[icomp]*x_d[icomp];
        temp[1]+=mass[type[i]]*x_d[icomp+1]*x_d[icomp+1];
        temp[2]+=mass[type[i]]*x_d[icomp+2]*x_d[icomp+2];
        temp[3]+=mass[type[i]]*x_d[icomp+1]*x_d[icomp+2];
        temp[4]+=mass[type[i]]*x_d[icomp]*x_d[icomp+2];
        temp[5]+=mass[type[i]]*x_d[icomp]*x_d[icomp+1];
        

    }
    
    for(int i=0;i<6;i++)
        ke_curr[i]=0.0;
    MPI_Allreduce(temp,ke_curr,6,MPI_TYPE0,MPI_SUM,world);
    ke_cur=(ke_curr[0]+ke_curr[1]+ke_curr[2]);

    t_cur=ke_cur/(boltz*3*atoms->tot_natms);

}
/*--------------------------------------------
 Nosé–Hoover thermostat chains 
 --------------------------------------------*/
void MD_NH::update_NH_T(TYPE0 dlt)
{
    TYPE0 dltm,dltm2,dltm4,exfac,velfac;
    int natoms=atoms->natms;

    
    for (int i=0;i<no_ch_eta;i++)
        eta_m[i]=boltz*t_tar/(t_freq*t_freq);
    eta_m[0]*=no_dof;
    
    dltm=dlt*(1.0/no_it_eta);
    dltm2=0.5*dltm;
    dltm4=0.25*dltm;
    velfac=1.0;
    
    eta_dd[0]=(ke_cur-ke_tar)/eta_m[0];
    for (int it=0;it<no_it_eta;it++)
    {
        exfac=1.0;
        for (int ich=no_ch_eta-1;ich>-1;ich--)
        {
            eta_d[ich]*=exfac;
            eta_d[ich]+=eta_dd[ich]*dltm2;
            eta_d[ich]*=exfac;
            exfac=exp(-dltm4*eta_d[ich]);
        }

        //rescale x_d dlt & claculate the new temperature
        exfac=exp(-dltm*eta_d[0]);
        velfac*=exfac;
        t_cur*=exfac*exfac;
        ke_cur*=exfac*exfac;
        ke_curr[0]*=exfac*exfac;
        ke_curr[1]*=exfac*exfac;
        ke_curr[2]*=exfac*exfac;
        ke_curr[3]*=exfac*exfac;
        ke_curr[4]*=exfac*exfac;
        ke_curr[5]*=exfac*exfac;
        
        exfac=exp(-dltm4*eta_d[1]);
        eta_d[0]*=exfac;
        eta_dd[0]=(ke_cur-ke_tar)/eta_m[0];
        eta_d[0]+=eta_dd[0]*dltm2;
        eta_d[0]*=exfac;
        
        for (int ich=1;ich<no_ch_eta;ich++)
        {
            if (ich==no_ch_eta-1) exfac=1.0;
            else exfac=exp(-dltm4*eta_d[ich+1]);
            
            eta_d[ich]*=exfac;
            eta_dd[ich]=(eta_m[ich-1]*eta_d[ich-1]*eta_d[ich-1]-boltz*t_tar)/eta_m[ich];
            eta_d[ich]+=eta_dd[ich]*dltm2;
            eta_d[ich]*=exfac;
        }
    }
    
    
    //TYPE0* x_d=(TYPE0*)atoms->vectors[x_d_n].ret_vec();
    TYPE0* x_d;
    atoms->vectors[x_d_n].ret(x_d);
    for (int i=0;i<3*natoms;i++)
        x_d[i]*=velfac;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_NH::update_NH_tau(TYPE0 dlt)
{
    TYPE0 dltm,dltm2,dltm4,exfac,kec;
    int dof=0;
    
    for (int i=0;i<6;i++)
        if (chk_tau[i]) dof++;
    
    for (int i=0;i<6;i++)
        if(chk_tau[i])
            omega_m[i]=boltz*t_tar/(tau_freq[i]*tau_freq[i]);
    
    for (int i=0;i<no_ch_peta;i++)
        peta_m[i]=boltz*t_tar/(tau_freq_m*tau_freq_m);
    
    for (int i=1;i<no_ch_peta;i++)
        peta_dd[i]=(peta_m[i-1]*peta_d[i-1]*peta_d[i-1]-boltz*t_tar)/peta_m[i];
    kec=0.0;
    for (int i=0;i<6;i++)
        kec+=omega_m[i]*omega_d[i]*omega_d[i];
    peta_dd[0]=(kec-boltz*t_tar)/peta_m[0];
    
    dltm=dlt*(1.0/no_it_peta);
    dltm2=0.5*dltm;
    dltm4=0.25*dltm;
    
    for (int it=0;it<no_ch_peta;it++)
    {
        exfac=1.0;
        for (int ich=no_ch_peta-1;ich>-1;ich--)
        {
            peta_d[ich]*=exfac;
            peta_d[ich]+=peta_dd[ich]*dltm2;
            peta_d[ich]*=exfac;
            exfac=exp(-dltm4*peta_d[ich]);
        }
        
        exfac=exp(-dltm*peta_d[0]);
        for (int i=0;i<6;i++)
            omega_d[i]*=exfac;
        
        kec*=exfac*exfac;
        
        exfac=exp(-dltm4*peta_d[1]);
        peta_d[0]*=exfac;
        peta_dd[0]=(kec-dof*boltz*t_tar)/peta_m[0];
        peta_d[0]+=peta_dd[0]*dltm2;
        peta_d[0]*=exfac;
        
        for (int ich=1;ich<no_ch_peta;ich++)
        {
            if (ich==no_ch_peta-1) exfac=1.0;
            else exfac=exp(-dltm4*peta_d[ich+1]);
            
            peta_d[ich]*=exfac;
            peta_dd[ich]=(peta_m[ich-1]*peta_d[ich-1]*peta_d[ich-1]-boltz*t_tar)/peta_m[ich];
            peta_d[ich]+=peta_dd[ich]*dltm2;
            peta_d[ich]*=exfac;
        }
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_NH::update_omega_d(TYPE0 dlt)
{
    TYPE0** H=atoms->H;
    TYPE0 vol= H[0][0]*H[1][1]*H[2][2];
    MTK_1=0.0;
    for(int i=0;i<3;i++)
        if (chk_tau[i])
            MTK_1+=ke_curr[i];
    MTK_1/=(chk_tau[0]+chk_tau[1]+chk_tau[2])
    *atoms->tot_natms;
    
    couple();
    for (int i=0;i<6;i++)
        if (chk_tau[i])
            omega_d[i]=((v_per_atm[i]+ke_curr[i]
                        -tau_tar[i]*vol)
                        -MTK_1)*dlt/omega_m[i];
    
    MTK_2=0.0;
    for (int i=0;i<3;i++)
        if (chk_tau[i])
            MTK_2+=omega_d[i];

    MTK_2/=(chk_tau[0]+chk_tau[1]+chk_tau[2])
    *atoms->tot_natms;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void MD_NH::update_x_d_xpnd(TYPE0 dlt)
{
    TYPE0 fac[3];
    TYPE0 temp[6];
    
    

    /*
    TYPE0* x_d=(TYPE0*)atoms->vectors[x_d_n].ret_vec();
    int* type=(int*)atoms->vectors[type_n].ret_vec();
     */
    TYPE0* x_d;
    atoms->vectors[x_d_n].ret(x_d);
    int* type;
    atoms->vectors[type_n].ret(type);
    TYPE0* mass=atom_types->mass;
    
    int natms=atoms->natms;
    int icomp;
    
    fac[0]=exp(0.5*dlt*(omega_d[0]+MTK_2));
    fac[1]=exp(0.5*dlt*(omega_d[1]+MTK_2));
    fac[2]=exp(0.5*dlt*(omega_d[2]+MTK_2));
    
    for(int i=0;i<6;i++)
        temp[i]=0.0;
    for (int i=0;i<natms;i++)
    {
        icomp=3*i;
        x_d[icomp]*=fac[0];
        x_d[icomp+1]*=fac[1];
        x_d[icomp+2]*=fac[2];
        
        x_d[icomp]-=dlt*(x_d[icomp+1]*omega_d[5]
                       +x_d[icomp+2]*omega_d[4]);
        x_d[icomp+1]-=dlt*x_d[icomp+2]*omega_d[3];
        
        x_d[icomp]*=fac[0];
        x_d[icomp+1]*=fac[1];
        x_d[icomp+2]*=fac[2];
        
        temp[0]+=mass[type[i]]*x_d[icomp]*x_d[icomp];
        temp[1]+=mass[type[i]]*x_d[icomp+1]*x_d[icomp+1];
        temp[2]+=mass[type[i]]*x_d[icomp+2]*x_d[icomp+2];
        temp[3]+=mass[type[i]]*x_d[icomp+1]*x_d[icomp+2];
        temp[4]+=mass[type[i]]*x_d[icomp]*x_d[icomp+2];
        temp[5]+=mass[type[i]]*x_d[icomp]*x_d[icomp+1];
    }
    
    for(int i=0;i<6;i++)
        ke_curr[i]=0.0;
    MPI_Allreduce(temp,ke_curr,6,MPI_TYPE0,MPI_SUM,world);
    ke_cur=(ke_curr[0]+ke_curr[1]+ke_curr[2]);
    t_cur=ke_cur/(boltz*3*atoms->tot_natms);
}
/*--------------------------------------------
 zero acceleration
 --------------------------------------------*/
void MD_NH::zero_f()
{
    //TYPE0* f=(TYPE0*)atoms->vectors[f_n].ret_vec();
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    for(int i=0;i<atoms->natms*3;i++)
        f[i]=0.0;
}
/*--------------------------------------------
 create initial velocity 
 --------------------------------------------*/
void MD_NH::create_vel(int seed,TYPE0 temperature)
{
    TYPE0* x_d;
    atoms->vectors[atoms->find("x_d")].ret(x_d);
    int* type;
    atoms->vectors[atoms->find("type")].ret(type);
    TYPE0* mass=atom_types->mass;
    
    int natms=atoms->natms;
    int icomp;
    TYPE0* temp;
    CREATE1D(temp,6);
    for(int i=0;i<6;i++)
        temp[i]=0.0;
    

    
    
    class Random* random=new Random(mapp,seed);
    for(int i=0;i<natms;i++)
    {
        icomp=3*i;
        for(int j=0;j<3;j++)
            x_d[icomp+j]=random->gaussian()/(sqrt(mass[type[i]]));

        temp[0]+=mass[type[i]]*x_d[icomp]*x_d[icomp];
        temp[1]+=mass[type[i]]*x_d[icomp+1]*x_d[icomp+1];
        temp[2]+=mass[type[i]]*x_d[icomp+2]*x_d[icomp+2];
        temp[3]+=mass[type[i]]*x_d[icomp+1]*x_d[icomp+2];
        temp[4]+=mass[type[i]]*x_d[icomp]*x_d[icomp+2];
        temp[5]+=mass[type[i]]*x_d[icomp]*x_d[icomp+1];
    }
    delete random;
    
    for(int i=0;i<6;i++)
        ke_curr[i]=0.0;

    
    MPI_Allreduce(&temp[0],&ke_curr[0],6,MPI_TYPE0,MPI_SUM,world);
    
    
    ke_cur=(ke_curr[0]+ke_curr[1]+ke_curr[2]);
    t_cur=ke_cur/(boltz*no_dof);
    
    TYPE0 ke_des=(boltz*3*atoms->tot_natms)*temperature;
    TYPE0 factor=sqrt(ke_des/ke_cur);
    TYPE0 facsq=ke_des/ke_cur;
    for(int i=0;i<3*natms;i++)
        x_d[i]*=factor;
    
    for(int i=0;i<6;i++)
        ke_curr[i]*=facsq;

    ke_cur*=facsq;
    t_cur*=facsq;
    delete [] temp;

}
/*--------------------------------------------
 create initial velocity
 --------------------------------------------*/
void MD_NH::init_vel(TYPE0 temperature)
{
    //TYPE0* x_d=(TYPE0*)atoms->vectors[atoms->find("x_d")].ret_vec();
    //int* type =(int*)atoms->vectors[atoms->find("type")].ret_vec();
    TYPE0* x_d;
    atoms->vectors[atoms->find("x_d")].ret(x_d);
    int* type;
    atoms->vectors[atoms->find("type")].ret(type);
    TYPE0* mass=atom_types->mass;
    
    int natms=atoms->natms;
    int icomp;
    TYPE0 temp[6];
    
    for(int i=0;i<6;i++)
        temp[i]=0.0;
    
    for(int i=0;i<natms;i++)
    {
        icomp=3*i;
        temp[0]+=mass[type[i]]*x_d[icomp]*x_d[icomp];
        temp[1]+=mass[type[i]]*x_d[icomp+1]*x_d[icomp+1];
        temp[2]+=mass[type[i]]*x_d[icomp+2]*x_d[icomp+2];
        temp[3]+=mass[type[i]]*x_d[icomp+1]*x_d[icomp+2];
        temp[4]+=mass[type[i]]*x_d[icomp]*x_d[icomp+2];
        temp[5]+=mass[type[i]]*x_d[icomp]*x_d[icomp+1];
    }
    
    
    for(int i=0;i<6;i++)
        ke_curr[i]=0.0;
    MPI_Allreduce(temp,ke_curr,6,MPI_TYPE0,MPI_SUM,world);
    
    ke_cur=(ke_curr[0]+ke_curr[1]+ke_curr[2]);
    t_cur=ke_cur/(boltz*3*atoms->tot_natms);
    
    if(ke_cur==0)
        error->abort("initial velocities are zero "
            "you should either use create_vel option "
            "or create the velocities in data file");
    
    TYPE0 ke_des=(boltz*no_dof)*temperature;
    TYPE0 factor=sqrt(ke_des/ke_cur);
    TYPE0 facsq=factor*factor;
    
    for(int i=0;i<3*natms;i++)
        x_d[i]*=factor;
    
    for(int i=0;i<6;i++)
        ke_curr[i]*=facsq;
    
    ke_cur*=facsq;
    t_cur*=facsq;
}


/*--------------------------------------------
 create initial velocity
 --------------------------------------------*/
void MD_NH::couple()
{
    TYPE0 tmp;
    if(chk_stress==XYZ)
    {
        tmp=v_per_atm[0]+v_per_atm[1]+v_per_atm[2];
        tmp=tmp/3.0;
        v_per_atm[0]=v_per_atm[1]=v_per_atm[2]=tmp;
        
        tmp=ke_curr[0]+ke_curr[1]+ke_curr[2];
        tmp=tmp/3.0;
        ke_curr[0]=ke_curr[1]=ke_curr[2]=tmp;
    }
    else if(chk_stress==YZ)
    {
        tmp=v_per_atm[1]+v_per_atm[2];
        tmp=tmp/2.0;
        v_per_atm[1]=v_per_atm[2]=tmp;

        tmp=ke_curr[1]+ke_curr[2];
        tmp=tmp/2.0;
        ke_curr[1]=ke_curr[2]=tmp;
    }
    else if(chk_stress==ZX)
    {
        tmp=v_per_atm[2]+v_per_atm[0];
        tmp=tmp/2.0;
        v_per_atm[2]=v_per_atm[0]=tmp;
        
        tmp=ke_curr[2]+ke_curr[0];
        tmp=tmp/2.0;
        ke_curr[2]=ke_curr[0]=tmp;

    }
    else if(chk_stress==XY)
    {
        tmp=v_per_atm[0]+v_per_atm[1];
        tmp=tmp/2.0;
        v_per_atm[0]=v_per_atm[1]=tmp;
        
        tmp=ke_curr[0]+ke_curr[1];
        tmp=tmp/2.0;
        ke_curr[0]=ke_curr[1]=tmp;

    }

}
