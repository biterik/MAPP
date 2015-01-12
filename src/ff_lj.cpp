/*--------------------------------------------
 Created by Sina on 07/15/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <cmath>
#include "ff_lj.h"
#include "memory.h"
#include "neighbor.h"
#include "atom_types.h"
using namespace MAPP_NS;
/*--------------------------------------------
 Fisher-Sinclair (FS) potential
 ref: 
 T. T. Lau, C. J. Forst, X. Lin, J. D. Gale,
 S. Yip, & K. J. Van Vliet
 Many-Body Potential for Point Defect Clusters
 in Fe-C Alloys
 Phys. Rev. Lett. Vol. 98, pp. 215501, 2007
 --------------------------------------------*/

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_lj::
ForceField_lj(MAPP* mapp) : ForceField(mapp)
{
    if(mapp->mode!=MD)
        error->abort("this forcefield works only with md mode");
    
    if(atoms->dimension!=3)
        error->abort("to use LJ potential, dimension of the box should be 3");
    arr_size=shift=0;
    
    int no_types=atom_types->no_types;
    
    int size=static_cast<int>((no_types+2)*(no_types+1)/2);
    GROW(cut_sq,arr_size,size);
    GROW(cut_sk_sq,arr_size,size);
    GROW(sigma,arr_size,size);
    GROW(epsilon,arr_size,size);
    GROW(offset,arr_size,size);
    GROW(chk_coef,arr_size,size);
    arr_size=size;
    
    for (int i=0;i<arr_size;i++)
        chk_coef[i]=0;
    
    CREATE1D(nrgy_strss,7);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_lj::~ForceField_lj()
{
    if(arr_size)
    {
        delete [] cut_sq;
        delete [] cut_sk_sq;
        delete [] sigma;
        delete [] epsilon;
        delete [] offset;
        delete [] chk_coef;
    }
    
    delete [] nrgy_strss;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_lj::coef(int narg,char** arg)
{
    //itype jtype epsilon sigma cuttoff
    
    if (narg!=6)
        error->abort("wrong coeff command for LJ FF");
    
    int ityp=atom_types->find_type(arg[1]);
    int jtyp=atom_types->find_type(arg[2]);
    int curs=COMP(ityp,jtyp);
    TYPE0 skin=atoms->skin;
    
    
    TYPE0 eps=epsilon[curs]=atof(arg[3]);
    TYPE0 sig=sigma[curs]=atof(arg[4]);
    TYPE0 cut=atof(arg[5]);

    if (eps<0.0)
    {
        error->abort("LJ Epsilon cannot be smaller than zero");
    }
    if (sig<0.0)
    {
        error->abort("LJ Sigma cannot be smaller than zero");
    }
    if (cut<0.0)
    {
        error->abort("LJ cutoff cannot be smaller than zero");
    }
    
    cut_sq[curs]=cut*cut;
    cut_sk_sq[curs]=(cut+skin)*(cut+skin);
    if(shift)
    {
        TYPE0 sig2=sig*sig/(cut*cut);
        TYPE0 sig6=sig2*sig2*sig2;
        TYPE0 sig12=sig6*sig6;
        
        offset[curs]=-4.0*eps*(sig12-sig6);
    }
    else
    {
        offset[curs]=0.0;
    }
    
    chk_coef[curs]=1;
    
    if (cut==0.0||sig==0.0||eps==0.0)
        chk_coef[curs]=0;
    else
        chk_coef[curs]=1;
}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_lj::init()
{
    TYPE0 skin=atoms->skin;
    TYPE0 ph_cut=0.0;
    
    for (int i=0;i<arr_size;i++)
        if (chk_coef[i])
        {
            cut_sk_sq[i]=cut_sq[i]+skin*skin
            +2*sqrt(cut_sq[i])*skin;
            ph_cut=MAX(ph_cut,sqrt(cut_sq[i]));
        }
    
    atoms->set_ph(ph_cut);
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    type_n=atoms->find("type");
    
    neighbor->pair_wise=1;

}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceField_lj::fin()
{
    
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceField_lj::
force_calc(int st_clc,TYPE0* en_st)
{
    /*
    TYPE0* x=(TYPE0*)atoms->vectors[x_n].ret_vec();
    TYPE0* f=(TYPE0*)atoms->vectors[f_n].ret_vec();
    int* type=(int*)atoms->vectors[type_n].ret_vec();
     */
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    int* type;
    atoms->vectors[type_n].ret(type);

    int natms=atoms->natms;
    int iatm,jatm;
    int itype,jtype,curs,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,sig,eps,csq;
    TYPE0 sig2,sig6,sig12,ft;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            curs=COMP(itype,jtype);
            if(chk_coef[curs])
            {
                icomp=3*iatm;
                jcomp=3*jatm;
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                
                rsq=dx0*dx0+dx1*dx1+dx2*dx2;
                
                csq=cut_sq[curs];
                if (rsq<csq)
                {
                    sig=sigma[curs];
                    eps=epsilon[curs];
                    sig2=sig*sig/rsq;
                    sig6=sig2*sig2*sig2;
                    sig12=sig6*sig6;
                    
                    ft=24.0*eps*(2.0*sig12-sig6)/rsq;
                    

                    f[icomp]+=ft*dx0;
                    f[icomp+1]+=ft*dx1;
                    f[icomp+2]+=ft*dx2;
                    if(jatm<natms)
                    {
                        nrgy_strss[0]+=4.0*eps*(sig12-sig6)
                        +offset[curs];
                        f[jcomp]-=ft*dx0;
                        f[jcomp+1]-=ft*dx1;
                        f[jcomp+2]-=ft*dx2;
                        if (st_clc)
                        {
                            nrgy_strss[1]-=ft*dx0*dx0;
                            nrgy_strss[2]-=ft*dx1*dx1;
                            nrgy_strss[3]-=ft*dx2*dx2;
                            nrgy_strss[4]-=ft*dx1*dx2;
                            nrgy_strss[5]-=ft*dx2*dx0;
                            nrgy_strss[6]-=ft*dx0*dx1;
                        }
                        
                    }
                    else
                    {
                        nrgy_strss[0]+=2.0*eps*(sig12-sig6)
                        +offset[curs]*0.5;
                        if (st_clc)
                        {
                            nrgy_strss[1]-=0.5*ft*dx0*dx0;
                            nrgy_strss[2]-=0.5*ft*dx1*dx1;
                            nrgy_strss[3]-=0.5*ft*dx2*dx2;
                            nrgy_strss[4]-=0.5*ft*dx1*dx2;
                            nrgy_strss[5]-=0.5*ft*dx2*dx0;
                            nrgy_strss[6]-=0.5*ft*dx0*dx1;
                        }
                    }                    
                }
            }
        }
    }
    
    if(st_clc)
    {
        for(int i=0;i<7;i++)
            en_st[i]=0.0;
        
        MPI_Allreduce(nrgy_strss,en_st,7,MPI_TYPE0,MPI_SUM,world);
    }
    else
    {
        en_st[0]=0.0;
        MPI_Allreduce(nrgy_strss,en_st,1,MPI_TYPE0,MPI_SUM,world);
    }
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
TYPE0 ForceField_lj::energy_calc()
{
    /*
    TYPE0* x=(TYPE0*)atoms->vectors[x_n].ret_vec();
    int* type=(int*)atoms->vectors[type_n].ret_vec();
     */
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    int* type;
    atoms->vectors[type_n].ret(type);
    
    int natms=atoms->natms;
    int iatm,jatm;
    
    int itype,jtype,curs,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,csq;
    TYPE0 eps,sig,sig2,sig6,sig12;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    TYPE0 en=0.0;
    TYPE0 en_tot=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            curs=COMP(itype,jtype);
            if(chk_coef[curs])
            {
                icomp=3*iatm;
                jcomp=3*jatm;
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                
                rsq=dx0*dx0+dx1*dx1+dx2*dx2;
                
                csq=cut_sq[curs];
                if (rsq<csq)
                {
                    sig=sigma[curs];
                    eps=epsilon[curs];
                    sig2=sig*sig/rsq;
                    sig6=sig2*sig2*sig2;
                    sig12=sig6*sig6;
                    
                    if(jatm<natms)
                        en+=4.0*eps*(sig12-sig6)
                        +offset[curs];
                    else
                        en+=2.0*eps*(sig12-sig6)
                        +offset[curs]*0.5;
                    
                }
            }
        }
    }
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}



