/*--------------------------------------------
 Created by Sina on 01/29/14.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <cmath>
#include "ff_fs.h"
#include "memory.h"
#include "neighbor.h"
#include "atoms.h"
#include "atom_types.h"
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 Finnis-Sinclair (FS) potential
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
ForceField_fs::
ForceField_fs(MAPP* mapp) : ForceField(mapp)
{
    max_pairs=0;
    
    if(atoms->dimension!=3)
        error->abort("to use FS potential, dimension of the box should be 3");
    
    int no_types=atom_types->no_types;
    CREATE2D(mat_t_1,no_types,no_types);
    CREATE2D(mat_t_2,no_types,no_types);
    CREATE1D(mat_A,no_types);
    for(int i=0;i<no_types;i++)
    {
        for(int j=0;j<no_types;j++)
            mat_t_1[i][j]=mat_t_2[i][j]=0.0;
        mat_A[i]=0.0;
    }
        
    
    int size=static_cast<int>((no_types+2)*(no_types+1)/2);
    
    arr_size=0;
    size=static_cast<int>((no_types+2)*(no_types+1)/2);
    GROW(cut_phi,arr_size,size);
    GROW(cut_rho,arr_size,size);
    
    GROW(mat_k_1,arr_size,size);
    GROW(mat_k_2,arr_size,size);
    GROW(mat_k_3,arr_size,size);
    
    GROW(cut_sq,arr_size,size);
    GROW(cut_sk_sq,arr_size,size);
    
    
    for(int i=0;i<size;i++)
        cut_phi[i]=cut_rho[i]
        =mat_k_1[i]=mat_k_2[i]=mat_k_3[i]
        =cut_sq[i]=cut_sk_sq[i]=0.0;
    
    GROW(chk_coef,arr_size,size);
    for(int i=0;i<size;i++)
        chk_coef[i]=0;
    
    
    arr_size=size;
    
    CREATE1D(nrgy_strss,7);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_fs::~ForceField_fs()
{
    if(atom_types->no_types)
    {
        for(int i=0;i<atom_types->no_types;i++)
        {
            delete [] mat_t_1[i];
            delete [] mat_t_2[i];
        }
        
        delete [] mat_t_1;
        delete [] mat_t_2;
        
        delete [] mat_A;
        
        delete [] cut_phi;
        delete [] cut_rho;
        
        delete [] mat_k_1;
        delete [] mat_k_2;
        delete [] mat_k_3;
        
        delete [] cut_sq;
        delete [] cut_sk_sq;
    }
    
    delete [] nrgy_strss;

}
/*--------------------------------------------
 read the force field file 
 --------------------------------------------*/
void ForceField_fs::coef(int narg,char** arg)
{
    if (narg!=2)
        error->abort("wrong coeff command "
        "for Finnis-Sinclair Force Field");
    
    /*A sample of the force field file:
     
     # The type of atoms in the force filed
     Fe C
     
     #section for A values it is equal to the number of types
     #the sequence must be kept like the above types
     #A_Fe
     1.8289905
     #A_C
     2.9588787
     
     #section for t_1 and t_2 values. no of entries=notyps^2
     #again sequence is important
     #Fe Fe
     1.0 0.504238
     #Fe C
     10.024001 1.638980
     #C Fe
     10.482408 3.782595
     #C C
     0.0 -7.329211
     
     #section for r_cut_rho r_cut_phi k_1 k_2 k_3. no of entries=notyps*(notyps+1)/2
     #Fe Fe
     3.569745 3.40 1.237115 -0.35921 -0.038560
     #Fe C
     2.545937 2.468801 8.972488 -4.086410 1.483233
     #C C
     2.892070 2.875598 22.061824 -17.468518 4.812639
     
     */
    
    if (atoms->my_p_no==0)
    {
        char* line;
        CREATE1D(line,MAXCHAR);
        int no_arg;
        int no_types=atom_types->no_types;
        char** argtrm = NULL;
        FILE* ff_file;
        ff_file=fopen(arg[1],"r");
        if(ff_file==NULL)
            error->abort("FS file %s not found",arg[1]);
        
        no_arg=0;
        while (no_arg==0)
        {
            fgets(line,MAXCHAR,ff_file);
            no_arg = mapp->parse_line(line,argtrm);
        }
        /*this should be the type of the atoms involved
        we should check wether we have such atomsw in
        our types
        */

        if(no_types!=no_arg)
            error->abort("number of types represented and forcefield do not match");
        
        
        int* types;
        CREATE1D(types,no_types);

                    /*
        for(int i=0;i<no_types;i++)
        {

            for(int j=0;j<atom_types->no_types;j++)
                if(!strcmp(argtrm[i],atom_types->atom_names[j]))
                    types[i]=j;
            if (types[i]==-1)
                error->abort("Atom type %s not found",argtrm[i]);

            
        }
       */
        for(int i=0;i<no_types;i++)
            types[i]=atom_types->find_type(argtrm[i]);
        
        for(int i=0;i<no_arg;i++)
            delete [] argtrm[i];
        if(no_arg)
            delete [] argtrm;
            
        
        for(int i=0;i<no_types;i++)
            for(int j=i+1;j<no_types;j++)
                if(types[i]==types[j])
                    error->abort("Atom type %s duplicated in the"
                    " force field file",argtrm[i]);
        
        
        // get the values for A matrix
        int no=0;
        while(no<no_types)
        {
            no_arg=0;
            while (no_arg==0)
            {
                fgets(line,MAXCHAR,ff_file);
                no_arg = mapp->parse_line(line,argtrm);
            }
            if(no_arg==1)
            {
                mat_A[types[no++]]=atof(argtrm[0]);
                
                for(int i=0;i<no_arg;i++)
                    delete [] argtrm[i];
                if(no_arg)
                    delete [] argtrm;
            }
            else
                error->abort("Wrong format in the force filed file");
        }
        
        
        //get the values for t_1 & t_2 matrix
        int no_0=0;
        while(no_0<no_types)
        {
            int no_1=0;
            while(no_1<no_types)
            {
                no_arg=0;
                while (no_arg==0)
                {
                    fgets(line,MAXCHAR,ff_file);
                    no_arg = mapp->parse_line(line,argtrm);
                }
                if(no_arg==2)
                {
                    mat_t_1[types[no_0]][types[no_1]]=atof(argtrm[0]);
                    mat_t_2[types[no_0]][types[no_1]]=atof(argtrm[1]);
                    no_1++;
                    
                    for(int i=0;i<no_arg;i++)
                        delete [] argtrm[i];
                    if(no_arg)
                        delete [] argtrm;
                }
                else
                    error->abort("Wrong format in the force filed file");
            }
            no_0++;
        }

        
        //get the values for cut_sq_rho, cut_sq_phi, k_1, k_2 & k_3 matrices
        int comp;
        no_0=0;
        while(no_0<no_types)
        {
            int no_1=no_0;
            while(no_1<no_types)
            {
                no_arg=0;
                while (no_arg==0)
                {
                    fgets(line,MAXCHAR,ff_file);
                    no_arg = mapp->parse_line(line,argtrm);
                }
                if(no_arg==5)
                {
                    comp=COMP(types[no_0],types[no_1]);
                    cut_rho[comp]=atof(argtrm[0]);
                    cut_phi[comp]=atof(argtrm[1]);
                    mat_k_1[comp]=atof(argtrm[2]);
                    mat_k_2[comp]=atof(argtrm[3]);
                    mat_k_3[comp]=atof(argtrm[4]);
                    no_1++;
                    
                    for(int i=0;i<no_arg;i++)
                        delete [] argtrm[i];
                    if(no_arg)
                        delete [] argtrm;
                }
                else
                    error->abort("Wrong format in the force filed file");
            }
            no_0++;
        }
        
        fclose(ff_file);
        
        
        
        for (int i=0;i<no_types;i++)
            for (int j=i;j<no_types;j++)
                chk_coef[COMP(types[i],types[j])]=1;
        
        delete [] line;
        delete [] types;
        
    }
    
    int no_types=atom_types->no_types;
    
    for(int i=0;i<no_types;i++)
    {
        MPI_Bcast(&mat_t_1[i][0],no_types,MPI_TYPE0,0,world);
        MPI_Bcast(&mat_t_2[i][0],no_types,MPI_TYPE0,0,world);
    }
    
    MPI_Bcast(&mat_A[0],no_types,MPI_TYPE0,0,world);

    int size=static_cast<int>((no_types+2)*(no_types+1)/2);
    
    MPI_Bcast(&cut_phi[0],size,MPI_TYPE0,0,world);
    MPI_Bcast(&cut_rho[0],size,MPI_TYPE0,0,world);
    MPI_Bcast(&mat_k_1[0],size,MPI_TYPE0,0,world);
    MPI_Bcast(&mat_k_2[0],size,MPI_TYPE0,0,world);
    MPI_Bcast(&mat_k_3[0],size,MPI_TYPE0,0,world);
    MPI_Bcast(&chk_coef[0],size,MPI_INT,0,world);
    


}
/*--------------------------------------------
 initiate before a run
 --------------------------------------------*/
void ForceField_fs::init()
{
    TYPE0 skin=atoms->skin;
    TYPE0 ph_cut=0.0;
    
    for (int i=0;i<arr_size;i++)
        if (chk_coef[i])
        {
            cut_sq[i]=MAX(cut_phi[i]*cut_phi[i],cut_rho[i]*cut_rho[i]);
            cut_sk_sq[i]=cut_sq[i]+(skin)*(skin)
            +2*sqrt(cut_sq[i])*(skin);
            ph_cut=MAX(ph_cut,sqrt(cut_sq[i]));
        }
    atoms->set_ph(ph_cut);
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    type_n=atoms->find("type");

    neighbor->pair_wise=1;
    
    rho_n=atoms->add<TYPE0>(1,1,"rho");
}
/*--------------------------------------------
 after a run
 --------------------------------------------*/
void ForceField_fs::fin()
{
    if(max_pairs)
    {
        delete [] drhoi_dr;
        delete [] drhoj_dr;
        max_pairs=0;
    }
    
    atoms->del(rho_n);
}
/*--------------------------------------------
 force and energy calculation
 --------------------------------------------*/
void ForceField_fs::
force_calc(int st_clc,TYPE0* en_st)
{
        if(max_pairs<neighbor->no_pairs)
    {
        if(max_pairs)
        {
            delete [] drhoi_dr;
            delete [] drhoj_dr;
        }
        
        max_pairs=neighbor->no_pairs;
        CREATE1D(drhoi_dr,max_pairs);
        CREATE1D(drhoj_dr,max_pairs);
    }

    
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    TYPE0* rho;
    atoms->vectors[rho_n].ret(rho);
    int* type;
    atoms->vectors[type_n].ret(type);
    
    int iatm,jatm;
    
    int itype,jtype,curs,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq;
    TYPE0 dr_rho,dr_phi,r,rho_coef,phi_coef,r_inv,rho_sqd;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    
    int istart;
    int natms=atoms->natms;
    
    for(iatm=0;iatm<natms;iatm++) rho[iatm]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            curs=COMP(itype,jtype);
            
            drhoi_dr[istart]=drhoj_dr[istart]=0.0;
            
            if(chk_coef[curs])
            {
                
                jcomp=3*jatm;
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                rsq=dx0*dx0+dx1*dx1+dx2*dx2;
                
                if(rsq<cut_sq[curs])
                {
                    r=sqrt(rsq);
                    r_inv=1.0/r;

                    if(r < cut_rho[curs])
                    {
                        dr_rho=r-cut_rho[curs];
                        rho[iatm]+=dr_rho*dr_rho*(mat_t_1[jtype][itype]
                                +mat_t_2[jtype][itype]*dr_rho);
                        
                        
                        drhoi_dr[istart]=-r_inv*dr_rho*(mat_t_1[jtype][itype]
                                    +1.5*mat_t_2[jtype][itype]*dr_rho);
                        
                        drhoj_dr[istart]=-r_inv*dr_rho*(mat_t_1[itype][jtype]
                                    +1.5*mat_t_2[itype][jtype]*dr_rho);
                        
                        if(jatm<natms)
                            rho[jatm]+=dr_rho*dr_rho*(mat_t_1[itype][jtype]
                                    +mat_t_2[itype][jtype]*dr_rho);
                        
                    }
                    
                    if(r < cut_phi[curs])
                    {
                        dr_phi=r-cut_phi[curs];
                        phi_coef=2.0*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq)
                        +dr_phi*dr_phi*(mat_k_2[curs]+2.0*mat_k_3[curs]*r);
                        phi_coef*=-r_inv;
                        
                        
                        f[icomp]+=dx0*phi_coef;
                        f[icomp+1]+=dx1*phi_coef;
                        f[icomp+2]+=dx2*phi_coef;
                        
                        
                        if(jatm<natms)
                        {
                            f[jcomp]-=dx0*phi_coef;
                            f[jcomp+1]-=dx1*phi_coef;
                            f[jcomp+2]-=dx2*phi_coef;
                            
                            nrgy_strss[0]+=dr_phi*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq);
                            if (st_clc)
                            {
                                nrgy_strss[1]+=phi_coef*dx0*dx0;
                                nrgy_strss[2]+=phi_coef*dx1*dx1;
                                nrgy_strss[3]+=phi_coef*dx2*dx2;
                                nrgy_strss[4]+=phi_coef*dx1*dx2;
                                nrgy_strss[5]+=phi_coef*dx2*dx0;
                                nrgy_strss[6]+=phi_coef*dx0*dx1;
                            }
                        }
                        else
                        {
                            nrgy_strss[0]+=0.5*dr_phi*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq);
                            if (st_clc)
                            {
                                nrgy_strss[1]+=0.5*phi_coef*dx0*dx0;
                                nrgy_strss[2]+=0.5*phi_coef*dx1*dx1;
                                nrgy_strss[3]+=0.5*phi_coef*dx2*dx2;
                                nrgy_strss[4]+=0.5*phi_coef*dx1*dx2;
                                nrgy_strss[5]+=0.5*phi_coef*dx2*dx0;
                                nrgy_strss[6]+=0.5*phi_coef*dx0*dx1;
                            }
                        }
                    }
                }
                
                
                
                
            }
            
            istart++;
        }
        
        rho_sqd=sqrt(rho[iatm]);
        nrgy_strss[0]-=mat_A[itype]*rho_sqd;
        rho[iatm]=-mat_A[itype]/rho_sqd;
    }
    
    atoms->update(rho_n);
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            if(drhoi_dr[istart]!=0.0 || drhoj_dr[istart]!=0.0)
            {
                jatm=neighbor_list[iatm][j];
                jcomp=3*jatm;
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                
                rho_coef=drhoi_dr[istart]*rho[iatm]
                +drhoj_dr[istart]*rho[jatm];
                f[icomp]+=dx0*rho_coef;
                f[icomp+1]+=dx1*rho_coef;
                f[icomp+2]+=dx2*rho_coef;
                
                if(jatm<natms)
                {
                    f[jcomp]-=dx0*rho_coef;
                    f[jcomp+1]-=dx1*rho_coef;
                    f[jcomp+2]-=dx2*rho_coef;
                    
                    if (st_clc)
                    {
                        nrgy_strss[1]+=rho_coef*dx0*dx0;
                        nrgy_strss[2]+=rho_coef*dx1*dx1;
                        nrgy_strss[3]+=rho_coef*dx2*dx2;
                        nrgy_strss[4]+=rho_coef*dx1*dx2;
                        nrgy_strss[5]+=rho_coef*dx2*dx0;
                        nrgy_strss[6]+=rho_coef*dx0*dx1;
                    }
                }
                else
                {
                    if (st_clc)
                    {
                        nrgy_strss[1]+=0.5*rho_coef*dx0*dx0;
                        nrgy_strss[2]+=0.5*rho_coef*dx1*dx1;
                        nrgy_strss[3]+=0.5*rho_coef*dx2*dx2;
                        nrgy_strss[4]+=0.5*rho_coef*dx1*dx2;
                        nrgy_strss[5]+=0.5*rho_coef*dx2*dx0;
                        nrgy_strss[6]+=0.5*rho_coef*dx0*dx1;
                    }
                }
                
            }
            
            istart++;

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
        MPI_Allreduce(nrgy_strss,en_st,1,MPI_TYPE0,MPI_SUM,world);
    }
    
    
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
TYPE0 ForceField_fs::energy_calc()
{
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* rho;
    atoms->vectors[rho_n].ret(rho);
    int* type;
    atoms->vectors[type_n].ret(type);
    
    int natms=atoms->natms;
    int iatm,jatm;
    
    int itype,jtype,curs,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,csq;
    TYPE0 dr_rho,dr_phi,r;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    for(int i=0;i<natms;i++)
        rho[i]=0.0;
    
    TYPE0 en=0.0;
    TYPE0 en_tot=0.0;
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            curs=COMP(itype,jtype);
            if(chk_coef[curs])
            {
                jcomp=3*jatm;
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                rsq=dx0*dx0+dx1*dx1+dx2*dx2;
                csq=cut_sq[curs];
                
                if(rsq < csq)
                {
                    r=sqrt(rsq);
                    if(r < cut_rho[curs])
                    {
                        dr_rho=r-cut_rho[curs];
                        rho[iatm]+=dr_rho*dr_rho*(mat_t_1[jtype][itype]
                        +mat_t_2[jtype][itype]*dr_rho);
                        if(jatm<natms)
                            rho[jatm]+=dr_rho*dr_rho*(mat_t_1[itype][jtype]
                            +mat_t_2[itype][jtype]*dr_rho);
                    }
                    
                    if(r < cut_phi[curs])
                    {
                        dr_phi=r-cut_phi[curs];
                        if(jatm<natms)
                        {
                            en+=dr_phi*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq);
                        }
                        else
                        {
                            en+=0.5*dr_phi*dr_phi*(mat_k_1[curs]+mat_k_2[curs]*r+mat_k_3[curs]*rsq);
                        }
                    }
                }
            }
        }
        en-=mat_A[itype]*sqrt(rho[iatm]);
    }

    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    return en_tot;
}


