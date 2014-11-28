#include "min_cg.h"
#include "neighbor.h"
#include "ff.h"
#include "thermo_dynamics.h"
#include "write.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min_CG::Min_CG(MAPP* mapp,int narg,char** arg):Min(mapp)
{
        
    CREATE2D(H_dof,dim,dim);
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            H_dof[i][j]=0;
    
    int icmp;
    int jcmp;
    int iarg=2;
    while(iarg<narg)
    {
        if(!strcmp(arg[iarg],"max_iteration"))
        {
            iarg++;
            if(iarg==narg)
                error->abort("maximum iteration not defined");
            max_iter=atoi(arg[iarg]);
            iarg++;
        }
        else if(!strcmp(arg[iarg],"energy_tol"))
        {
            iarg++;
            if(iarg==narg)
                error->abort("energy tolerance not defined");
            energy_tolerance=atof(arg[iarg]);
            iarg++;
        }
        else if(sscanf(arg[iarg],"H[%d][%d]",&icmp,&jcmp)==2)
        {
            if(icmp<0 || icmp>=atoms->dimension)
                error->abort("wrong command %s",arg[iarg]);
            if(jcmp<0 || jcmp>=atoms->dimension)
                error->abort("wrong command %s",arg[iarg]);
            
            if(icmp<=jcmp)
                H_dof[jcmp][icmp]=1;
            else
                H_dof[icmp][jcmp]=1;
            iarg++;
        }
        else
            error->abort("wrong command %s",arg[iarg]);
    }
    
    if(max_iter<=0)
        error->abort("maximum iteration cannot be equal or less than zero");
    if(energy_tolerance<=0.0)
        error->abort("energy tolerance cannot be equal or less than zero");
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Min_CG::~Min_CG()
{
    delete thermo;
    for(int i=0;i<dim;i++)
        delete [] H_dof[i];
    delete [] H_dof;
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void Min_CG::init()
{

    x_dim=atoms->vectors[0].dim;
    f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<TYPE0>(0,x_dim,"f");
    
    if(mapp->mode==DMD)
        type_n=atoms->find("c");
    else
        type_n=atoms->find("type");
    
    x_prev_n=atoms->add<TYPE0>(0,x_dim,"x_prev");
    f_prev_n=atoms->add<TYPE0>(0,x_dim,"f_prev");
    h_n=atoms->add<TYPE0>(0,x_dim,"h");

    id_n=atoms->find("id");
    
    dof_n=atoms->find_exist("dof");
    if(dof_n<0)
        vecs_comm=new VecLst(mapp,7,0,type_n,f_n,x_prev_n,f_prev_n,h_n,id_n);
    else
        vecs_comm=new VecLst(mapp,8,0,type_n,f_n,x_prev_n,f_prev_n,h_n,dof_n,id_n);

    vecs_comm->add_update(0);
    atoms->reset_comm(vecs_comm);

    
    forcefield->init();
    atoms->ph_setup(1,vecs_comm);
    neighbor->init();
    neighbor->create_list(0,1);
    atoms->store_0();

    line_search=new LineSearch_BackTrack(mapp,vecs_comm);
    line_search->h_n=h_n;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            if(H_dof[i][j])
                chng_box=1;
    
    
    if(chng_box)
    {
        CREATE2D(h_H,dim,dim);
        CREATE2D(f_H,dim,dim);
        CREATE2D(f_H_prev,dim,dim);
        CREATE2D(H_prev,dim,dim);
        CREATE2D(B_prev,dim,dim);
        
        line_search->h_H=h_H;
        line_search->f_H=f_H;
        line_search->H_prev=H_prev;
        line_search->B_prev=B_prev;
    }
    line_search->chng_box=chng_box;
    TYPE0* f;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    if(chng_box)
    {
        TYPE0** stress;
        CREATE2D(stress,dim,dim);
        TYPE0** H=atoms->H;
        TYPE0** B=atoms->B;
        
        atoms->vectors[f_n].ret(f);
        for(int i=0;i<x_dim*atoms->natms;i++)
            f[i]=0.0;
        
        forcefield->force_calc(1,energy_stress);
        rectify_f(f);
        curr_energy=energy_stress[0];
        thermo->update(pe_idx,energy_stress[0]);
        thermo->update(stress_idx,6,&energy_stress[1]);
        
        stress[0][0]=energy_stress[1];
        stress[1][1]=energy_stress[2];
        stress[2][2]=energy_stress[3];
        stress[1][2]=stress[2][1]=energy_stress[4];
        stress[0][2]=stress[2][0]=energy_stress[5];
        stress[0][1]=stress[1][0]=energy_stress[6];
        
        for(int i=0;i<dim;i++)
        {
            for(int j=0;j<dim;j++)
            {
                f_H[j][i]=0.0;
                for(int k=0;k<dim;k++)
                    f_H[j][i]+=stress[i][k]*B[k][j];
            }
        }
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                if(H_dof[i][j]==0)
                    f_H[i][j]=0.0;
                
        int icomp=0;
        
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<dim;j++)
            {
                f[icomp+j]=f[icomp+j]*H[j][j];
                for(int k=j+1;k<dim;k++)
                    f[icomp+j]+=H[k][j]*f[icomp+k];
            }
            icomp+=x_dim;
        }
        
        for(int i=0;i<dim;i++)
            delete [] stress[i];
        delete [] stress;
    }
    else
    {
        atoms->vectors[f_n].ret(f);
        for(int i=0;i<x_dim*atoms->natms;i++)
            f[i]=0.0;
        
        forcefield->force_calc(1,energy_stress);
        rectify_f(f);
        curr_energy=energy_stress[0];
        thermo->update(pe_idx,energy_stress[0]);
        thermo->update(stress_idx,6,&energy_stress[1]);
        
    }
    if(write!=NULL)
        write->init();
    thermo->init();
    delete [] energy_stress;


}
/*--------------------------------------------
 min
 --------------------------------------------*/
void Min_CG::run()
{
    
    TYPE0* x;
    TYPE0* f;
    TYPE0* x_0;
    TYPE0* f_0;
    TYPE0* h;
    TYPE0 prev_energy;
    TYPE0 alpha;
    int size;
    int istp=0;
    TYPE0 f0_f0;
    TYPE0 f_f;
    TYPE0 f_f0;
    TYPE0 ratio;
    TYPE0 inner;
    TYPE0** H;
    TYPE0** B;
    err=LS_S;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    
    if(chng_box)
    {
        H=atoms->H;
        B=atoms->B;
        TYPE0** stress;
        CREATE2D(stress,dim,dim);
        
        atoms->vectors[f_n].ret(f);
        atoms->vectors[h_n].ret(h);
        memcpy(h,f,x_dim*atoms->natms*sizeof(TYPE0));
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                h_H[i][j]=f_H[i][j];
        
        inner=0.0;
        for(int i=0;i<atoms->natms*x_dim;i++)
            inner+=f[i]*f[i];
        f0_f0=0.0;
        MPI_Allreduce(&inner,&f0_f0,1,MPI_TYPE0,MPI_SUM,world);
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                if(H_dof[i][j])
                    f0_f0+=f_H[i][j]*f_H[i][j];
        
        while(err==LS_S)
        {
            
            if(f0_f0==0.0)
            {
                err=LS_F_GRAD0;
                continue;
            }
            
            atoms->vectors[0].ret(x);
            atoms->vectors[f_n].ret(f);
            atoms->vectors[x_prev_n].ret(x_0);
            atoms->vectors[f_prev_n].ret(f_0);
            size=atoms->natms*x_dim*sizeof(TYPE0);
            
            memcpy(x_0,x,size);
            memcpy(f_0,f,size);
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<dim;j++)
                {
                    H_prev[i][j]=H[i][j];
                    B_prev[i][j]=B[i][j];
                    f_H_prev[i][j]=f_H[i][j];
                }
            }

            prev_energy=curr_energy;
            
            if(write!=NULL)
                write->write();
            thermo->thermo_print();
            err=line_search->line_min(curr_energy,alpha);
            if(err==LS_S)
                if(prev_energy-curr_energy<energy_tolerance)
                    err=MIN_F_TOLERANCE;
            if(err==LS_S)
                if(istp==max_iter)
                    err=MIN_F_MAX_ITER;
            
            atoms->vectors[f_n].ret(f);
            
            for(int i=0;i<x_dim*atoms->natms;i++)
                f[i]=0.0;
            
            forcefield->force_calc(1,energy_stress);
            rectify_f(f);
            if(thermo->test_prev_step() || err)
            {
                thermo->update(pe_idx,energy_stress[0]);
                thermo->update(stress_idx,6,&energy_stress[1]);
            }
            
            if(err)
                continue;
            
            stress[0][0]=energy_stress[1];
            stress[1][1]=energy_stress[2];
            stress[2][2]=energy_stress[3];
            stress[1][2]=stress[2][1]=energy_stress[4];
            stress[0][2]=stress[2][0]=energy_stress[5];
            stress[0][1]=stress[1][0]=energy_stress[6];
            
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<dim;j++)
                {
                    f_H[j][i]=0.0;
                    for(int k=0;k<dim;k++)
                        f_H[j][i]+=stress[i][k]*B[k][j];
                }
            }
            
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                    if(H_dof[i][j]==0)
                        f_H[i][j]=0.0;
            
            atoms->vectors[f_n].ret(f);

            int icomp=0;
            
            //change f to fractional coordinates
            for(int i=0;i<atoms->natms;i++)
            {
                for(int j=0;j<dim;j++)
                {
                    f[icomp+j]=f[icomp+j]*H[j][j];
                    for(int k=j+1;k<dim;k++)
                        f[icomp+j]+=H[k][j]*f[icomp+k];
                }
                icomp+=x_dim;
            }
            
            
            atoms->vectors[h_n].ret(h);
            atoms->vectors[f_prev_n].ret(f_0);
            inner=0.0;
            for(int i=0;i<x_dim*atoms->natms;i++)
                inner+=f[i]*f[i];
            
            f_f=0.0;
            MPI_Allreduce(&inner,&f_f,1,MPI_TYPE0,MPI_SUM,world);
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                    if(H_dof[i][j])
                        f_f+=f_H[i][j]*f_H[i][j];
            
            inner=0.0;
            for(int i=0;i<x_dim*atoms->natms;i++)
                inner+=f[i]*f_0[i];
            f_f0=0.0;
            MPI_Allreduce(&inner,&f_f0,1,MPI_TYPE0,MPI_SUM,world);
            for(int i=0;i<dim;i++)
                for(int j=0;j<dim;j++)
                    if(H_dof[i][j])
                        f_f0+=f_H_prev[i][j]*f_H[i][j];
            
            ratio=(f_f-f_f0)/(f0_f0);
            
            for(int i=0;i<x_dim*atoms->natms;i++)
            {
                h[i]*=ratio;
                h[i]+=f[i];
            }
            
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<dim;j++)
                {
                    if(H_dof[i][j])
                    {
                        h_H[i][j]*=ratio;
                        h_H[i][j]+=f_H[i][j];
                    }
                }
            }
            
            f0_f0=f_f;
            istp++;
            step_no++;
        }
        
        for(int i=0;i<dim;i++)
            delete [] stress[i];
        delete [] stress;
        
        
    }
    else
    {
        
        atoms->vectors[h_n].ret(h);
        atoms->vectors[f_n].ret(f);
        memcpy(h,f,x_dim*atoms->natms*sizeof(TYPE0));
        
        atoms->vectors[f_n].ret(f);
        inner=0.0;
        for(int i=0;i<x_dim*atoms->natms;i++)
            inner+=f[i]*f[i];
        f0_f0=0.0;
        MPI_Allreduce(&inner,&f0_f0,1,MPI_TYPE0,MPI_SUM,world);
        
        while(err==LS_S)
        {
            
            if(f0_f0==0.0)
            {
                err=LS_F_GRAD0;
                continue;
            }
            
            atoms->vectors[0].ret(x);
            atoms->vectors[f_n].ret(f);
            atoms->vectors[x_prev_n].ret(x_0);
            atoms->vectors[f_prev_n].ret(f_0);
            size=atoms->natms*x_dim*sizeof(TYPE0);
            
            memcpy(x_0,x,size);
            memcpy(f_0,f,size);
            
            prev_energy=curr_energy;
            
            if(write!=NULL)
                write->write();
            
            thermo->thermo_print();
            err=line_search->line_min(curr_energy,alpha);
            
            if(err==LS_S)
                if(prev_energy-curr_energy<energy_tolerance)
                    err=MIN_F_TOLERANCE;
            if(err==LS_S)
                if(istp==max_iter)
                    err=MIN_F_MAX_ITER;
            
            
            atoms->vectors[f_n].ret(f);
            atoms->vectors[f_prev_n].ret(f_0);
            atoms->vectors[h_n].ret(h);
            
            for(int i=0;i<x_dim*atoms->natms;i++)
                f[i]=0.0;
            
            if(thermo->test_prev_step() || err)
            {
                forcefield->force_calc(1,energy_stress);
                rectify_f(f);
                thermo->update(pe_idx,energy_stress[0]);
                thermo->update(stress_idx,6,&energy_stress[1]);
                curr_energy=energy_stress[0];
            }
            else
            {
                forcefield->force_calc(0,&curr_energy);
                rectify_f(f);
            }

            if(err)
                continue;
            
            inner=0.0;
            for(int i=0;i<x_dim*atoms->natms;i++)
                inner+=f[i]*f[i];
            f_f=0.0;
            MPI_Allreduce(&inner,&f_f,1,MPI_TYPE0,MPI_SUM,world);
            
            inner=0.0;
            for(int i=0;i<x_dim*atoms->natms;i++)
                inner+=f[i]*f_0[i];
            f_f0=0.0;
            MPI_Allreduce(&inner,&f_f0,1,MPI_TYPE0,MPI_SUM,world);
            
            ratio=(f_f-f_f0)/(f0_f0);
            
            for(int i=0;i<x_dim*atoms->natms;i++)
            {
                h[i]*=ratio;
                h[i]+=f[i];
            }
            
            f0_f0=f_f;
            istp++;
            step_no++;
            
        }
    }
    
    delete [] energy_stress;
   
}
/*--------------------------------------------
 finishing minimization
 --------------------------------------------*/
void Min_CG::fin()
{
    forcefield->fin();
    neighbor->fin();
    
    atoms->del(h_n);
    atoms->del(f_prev_n);
    atoms->del(x_prev_n);

    delete vecs_comm;
    delete line_search;
    
    if(chng_box)
    {
        for(int i=0;i<dim;i++)
            delete [] f_H_prev[i];
        delete [] f_H_prev;
        
        for(int i=0;i<dim;i++)
            delete [] f_H[i];
        delete [] f_H;
        
        for(int i=0;i<dim;i++)
            delete [] h_H[i];
        delete [] h_H;
        
        for(int i=0;i<dim;i++)
            delete [] H_prev[i];
        delete [] H_prev;
        
        for(int i=0;i<dim;i++)
            delete [] B_prev[i];
        delete [] B_prev;

    }
    
    if(write!=NULL)
        write->fin();
    
    thermo->fin();
    errors();
    
    atoms->x2s(atoms->natms);
}
