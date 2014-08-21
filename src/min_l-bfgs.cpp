/*--------------------------------------------
 Created by Sina on 07/02/14.
 Copyright (c) 2013 MIT. All rights reserved.
 
 L-BFGS minimization is written based on
 Numerical Optimization written by Nocedal & 
 Wright, second edition, pages 177-179, 
 Algorithm 7.4 & 7.5 Equation (7.20)
 
 with respect to notations in Nocedal:
 new_y_i=-y_i
 new_rho_i=-rho_i
 new_alpha_i=-alpha_i
 new_beta=-beta
 --------------------------------------------*/
#include "min_l-bfgs.h"
#include "neighbor.h"
#include "ff.h"
#include "thermo_dynamics.h"
#include "write.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min_LBFGS::Min_LBFGS(MAPP* mapp,int narg,char** arg):Min(mapp)
{

    m_it=2;
    
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
        else if(!strcmp(arg[iarg],"m"))
        {
            iarg++;
            if(iarg==narg)
                error->abort("wrong command");
            m_it=atoi(arg[iarg]);
            iarg++;
        }
        else if(sscanf(arg[iarg],"H[%d][%d]",&icmp,&jcmp)==2)
        {
            if(icmp<0 || icmp>=dim)
                error->abort("wrong command");
            if(jcmp<0 || jcmp>=dim)
                error->abort("wrong command");
            
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
Min_LBFGS::~Min_LBFGS()
{
    delete thermo;
    for(int i=0;i<dim;i++)
        delete [] H_dof[i];
    delete [] H_dof;
}
/*--------------------------------------------
 init before a run
 --------------------------------------------*/
void Min_LBFGS::init()
{
    char* vec_name;
    
    CREATE1D(s_list,m_it);
    CREATE1D(y_list,m_it);
    x_dim=atoms->vectors[0].dim;
    type_n=atoms->find("type");
    f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<TYPE0>(0, 3,"f");
    
    for (int i=0;i<m_it;i++)
    {
        CREATE1D(vec_name,20);
        sprintf(vec_name,"s_%d",i);
        s_list[i]=atoms->add<TYPE0>(0,x_dim,vec_name);
        delete [] vec_name;
        
        CREATE1D(vec_name,20);
        sprintf(vec_name,"y_%d",i);
        y_list[i]=atoms->add<TYPE0>(0,x_dim,vec_name);
        delete [] vec_name;
    }
    
    x_prev_n=atoms->add<TYPE0>(0,x_dim,"x_prev");
    f_prev_n=atoms->add<TYPE0>(0,x_dim,"f_prev");
    h_n=atoms->add<TYPE0>(0,x_dim,"h");
    
    
    int* tmp_lst;
    CREATE1D(tmp_lst,2*m_it+6);
    int icurs=0;
    tmp_lst[icurs++]=0;
    tmp_lst[icurs++]=type_n;
    tmp_lst[icurs++]=f_n;
    for (int i=0;i<m_it;i++)
    {
        tmp_lst[icurs++]=s_list[i];
        tmp_lst[icurs++]=y_list[i];
    }
    tmp_lst[icurs++]=x_prev_n;
    tmp_lst[icurs++]=f_prev_n;
    tmp_lst[icurs++]=h_n;
    
    vecs_comm=new VecLst(mapp,tmp_lst,icurs);
    delete [] tmp_lst;
    vecs_comm->add_update(0);
    
    atoms->reset_comm(vecs_comm);
    
    forcefield->init();
    atoms->ph_setup(1,vecs_comm);
    neighbor->init();
    neighbor->create_list(0,1);
    atoms->store_0();
    CREATE1D(rho,m_it);
    CREATE1D(alpha,m_it);
    line_search=new LineSearch_BackTrack(mapp,vecs_comm);
    line_search->h_n=h_n;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            if(H_dof[i][j])
                chng_box=1;
    
    
    if(chng_box)
    {
        CREATE1D(H_y,m_it);
        CREATE1D(H_s,m_it);
        for(int i=0;i<m_it;i++)
        {
            CREATE2D(H_y[i],dim,dim);
            CREATE2D(H_s[i],dim,dim);
        }
        
        CREATE1D(H_s_y_list,m_it);
        for(int i=0;i<m_it;i++)
            H_s_y_list[i]=i;
        
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
 run
 --------------------------------------------*/
void Min_LBFGS::run()
{
    TYPE0* f;
    TYPE0* x;
    TYPE0* x_0;
    TYPE0* f_0;
    TYPE0* h;
    TYPE0* tmp_vec0;
    TYPE0* tmp_vec1;
    TYPE0** H;
    TYPE0** B;
    TYPE0 inner;
    TYPE0 tot_inner;
    TYPE0 prev_energy;
    TYPE0 inner_tmp;
    TYPE0 tot_inner_tmp;
    TYPE0 ratio=1.0;
    TYPE0 alpha_m;
    int s_list_tmp,y_list_tmp,s_y_list_tmp;
    int k_it=0;
    int size;
    int istp=0;
    err=LS_S;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    
    if(chng_box)
    {

        H=atoms->H;
        B=atoms->B;
        TYPE0** stress;
        CREATE2D(stress,dim,dim);
        while(err==LS_S)
        {
            atoms->vectors[h_n].ret(h);
            atoms->vectors[0].ret(x);
            atoms->vectors[f_n].ret(f);
            atoms->vectors[x_prev_n].ret(x_0);
            atoms->vectors[f_prev_n].ret(f_0);
            size=atoms->natms*x_dim*sizeof(TYPE0);
            
            memcpy(h,f,size);
            memcpy(x_0,x,size);
            memcpy(f_0,f,size);
            
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<dim;j++)
                {
                    h_H[i][j]=f_H[i][j];
                    f_H_prev[i][j]=f_H[i][j];
                    H_prev[i][j]=H[i][j];
                    B_prev[i][j]=B[i][j];
                }
            }
            
            inner=0.0;
            for(int i=k_it-1;i>-1;i--)
            {
                
                atoms->vectors[s_list[i]].ret(tmp_vec0);
                inner=0.0;
                for(int j=0;j<atoms->natms*x_dim;j++)
                    inner+=h[j]*tmp_vec0[j];
                
                alpha[i]=0.0;
                MPI_Allreduce(&inner,&(alpha[i]),1,MPI_TYPE0,MPI_SUM,world);
                
                for(int j=0;j<dim;j++)
                {
                    for(int k=0;k<dim;k++)
                    {
                        if(H_dof[j][k])
                            alpha[i]+=h_H[j][k]*H_s[H_s_y_list[i]][j][k];
                    }
                }
                
                alpha[i]*=rho[i];
                
                atoms->vectors[y_list[i]].ret(tmp_vec0);
                for(int j=0;j<atoms->natms*x_dim;j++)
                    h[j]-=alpha[i]*tmp_vec0[j];
                
                for(int j=0;j<dim;j++)
                {
                    for(int k=0;k<dim;k++)
                    {
                        if(H_dof[j][k])
                            h_H[j][k]-=alpha[i]*H_y[H_s_y_list[i]][j][k];
                    }
                }
            }
            
            
            for(int j=0;j<atoms->natms*x_dim;j++)
                h[j]*=ratio;
            for(int j=0;j<dim;j++)
            {
                for(int k=0;k<dim;k++)
                {
                    if(H_dof[j][k])
                        h_H[j][k]*=ratio;
                }
            }
            
            for(int i=0;i<k_it;i++)
            {
                
                atoms->vectors[y_list[i]].ret(tmp_vec0);
                inner=0.0;
                for(int j=0;j<atoms->natms*x_dim;j++)
                    inner+=h[j]*tmp_vec0[j];
                
                tot_inner=0.0;
                MPI_Allreduce(&inner,&tot_inner,1,MPI_TYPE0,MPI_SUM,world);
                
                for(int j=0;j<dim;j++)
                {
                    for(int k=0;k<dim;k++)
                    {
                        if(H_dof[j][k])
                            tot_inner+=h_H[j][k]*H_y[H_s_y_list[i]][j][k];
                    }
                }
                
                tot_inner*=rho[i];
                
                atoms->vectors[s_list[i]].ret(tmp_vec0);
                for(int j=0;j<atoms->natms*x_dim;j++)
                    h[j]-=(tot_inner+alpha[i])*tmp_vec0[j];
                
                for(int j=0;j<dim;j++)
                    for(int k=0;k<dim;k++)
                        if(H_dof[j][k])
                            h_H[j][k]-=(tot_inner+alpha[i])*H_s[H_s_y_list[i]][j][k];
                
            }
            
            

            prev_energy=curr_energy;
            
            if(write!=NULL)
                write->write();
            thermo->thermo_print();
            
            err=line_search->line_min(curr_energy,alpha_m);
            
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
            
            
            
            if(m_it)
            {
                if(k_it==m_it)
                {
                    s_list_tmp=s_list[0];
                    y_list_tmp=y_list[0];
                    s_y_list_tmp=H_s_y_list[0];
                    for(int i=0;i<m_it-1;i++)
                    {
                        s_list[i]=s_list[i+1];
                        y_list[i]=y_list[i+1];
                        H_s_y_list[i]=H_s_y_list[i+1];
                        rho[i]=rho[i+1];
                    }
                    s_list[m_it-1]=s_list_tmp;
                    y_list[m_it-1]=y_list_tmp;
                    H_s_y_list[m_it-1]=s_y_list_tmp;
                    
                }
                else
                    k_it++;
                
                
                atoms->vectors[0].ret(x);
                atoms->vectors[f_n].ret(f);
                atoms->vectors[h_n].ret(h);
                atoms->vectors[x_prev_n].ret(x_0);
                atoms->vectors[f_prev_n].ret(f_0);
                atoms->vectors[s_list[k_it-1]].ret(tmp_vec0);
                atoms->vectors[y_list[k_it-1]].ret(tmp_vec1);
                
                inner=0.0;
                inner_tmp=0.0;
                for(int i=0;i<x_dim*atoms->natms;i++)
                {
                    tmp_vec0[i]=alpha_m*h[i];
                    tmp_vec1[i]=f[i]-f_0[i];
                    inner+=tmp_vec0[i]*tmp_vec1[i];
                    inner_tmp+=tmp_vec1[i]*tmp_vec1[i];
                }
                tot_inner=0.0;
                MPI_Allreduce(&inner,&tot_inner,1,MPI_TYPE0,MPI_SUM,world);
                tot_inner_tmp=0.0;
                MPI_Allreduce(&inner_tmp,&tot_inner_tmp,1,MPI_TYPE0,MPI_SUM,world);
                
                for(int i=0;i<dim;i++)
                {
                    for(int j=0;j<dim;j++)
                    {
                        H_s[H_s_y_list[k_it-1]][i][j]=alpha_m*h_H[i][j];
                        H_y[H_s_y_list[k_it-1]][i][j]=f_H[i][j]-f_H_prev[i][j];
                        tot_inner+=H_s[H_s_y_list[k_it-1]][i][j]*H_y[H_s_y_list[k_it-1]][i][j];
                        tot_inner_tmp+=H_s[H_s_y_list[k_it-1]][i][j]*H_s[H_s_y_list[k_it-1]][i][j];
                    }
                }
                
                ratio=-tot_inner/tot_inner_tmp;
                rho[k_it-1]=1.0/tot_inner;
            }
            else
            {
                atoms->vectors[f_n].ret(f);
                atoms->vectors[h_n].ret(h);
                atoms->vectors[f_prev_n].ret(f_0);
                inner=0.0;
                inner_tmp=0.0;
                for(int i=0;i<x_dim*atoms->natms;i++)
                {
                    inner+=(alpha_m*h[i])*(f[i]-f_0[i]);
                    inner_tmp+=(f[i]-f_0[i])*(f[i]-f_0[i]);
                }
                tot_inner=0.0;
                MPI_Allreduce(&inner,&tot_inner,1,MPI_TYPE0,MPI_SUM,world);
                tot_inner_tmp=0.0;
                MPI_Allreduce(&inner_tmp,&tot_inner_tmp,1,MPI_TYPE0,MPI_SUM,world);
                
                for(int i=0;i<dim;i++)
                {
                    for(int j=0;j<dim;j++)
                    {
                        tot_inner+=alpha_m*h_H[i][j]*(f_H[i][j]-f_H_prev[i][j]);
                        tot_inner_tmp+=(f_H[i][j]-f_H_prev[i][j])*(f_H[i][j]-f_H_prev[i][j]);
                    }
                }
                
                ratio=-tot_inner/tot_inner_tmp;
            }
            
            istp++;
            step_no++;
        }
        
        
        for(int i=0;i<dim;i++)
            delete [] stress[i];
        delete [] stress;
        
        
    }
    else
    {
        while(err==LS_S)
        {
            atoms->vectors[h_n].ret(h);
            atoms->vectors[0].ret(x);
            atoms->vectors[f_n].ret(f);
            atoms->vectors[x_prev_n].ret(x_0);
            atoms->vectors[f_prev_n].ret(f_0);
            size=atoms->natms*x_dim*sizeof(TYPE0);
            
            memcpy(h,f,size);
            memcpy(x_0,x,size);
            memcpy(f_0,f,size);

            
            inner=0.0;
            for(int i=k_it-1;i>-1;i--)
            {
                
                atoms->vectors[s_list[i]].ret(tmp_vec0);
                inner=0.0;
                for(int j=0;j<atoms->natms*x_dim;j++)
                    inner+=h[j]*tmp_vec0[j];
                
                alpha[i]=0.0;
                MPI_Allreduce(&inner,&(alpha[i]),1,MPI_TYPE0,MPI_SUM,world);
                alpha[i]*=rho[i];
                
                
                atoms->vectors[y_list[i]].ret(tmp_vec0);
                for(int j=0;j<atoms->natms*x_dim;j++)
                    h[j]-=alpha[i]*tmp_vec0[j];
            }
            
            
            for(int j=0;j<atoms->natms*x_dim;j++)
                h[j]*=ratio;
            
            
            for(int i=0;i<k_it;i++)
            {
                
                atoms->vectors[y_list[i]].ret(tmp_vec0);
                inner=0.0;
                for(int j=0;j<atoms->natms*x_dim;j++)
                    inner+=h[j]*tmp_vec0[j];
                
                tot_inner=0.0;
                MPI_Allreduce(&inner,&tot_inner,1,MPI_TYPE0,MPI_SUM,world);
                tot_inner*=rho[i];
                
                atoms->vectors[s_list[i]].ret(tmp_vec0);
                for(int j=0;j<atoms->natms*x_dim;j++)
                    h[j]-=(tot_inner+alpha[i])*tmp_vec0[j];
            }
            
            prev_energy=curr_energy;
            
            if(write!=NULL)
                write->write();
            thermo->thermo_print();
            
            err=line_search->line_min(curr_energy,alpha_m);
            
            if(err==LS_S)
                if(prev_energy-curr_energy<energy_tolerance)
                    err=MIN_F_TOLERANCE;
            if(err==LS_S)
                if(istp==max_iter)
                    err=MIN_F_MAX_ITER;
            
            atoms->vectors[f_n].ret(f);
            
            for(int i=0;i<x_dim*atoms->natms;i++)
                f[i]=0.0;
            
            if(thermo->test_prev_step() || err)
            {
                forcefield->force_calc(1,energy_stress);
                thermo->update(pe_idx,energy_stress[0]);
                thermo->update(stress_idx,6,&energy_stress[1]);
                curr_energy=energy_stress[0];
            }
            else
                forcefield->force_calc(0,&curr_energy);
            
            if(err)
                continue;
            
            if(m_it)
            {
                if(k_it==m_it)
                {
                    s_list_tmp=s_list[0];
                    y_list_tmp=y_list[0];
                    for(int i=0;i<m_it-1;i++)
                    {
                        s_list[i]=s_list[i+1];
                        y_list[i]=y_list[i+1];
                        rho[i]=rho[i+1];
                    }
                    s_list[m_it-1]=s_list_tmp;
                    y_list[m_it-1]=y_list_tmp;
                    
                }
                else
                    k_it++;
                
                
                atoms->vectors[0].ret(x);
                atoms->vectors[f_n].ret(f);
                atoms->vectors[h_n].ret(h);
                atoms->vectors[x_prev_n].ret(x_0);
                atoms->vectors[f_prev_n].ret(f_0);
                atoms->vectors[s_list[k_it-1]].ret(tmp_vec0);
                atoms->vectors[y_list[k_it-1]].ret(tmp_vec1);
                
                inner=0.0;
                inner_tmp=0.0;
                for(int i=0;i<x_dim*atoms->natms;i++)
                {
                    tmp_vec0[i]=alpha_m*h[i];
                    tmp_vec1[i]=f[i]-f_0[i];
                    inner+=tmp_vec0[i]*tmp_vec1[i];
                    inner_tmp+=tmp_vec1[i]*tmp_vec1[i];
                }
                tot_inner=0.0;
                MPI_Allreduce(&inner,&tot_inner,1,MPI_TYPE0,MPI_SUM,world);
                tot_inner_tmp=0.0;
                MPI_Allreduce(&inner_tmp,&tot_inner_tmp,1,MPI_TYPE0,MPI_SUM,world);
                
                ratio=-tot_inner/tot_inner_tmp;
                rho[k_it-1]=1.0/tot_inner;
            }
            else
            {
                atoms->vectors[f_n].ret(f);
                atoms->vectors[h_n].ret(h);
                atoms->vectors[f_prev_n].ret(f_0);
                inner=0.0;
                inner_tmp=0.0;
                for(int i=0;i<x_dim*atoms->natms;i++)
                {
                    inner+=(alpha_m*h[i])*(f[i]-f_0[i]);
                    inner_tmp+=(f[i]-f_0[i])*(f[i]-f_0[i]);
                }
                tot_inner=0.0;
                MPI_Allreduce(&inner,&tot_inner,1,MPI_TYPE0,MPI_SUM,world);
                tot_inner_tmp=0.0;
                MPI_Allreduce(&inner_tmp,&tot_inner_tmp,1,MPI_TYPE0,MPI_SUM,world);
                
                ratio=-tot_inner/tot_inner_tmp;
            }
            
            istp++;
            step_no++;
        }
    }
    
    delete [] energy_stress;
}
/*--------------------------------------------
 fin after a run
 --------------------------------------------*/
void Min_LBFGS::fin()
{
    forcefield->fin();
    neighbor->fin();
    atoms->del(h_n);
    atoms->del(f_prev_n);
    atoms->del(x_prev_n);
    
    delete line_search;
    delete vecs_comm;
    
    int tmp;
    
    for(int i=0;i<m_it;i++)
    {
        for(int j=i+1;j<m_it;j++)
        {
            if(s_list[j]<s_list[i])
            {
                tmp=s_list[j];
                s_list[j]=s_list[i];
                s_list[i]=tmp;
            }
            
            if(y_list[j]<y_list[i])
            {
                tmp=y_list[j];
                y_list[j]=y_list[i];
                y_list[i]=tmp;
            }
        }
    }
    
    
    if(m_it)
    {
        for(int i=m_it-1;i>-1;i--)
        {
            atoms->del(y_list[i]);
            atoms->del(s_list[i]);
        }
        
        delete [] s_list;
        delete [] y_list;
        
        delete [] rho;
        delete [] alpha;
    }

 
    if(chng_box)
    {
        for(int i=0;i<m_it;i++)
        {
            for(int j=0;j<dim;j++)
            {
                delete [] H_y[i][j];
                delete [] H_s[i][j];
            }
            
            delete [] H_y[i];
            delete [] H_s[i];
        }
        delete [] H_y;
        delete [] H_s;
        
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
        
        if(m_it)
            delete [] H_s_y_list;
    }
    
    if(write!=NULL)
        write->fin();
    thermo->fin();
    errors();
    atoms->x2s(atoms->natms);
}

