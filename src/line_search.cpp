

#include "line_search.h"
#include "ff.h"
using namespace std;
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch::LineSearch(MAPP* mapp,VecLst*vec_list)
: InitPtrs(mapp)
{
    x_n=atoms->find("x");
    x_prev_n=atoms->find("x_prev");
    f_n=atoms->find("f");
    dim=atoms->dimension;
    x_dim=atoms->vectors[0].dim;
    vecs_comm=vec_list;
    dof_n=atoms->find_exist("dof");
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearch::~LineSearch()
{

}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
TYPE0 LineSearch::inner_f_h()
{
    TYPE0* f;
    TYPE0* h;
    atoms->vectors[f_n].ret(f);
    atoms->vectors[h_n].ret(h);
    TYPE0 inner=0.0;
    TYPE0 inner_tot=0.0;
    for(int i=0;i<x_dim*atoms->natms;i++)
        inner+=f[i]*h[i];
    
    MPI_Allreduce(&inner,&inner_tot,1,MPI_TYPE0,MPI_SUM,world);
    
    return inner_tot;
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
TYPE0 LineSearch::inner_f_h_s()
{
    TYPE0* f;
    TYPE0* h;
    atoms->vectors[f_n].ret(f);
    atoms->vectors[h_n].ret(h);
    TYPE0 inner=0.0;
    TYPE0 inner_tot=0.0;
    for(int i=0;i<x_dim*atoms->natms;i++)
        inner+=f[i]*h[i];
    
    MPI_Allreduce(&inner,&inner_tot,1,MPI_TYPE0,MPI_SUM,world);
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            inner_tot+=h_H[i][j]*f_H[i][j];
    
    return inner_tot;
}

/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
void LineSearch::normalize_h()
{
    TYPE0* h;
    atoms->vectors[h_n].ret(h);
    TYPE0 inner=0.0;
    TYPE0 inner_tot=0.0;
    for(int i=0;i<x_dim*atoms->natms;i++)
        inner+=h[i]*h[i];
    MPI_Allreduce(&inner,&inner_tot,1,MPI_TYPE0,MPI_SUM,world);
        
    inner_tot=d_max/sqrt(inner_tot);
    
    for(int i=0;i<x_dim*atoms->natms;i++)
        h[i]*=inner_tot;
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
void LineSearch::normalize_h_s()
{
    TYPE0* h;
    atoms->vectors[h_n].ret(h);
    TYPE0 inner=0.0;
    TYPE0 inner_tot=0.0;
    for(int i=0;i<x_dim*atoms->natms;i++)
        inner+=h[i]*h[i];
    MPI_Allreduce(&inner,&inner_tot,1,MPI_TYPE0,MPI_SUM,world);
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            inner_tot+=h_H[i][j]*h_H[i][j];
    
    inner_tot=s_max/sqrt(inner_tot);
    
    for(int i=0;i<x_dim*atoms->natms;i++)
        h[i]*=inner_tot;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            h_H[i][j]*=inner_tot;
}
/*--------------------------------------------
 inner product of f and h
 --------------------------------------------*/
TYPE0 LineSearch::energy(TYPE0 alpha)
{
    TYPE0* x;
    TYPE0* x_prev;
    TYPE0* h;
    char* dof=NULL;
    
    atoms->vectors[x_n].ret(x);
    atoms->vectors[x_prev_n].ret(x_prev);
    atoms->vectors[h_n].ret(h);
    if(dof_n>-1)
        atoms->vectors[dof_n].ret(dof);

    if(chng_box)
    {
        TYPE0** H=atoms->H;
        
        for(int i=0;i<dim;i++)
        {
            for(int j=0;j<dim;j++)
            {
                H[i][j]=H_prev[i][j]+alpha*h_H[i][j];
                N[i][j]=0.0;
                for(int k=0;k<dim;k++)
                {
                    N[i][j]+=alpha*h_H[i][k]*B_prev[k][j];
                }
            }
        }
        for(int i=0;i<dim;i++)
            N[i][i]++;

        
        if(dim==3)
            M3INV_TRI_LOWER(atoms->H,atoms->B);
        else
            atoms->invert_lower_triangle(atoms->H,atoms->B,dim);
       
        int icomp=0;
        for(int i=0;i<atoms->natms;i++)
        {
            for(int j=0;j<dim;j++)
            {
                x[icomp+j]=0.0;
                for(int k=j;k<dim;k++)
                    x[icomp+j]+=N[k][j]*x_prev[icomp+k]+alpha*H[k][j]*h[icomp+k];
            }
            icomp+=x_dim;
        }
        
        if(dof_n>-1)
            for(int i=0;i<(atoms->natms)*x_dim;i++)
                if(dof[i]==1) x[i]=x_prev[i];
        
        atoms->update_0(1,1,vecs_comm);
        
        return forcefield->energy_calc();
    }
    else
    {
        
        for(int i=0;i<x_dim*atoms->natms;i++)
            x[i]=x_prev[i]+alpha*h[i];
        
        atoms->update_0(0,1,vecs_comm);
        
        return forcefield->energy_calc();
    }
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
LineSearch_BackTrack::LineSearch_BackTrack(MAPP* mapp,VecLst*vec_list)
:LineSearch(mapp,vec_list)
{
    c=0.4;
    rho=0.5;
    alpha_max=1.0;
    alpha_min=0.0;
    d_max=1.0;
    s_max=0.1;
    CREATE2D(N,dim,dim);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
LineSearch_BackTrack::~LineSearch_BackTrack()
{
    for(int i=0;i<dim;i++)
        delete [] N[i];
    delete [] N;
}
/*--------------------------------------------
 minimize line
 --------------------------------------------*/
int LineSearch_BackTrack::line_min(TYPE0& nrgy,TYPE0& alph)
{
    TYPE0 inner;
    TYPE0 alpha_m;
    TYPE0 max_h;
    TYPE0 max_h_tot;
    TYPE0 current_energy,ideal_energy;
    TYPE0* h;
    
    
    if(chng_box)
    {
        alph=0.0;
        normalize_h_s();
        inner=inner_f_h_s();
        
        if(inner<=0)
            return LS_F_DOWNHILL;
        
        
        atoms->vectors[h_n].ret(h);
        
        max_h=0.0;
        max_h_tot=0.0;
        for(int i=0;i<atoms->natms*x_dim;i++)
            max_h=MAX(max_h,fabs(h[i]));
        MPI_Allreduce(&max_h,&max_h_tot,1,MPI_TYPE0,MPI_MAX,world);
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                max_h_tot=MAX(max_h_tot,fabs(h_H[i][j]));
        
        if(max_h_tot==0)
            return LS_F_GRAD0;
        
        alpha_m=MIN(alpha_max,max_h_tot);
        
        if(alpha_m<=alpha_min)
            return LS_F_ALPHAMIN;
        
        while (1)
        {
            
            ideal_energy=nrgy-alpha_m*c*inner;
            current_energy=energy(alpha_m);
            if(current_energy<=ideal_energy)
            {
                nrgy=current_energy;
                alph=alpha_m;
                return LS_S;
            }
            alpha_m*=rho;
            
            if(alpha_m<=alpha_min)
            {
                nrgy=energy(0);
                return LS_F_ALPHAMIN;
            }
        }
    }
    else
    {
        alph=0.0;
        normalize_h();
        inner=inner_f_h();
        
        if(inner<=0)
            return LS_F_DOWNHILL;
        
        
        atoms->vectors[h_n].ret(h);
        
        max_h=0.0;
        max_h_tot=0.0;
        for(int i=0;i<atoms->natms*x_dim;i++)
            max_h=MAX(max_h,fabs(h[i]));
        MPI_Allreduce(&max_h,&max_h_tot,1,MPI_TYPE0,MPI_MAX,world);
        if(max_h_tot==0)
            return LS_F_GRAD0;
        
        alpha_m=MIN(alpha_max,max_h_tot);

        if(alpha_m<=alpha_min)
            return LS_F_ALPHAMIN;
        
        while (1)
        {
            ideal_energy=nrgy-alpha_m*c*inner;
            current_energy=energy(alpha_m);
            if(current_energy<=ideal_energy)
            {
                nrgy=current_energy;
                alph=alpha_m;
                return LS_S;
            }
            alpha_m*=rho;
            
            if(alpha_m<=alpha_min)
            {
                nrgy=energy(0);
                return LS_F_ALPHAMIN;
            }
        }
    }
    
}
