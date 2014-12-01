

#include "ff_eam_dmd_full.h"
#include "neighbor.h"
#include "atom_types.h"
using namespace MAPP_NS;
#define PI_IN_SQ 0.564189583547756286948079451561
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam_dmd_full::
ForceField_eam_dmd_full(MAPP* mapp) : ForceField(mapp)
{
    no_i=0;
    allocated=0;
    eam_mode=NOT_SET;
    max_pairs=0;
    neigh_lst_sz_sz=0;
    
    
    no_types=atom_types->no_types;
    if(atoms->vectors[0].dim!=3+no_types)
        error->abort("the dimension of x"
        " vector should be 3 + no of types");
    
    
    
    CREATE1D(nrgy_strss,7);
    CREATE1D(cut_sk_sq,no_types*(no_types+1)/2);
    CREATE1D(c_0,no_types);
    
    c_d_n=atoms->find_exist("c_d");
    
    if(c_d_n<0)
    {
        int tmp0=atoms->find("c");
        tmp0=atoms->vectors[tmp0].dim;
        c_d_n=atoms->add<TYPE0>(0,tmp0,"c_d");
    }
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam_dmd_full::~ForceField_eam_dmd_full()
{
    
    if(no_types)
    {
        delete [] cut_sk_sq;
        delete [] c_0;
    }
    
    delete [] nrgy_strss;
    
    if(allocated) clean_up();
    
    if(no_i)
    {
        delete [] xi;
        delete [] wi_0;
        delete [] wi_1;
        delete [] wi_2;
    }
    
    if(max_pairs)
    {
        delete [] rho;
        delete [] drho_dr;
        delete [] drho_dalpha;
        delete [] phi;
        delete [] dphi_dr;
        delete [] dphi_dalpha;
    }
    
    if(neigh_lst_sz_sz)
    {
        for(int i=0;i<neigh_lst_sz_sz;i++)
            if(neigh_lst_sz[i])
                delete [] neigh_lst[i];
        
        delete [] neigh_lst_sz;
        delete [] neigh_lst;
    }
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam_dmd_full::
force_calc(int st_clc,TYPE0* en_st)
{
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    TYPE0* E;
    atoms->vectors[E_n].ret(E);
    TYPE0* dE;
    atoms->vectors[dE_n].ret(dE);
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
       
    int iatm,jatm;
    
    int icomp,jcomp,phi_start,rho_start;
    TYPE0 dx0,dx1,dx2,rsq,z2p,z2;
    TYPE0 alpha;
    TYPE0 r=0.0,r_inv=0.0;
    int m;
    TYPE0* coef;
    TYPE0 alpha_sq,inv_alph_sq,tmp0,tmp1,tmp2;
    TYPE0 upper,lower;
    TYPE0 fpair,apair;
    
    TYPE0 p,p2,p3,p4;
    TYPE0 f_0,f_1,f_2,g_0,g_1,g_2;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    if(max_pairs<neighbor->no_pairs)
    {
        if(max_pairs)
        {
            delete [] rho;
            delete [] drho_dr;
            delete [] drho_dalpha;
            delete [] phi;
            delete [] dphi_dr;
            delete [] dphi_dalpha;
        }
        max_pairs=neighbor->no_pairs;
        int no_0=max_pairs*no_types*no_types;
        int no_1=max_pairs*no_types*(no_types+1)/2;
        CREATE1D(rho,no_0);
        CREATE1D(drho_dr,no_0);
        CREATE1D(drho_dalpha,no_0);
        for(int i=0;i<no_0;i++)
            rho[i]=drho_dr[i]=drho_dalpha[i]=0.0;
        CREATE1D(phi,no_1);
        CREATE1D(dphi_dr,no_1);
        CREATE1D(dphi_dalpha,no_1);
        for(int i=0;i<no_1;i++)
            phi[i]=dphi_dr[i]=dphi_dalpha[i]=0.0;
    }
    
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++) E[i]=0.0;
    for(int i=0;i<natms*(3+no_types);i++) f[i]=0.0;
    
    phi_start=0;
    rho_start=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+no_types)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int itype=0;itype<no_types; itype++)
                {
                    for(int jtype=0;jtype<no_types; jtype++)
                    {
                        alpha=(x[icomp+3+itype]*x[jcomp+3+jtype])/(x[icomp+3+itype]+x[jcomp+3+jtype]);

                        if(alpha_min<alpha && alpha<alpha_max)
                        {
                            alpha_sq=sqrt(alpha);
                            inv_alph_sq=1.0/alpha_sq;
                            upper=(r+rc)*alpha_sq;
                            lower=(r-rc)*alpha_sq;
                            
                            if(lower<xi[no_i-1])
                            {
                                f_0=f_1=f_2=g_0=g_1=g_2=0.0;
                                for(int i=0;i<no_i;i++)
                                {
                                    if(xi[i]>lower && xi[i]<upper)
                                    {
                                        tmp0=r-xi[i]*inv_alph_sq;
                                        
                                        p=fabs(tmp0)*dr_inv;
                                        m=static_cast<int>(p);
                                        m=MIN(m,nr-2);
                                        
                                        p-=m;
                                        p=MIN(p,1.0);
                                        p2=p*p;
                                        p3=p2*p;
                                        p4=p3*p;
                                        
                                        coef=rho_arr[type2rho[itype][jtype]][m];
                                        tmp2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                        tmp2*=tmp0;
                                        g_0+=wi_0[i]*tmp2;
                                        g_1+=wi_1[i]*tmp2;
                                        g_2+=wi_2[i]*tmp2;
                                        
                                        if(itype>=jtype)
                                        {
                                            coef=phi_r_arr[type2phi[itype][jtype]][m];
                                            tmp1=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            if(tmp0<0.0)
                                                tmp1*=-1.0;
                                            f_0+=wi_0[i]*tmp1;
                                            f_1+=wi_1[i]*tmp1;
                                            f_2+=wi_2[i]*tmp1;
                                        }
                                    }
                                }

                                
                                if(itype>=jtype)
                                {
                                    f_0*=PI_IN_SQ*r_inv;
                                    f_1*=PI_IN_SQ*r_inv;
                                    f_2*=PI_IN_SQ*r_inv;
                                    phi[phi_start+type2phi[itype][jtype]]=f_0;
                                    dphi_dr[phi_start+type2phi[itype][jtype]]=-f_0*r_inv-2.0*alpha_sq*f_1;
                                    dphi_dalpha[phi_start+type2phi[itype][jtype]]=(0.5*f_0-f_2)/alpha;
                                }
                                
                                g_0*=PI_IN_SQ*r_inv;
                                g_1*=PI_IN_SQ*r_inv;
                                g_2*=PI_IN_SQ*r_inv;
                                rho[rho_start+type2rho[itype][jtype]]=g_0;
                                drho_dr[rho_start+type2rho[itype][jtype]]=-g_0*r_inv-2.0*alpha_sq*g_1;
                                drho_dalpha[rho_start+type2rho[itype][jtype]]=(0.5*g_0-g_2)/alpha;

                                E[iatm*no_types+itype]+=g_0*c[jatm*no_types+jtype];
                                
                                if (jatm<natms)
                                {
                                    E[jatm*no_types+jtype]+=g_0*c[iatm*no_types+itype];
                                    if(itype>=jtype)
                                        nrgy_strss[0]+=0.5*(c[jatm*no_types+jtype]*c[iatm*no_types+itype]
                                        +c[jatm*no_types+itype]*c[iatm*no_types+jtype])*f_0;
                                }
                                else
                                {
                                    if(itype>=jtype)
                                        nrgy_strss[0]+=0.25*(c[jatm*no_types+jtype]*c[iatm*no_types+itype]
                                        +c[jatm*no_types+itype]*c[iatm*no_types+jtype])*f_0;
                                }
                            }
                        }
                        else if (alpha_max<=alpha)
                        {
                            if(rsq<cut_sq)
                            {
                                r=sqrt(rsq);
                                p=r*dr_inv;
                                m=static_cast<int>(p);
                                m=MIN(m,nr-2);
                                p-=m;
                                p=MIN(p,1.0);
                                p2=p*p;
                                p3=p2*p;
                                p4=p3*p;
                                
                                if(itype>=jtype)
                                {
                                    coef=phi_r_arr[type2phi[itype][jtype]][m];
                                    z2p=coef[6]*p2+coef[5]*p+coef[4];
                                    z2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                    phi[phi_start+type2phi[itype][jtype]]=z2*r_inv;
                                    dphi_dalpha[phi_start+type2phi[itype][jtype]]=0.0;
                                    dphi_dr[phi_start+type2phi[itype][jtype]]=(z2p-z2*r_inv)*r_inv;
                                }
                                
                                coef=rho_arr[type2rho[itype][jtype]][m];
                                rho[rho_start+type2rho[itype][jtype]]=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                drho_dalpha[rho_start+type2rho[itype][jtype]]=0.0;
                                drho_dr[rho_start+type2rho[itype][jtype]]=coef[6]*p2+coef[5]*p+coef[4];
                                
                                E[iatm*no_types+itype]+=rho[rho_start+type2rho[itype][jtype]]*c[jatm*no_types+jtype];
                                
                                if (jatm<natms)
                                {
                                    E[jatm*no_types+jtype]+=rho[rho_start+type2rho[itype][jtype]]*c[iatm*no_types+itype];
                                    if(itype>=jtype)
                                    {
                                        nrgy_strss[0]+=0.5*(c[jatm*no_types+jtype]*c[iatm*no_types+itype]
                                        +c[jatm*no_types+itype]*c[iatm*no_types+jtype])*phi[phi_start+type2phi[itype][jtype]];
                                    }
                                }
                                else
                                {
                                    if(itype>=jtype)
                                    {
                                        nrgy_strss[0]+=0.25*(c[jatm*no_types+jtype]*c[iatm*no_types+itype]
                                        +c[jatm*no_types+itype]*c[iatm*no_types+jtype])*phi[phi_start+type2phi[itype][jtype]];
                                    }
                                }
                            }
                        }
                    }
                    
                }
            }

            
            phi_start+=no_types*(no_types+1)/2;
            rho_start+=no_types*no_types;
        }
        
        
        for(int itype=0;itype<no_types; itype++)
        {
            p=E[iatm*no_types+itype]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            p4=p3*p;
            coef=F_arr[itype][m];
            
            tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            tmp1=coef[8]*p3+coef[7]*p2+coef[6]*p+coef[5];
            
            if(E[iatm*no_types+itype]>rho_max)
                tmp0+=tmp1*(E[iatm*no_types+itype]-rho_max);
            
            E[iatm*no_types+itype]=tmp0;
            dE[iatm*no_types+itype]=tmp1;
            
            nrgy_strss[0]+=c[iatm*no_types+itype]*E[iatm*no_types+itype];

            nrgy_strss[0]+=kbT*calc_ent(c[iatm*no_types+itype]);
            
            nrgy_strss[0]+=1.5*kbT*c[iatm*no_types+itype]*log(x[icomp+3+itype])
            +c[iatm*no_types+itype]*c_0[itype];
            
            f[icomp+3+itype]-=1.5*kbT*c[iatm*no_types+itype]/x[icomp+3+itype];
        }
        
    }
    
    atoms->update(dE_n);
    
    phi_start=0;
    rho_start=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+no_types)*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            r_inv=1.0/sqrt(dx0*dx0+dx1*dx1+dx2*dx2);
            
            for(int itype=0;itype<no_types; itype++)
            {
                for(int jtype=0;jtype<no_types; jtype++)
                {
                    fpair=-(drho_dr[rho_start+type2rho[jtype][itype]]*dE[iatm*no_types+itype]
                    +drho_dr[rho_start+type2rho[itype][jtype]]*dE[jatm*no_types+jtype]
                    +dphi_dr[phi_start+type2phi[itype][jtype]])
                    *(c[iatm*no_types+itype]*c[jatm*no_types+jtype])*r_inv;
                    
                    apair=-(drho_dalpha[rho_start+type2rho[jtype][itype]]*dE[iatm*no_types+itype]
                    +drho_dalpha[rho_start+type2rho[itype][jtype]]*dE[jatm*no_types+jtype]
                    +dphi_dalpha[phi_start+type2phi[itype][jtype]])
                    *c[iatm*no_types+itype]*c[jatm*no_types+jtype]/
                    ((x[icomp+3+itype]+x[jcomp+3+jtype])
                     *(x[icomp+3+itype]+x[jcomp+3+jtype]));
                    
                    
                    if(apair!=0.0 || fpair!=0.0)
                    {
                        f[icomp]+=dx0*fpair;
                        f[icomp+1]+=dx1*fpair;
                        f[icomp+2]+=dx2*fpair;
                        f[icomp+3+itype]+=apair*x[jcomp+3+jtype]*x[jcomp+3+jtype];
                        
                        if (jatm<natms)
                        {
                            f[jcomp]-=dx0*fpair;
                            f[jcomp+1]-=dx1*fpair;
                            f[jcomp+2]-=dx2*fpair;
                            f[jcomp+3+jtype]+=apair*x[icomp+3+itype]*x[icomp+3+itype];
                            
                            if (st_clc)
                            {
                                nrgy_strss[1]+=fpair*dx0*dx0;
                                nrgy_strss[2]+=fpair*dx1*dx1;
                                nrgy_strss[3]+=fpair*dx2*dx2;
                                nrgy_strss[4]+=fpair*dx1*dx2;
                                nrgy_strss[5]+=fpair*dx2*dx0;
                                nrgy_strss[6]+=fpair*dx0*dx1;
                            }
                        }
                        else
                        {
                            if (st_clc)
                            {
                                nrgy_strss[1]+=0.5*fpair*dx0*dx0;
                                nrgy_strss[2]+=0.5*fpair*dx1*dx1;
                                nrgy_strss[3]+=0.5*fpair*dx2*dx2;
                                nrgy_strss[4]+=0.5*fpair*dx1*dx2;
                                nrgy_strss[5]+=0.5*fpair*dx2*dx0;
                                nrgy_strss[6]+=0.5*fpair*dx0*dx1;
                            }
                        }
                            
                        
                    }
                    
                }
            }
            

            phi_start+=no_types*(no_types+1)/2;
            rho_start+=no_types*no_types;
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
 energy calculation
 --------------------------------------------*/
TYPE0 ForceField_eam_dmd_full::energy_calc()
{
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    TYPE0* E;
    atoms->vectors[E_n].ret(E);
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    
    int iatm,jatm;
    
    int icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,z2;
    TYPE0 alpha;
    TYPE0 r,r_inv=0.0;
    int m;
    TYPE0* coef;
    TYPE0 alpha_sq,inv_alph_sq,tmp0,tmp1,tmp2;
    TYPE0 upper,lower;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++) E[i]=0.0;
    TYPE0 en=0.0;
    TYPE0 en_tot=0.0;
    
    
    TYPE0 f_0,g_0;
    TYPE0 p,p2,p3,p4;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+no_types)*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            if(rsq<cut_sq_mod)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                
                for(int itype=0;itype<no_types; itype++)
                {
                    for(int jtype=itype;jtype<no_types; jtype++)
                    {
                        alpha=(x[icomp+3+itype]*x[jcomp+3+jtype])/(x[icomp+3+itype]+x[jcomp+3+jtype]);

                        if(alpha_min<alpha && alpha<alpha_max)
                        {
                            alpha_sq=sqrt(alpha);
                            inv_alph_sq=1.0/alpha_sq;
                            upper=(r+rc)*alpha_sq;
                            lower=(r-rc)*alpha_sq;
                            
                            f_0=g_0=0.0;
                            if(lower<xi[no_i-1])
                            {
                                
                                for(int i=0;i<no_i;i++)
                                {
                                    if(xi[i]>lower && xi[i]<upper)
                                    {
                                        tmp0=r-xi[i]*inv_alph_sq;
                                        
                                        p=fabs(tmp0)*dr_inv;
                                        m=static_cast<int>(p);
                                        m=MIN(m,nr-2);
                                        
                                        p-=m;
                                        p=MIN(p,1.0);
                                        p2=p*p;
                                        p3=p2*p;
                                        p4=p3*p;
                                        
                                        if(itype>=jtype)
                                        {
                                            coef=phi_r_arr[type2phi[itype][jtype]][m];
                                            tmp1=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                            if(tmp0<0.0)
                                                tmp1*=-1.0;
                                            f_0+=wi_0[i]*tmp1;
                                        }
                                        coef=rho_arr[type2rho[itype][jtype]][m];
                                        tmp2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                        tmp2*=tmp0;
                                        g_0+=wi_0[i]*tmp2;
                                    }
                                }
                                if(itype>=jtype)
                                    f_0*=PI_IN_SQ*r_inv;
                                
                                g_0*=PI_IN_SQ*r_inv;

                                E[iatm*no_types+itype]+=g_0*c[jatm*no_types+jtype];
                                
                                if (jatm<natms)
                                {
                                    E[jatm*no_types+jtype]+=g_0*c[iatm*no_types+itype];
                                    if(itype>=jtype)
                                    {
                                        en+=0.5*(c[jatm*no_types+jtype]*c[iatm*no_types+itype]
                                        +c[jatm*no_types+itype]*c[iatm*no_types+jtype])*f_0;
                                    }
                                }
                                else
                                {
                                    if(itype>=jtype)
                                    {
                                        en+=0.25*(c[jatm*no_types+jtype]*c[iatm*no_types+itype]
                                        +c[jatm*no_types+itype]*c[iatm*no_types+jtype])*f_0;
                                    }
                                    
                                }

                            }
                            
                        }
                        else if (alpha_max<=alpha)
                        {
                            f_0=0.0;
                            if(rsq<cut_sq)
                            {
                                r=sqrt(rsq);
                                p=r*dr_inv;
                                m=static_cast<int>(p);
                                m=MIN(m,nr-2);
                                p-=m;
                                p=MIN(p,1.0);
                                p2=p*p;
                                p3=p2*p;
                                p4=p3*p;
                                
                                
                                if(itype>=jtype)
                                {
                                    coef=phi_r_arr[type2phi[itype][jtype]][m];
                                    z2=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                    f_0=z2*r_inv;
                                }

                                coef=rho_arr[type2rho[itype][jtype]][m];
                                g_0=coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
                                
                                E[iatm*no_types+itype]+=g_0*c[jatm*no_types+jtype];
                                
                                if (jatm<natms)
                                {
                                    E[jatm*no_types+jtype]+=g_0*c[iatm*no_types+itype];
                                    if(itype>=jtype)
                                    {
                                        en+=0.5*(c[jatm*no_types+jtype]*c[iatm*no_types+itype]
                                        +c[jatm*no_types+itype]*c[iatm*no_types+jtype])*f_0;
                                    }
                                }
                                else
                                {
                                    if(itype>=jtype)
                                    {
                                        en+=0.25*(c[jatm*no_types+jtype]*c[iatm*no_types+itype]
                                        +c[jatm*no_types+itype]*c[iatm*no_types+jtype])*f_0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        
        for(int itype=0;itype<no_types; itype++)
        {
            p=E[iatm*no_types+itype]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            p4=p3*p;
            coef=F_arr[itype][m];
            
            tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            tmp1=coef[8]*p3+coef[7]*p2+coef[6]*p+coef[5];
            
            if(E[iatm*no_types+itype]>rho_max)
                tmp0+=tmp1*(E[iatm*no_types+itype]-rho_max);
            
            E[iatm*no_types+itype]=tmp0;
            
            en+=c[iatm*no_types+itype]*E[iatm*no_types+itype];
            
            /*
            en+=kbT*(c[iatm*no_types+itype]*log(c[iatm*no_types+itype])
            +(1.0-c[iatm*no_types+itype])*log(1.0-c[iatm*no_types+itype]));
            */
            
            en+=kbT*calc_ent(c[iatm*no_types+itype]);
            
            en+=1.5*kbT*c[iatm*no_types+itype]*log(x[icomp+3+itype])
            +c[iatm*no_types+itype]*c_0[itype];
        }
        
    }
    
    
    
    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    
    return en_tot;
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_eam_dmd_full::init()
{
    TYPE0 skin=atoms->skin;
    TYPE0 ph_cut=0.0;
    for (int i=0;i<no_types*(no_types+1)/2;i++)
        cut_sk_sq[i]=cut_sq_mod+(skin)*(skin)
        +2.0*sqrt(cut_sq_mod)*(skin);
    
    ph_cut=sqrt(cut_sq_mod);
    
    atoms->set_ph(ph_cut);
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    c_n=atoms->find("c");
    
    neighbor->pair_wise=1;
    
    E_n=atoms->add<TYPE0>(1,no_types,"E");
    dE_n=atoms->add<TYPE0>(1,no_types,"dE");
    ddE_n=atoms->add<TYPE0>(1,no_types,"ddE");
    n_n=atoms->add<TYPE0>(1,no_types,"n");
    s_n=atoms->add<TYPE0>(1,no_types,"s");
    t_n=atoms->add<TYPE0>(1,no_types,"t");
    v_n=atoms->add<TYPE0>(1,no_types,"v");
    crd_n=atoms->add<int>(1,1,"crd");
    
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_eam_dmd_full::fin()
{
    atoms->del(crd_n);
    atoms->del(v_n);
    atoms->del(t_n);
    atoms->del(s_n);
    atoms->del(n_n);
    atoms->del(ddE_n);
    atoms->del(dE_n);
    atoms->del(E_n);
    
    if(max_pairs)
    {
        delete [] rho;
        delete [] drho_dr;
        delete [] drho_dalpha;
        delete [] phi;
        delete [] dphi_dr;
        delete [] dphi_dalpha;
        max_pairs=0;
    }

}
/*--------------------------------------------
 ff_coef no_i alpha_min alpha_max T kB hbar
 SetFL Ni_u3.eam
 hbar 6.5821192815e-16 eVs
 kB 8.617332478e-5 eV/K
 T 300.0 K
 mass conversion from amu to eVs^2/A^2:
 1.0364269184093291236e-28
 --------------------------------------------*/
void ForceField_eam_dmd_full::coef(int narg,char** arg)
{
    TYPE0 kb,T,hbar;
    if (narg<10)
        error->abort("wrong coeff command "
                     "for eam Force Field");
    rsq_crd=atof(arg[1]);
    rsq_crd*=rsq_crd;
    set_weight_abs(atoi(arg[2]));
    alpha_min=atof(arg[3]);
    alpha_max=atof(arg[4]);
    T=atof(arg[5]);
    kb=atof(arg[6]);
    hbar=atof(arg[7]);
    
    clean_up();
    if(strcmp(arg[8],"FS")==0)
    {
        eam_mode=FINNIS_FL;
        set_funcfl(narg-9,&arg[9]);
    }
    else if(strcmp(arg[8],"SetFL")==0)
    {
        eam_mode=SET_FL;
        set_setfl(narg-9,&arg[9]);
    }
    else if(strcmp(arg[8],"FuncFL")==0)
    {
        eam_mode=FUNC_FL;
        set_funcfl(narg-9,&arg[9]);
    }
    else
        error->abort("wrong coeff command "
                     "for eam Force Field");
    set_arrays();
    
    TYPE0 mass;
    TYPE0 deb_l;

    for(int i=0;i<no_types;i++)
    {
        mass=atom_types->mass[i];
        mass*=1.0364269184093291236e-28;
        deb_l=hbar*hbar*2.0/(mass*kb*T);
        c_0[i]=1.5*kb*T*(log(deb_l)-1.0);
    }

    kbT=kb*T;
    beta=1.0/kbT;
    
    if(alpha_min==0.0)
        error->abort("minimum alpha cannot be zero");
    
    rc=(static_cast<TYPE0>(nr)-1.0)*dr;
    rho_max=(static_cast<TYPE0>(nrho)-1.0)*drho;
    cut_sq=rc*rc;
    mod_rc=rc+xi[no_i-1]/sqrt(alpha_min);
    cut_sq_mod=mod_rc*mod_rc;
    

    
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam_dmd_full::set_setfl(int no_files
                                   ,char** file_names)
{
    if(no_files!=no_types)
        error->abort("no of types and number "
         "of files should be equal in eam");
    
    TYPE0* drs;
    TYPE0* drhos;
    int* nrs;
    int* nrhos;
    CREATE1D(drs,no_types);
    CREATE1D(drhos,no_types);
    CREATE1D(nrs,no_types);
    CREATE1D(nrhos,no_types);
    
    TYPE0** tmp_F;
    TYPE0** tmp_rho;
    TYPE0** tmp_zi;
    CREATE1D(tmp_F,no_types);
    CREATE1D(tmp_rho,no_types);
    CREATE1D(tmp_zi,no_types);
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    int narg;
    char** arg;
    int nlines;
    int err;
    int no;
    for(int ityp=0;ityp<no_types;ityp++)
    {
        fp=NULL;
        if(atoms->my_p_no==0)
        {
            fp=fopen(file_names[ityp],"r");
            if(fp==NULL)
                error->abort("file %s not found",file_names[ityp]);
        }
        
        for(int i=0;i<3;i++)
            if(line_read(fp,line)==-1)
                error->abort("eam potential file ended immaturely");
        
        narg=mapp->parse_line(line,arg);
        if(narg!=5)
            error->abort("wrong line in eam file: %s",line);
        
        nrhos[ityp]=atoi(arg[0]);
        nrs[ityp]=atoi(arg[2]);
        drhos[ityp]=atof(arg[1]);
        drs[ityp]=atof(arg[3]);
        
        
        CREATE1D(tmp_F[ityp],nrhos[ityp]+1);
        CREATE1D(tmp_rho[ityp],nrs[ityp]+1);
        CREATE1D(tmp_zi[ityp],nrs[ityp]+1);
        
        for(int i=0;i<5;i++)
            delete [] arg[i];
        delete [] arg;
        
        if(nrhos[ityp]%5!=0 || nrs[ityp]%5!=0)
            
            error->abort("nro and nrho remainder by 5 should be zero");
        
        nlines=nrhos[ityp]/5;
        no=1;
        for(int i=0;i<nlines;i++)
        {
            err=line_read(fp,line);
            if(err==-1)
                error->abort("eam potential file ended immaturely");
            
            if(sscanf(line,"%lf %lf %lf %lf %lf"
                      ,&tmp_F[ityp][no],&tmp_F[ityp][no+1]
                      ,&tmp_F[ityp][no+2],&tmp_F[ityp][no+3]
                      ,&tmp_F[ityp][no+4])!=5)
                error->abort("wrong line in eam file: %s",line);
            no+=5;
            
        }
        
        nlines=nrs[ityp]/5;
        no=1;
        for(int i=0;i<nlines;i++)
        {
            err=line_read(fp,line);
            if(err==-1)
                error->abort("eam potential file ended immaturely");
            
            if(sscanf(line,"%lf %lf %lf %lf %lf"
                      ,&tmp_zi[ityp][no],&tmp_zi[ityp][no+1]
                      ,&tmp_zi[ityp][no+2],&tmp_zi[ityp][no+3]
                      ,&tmp_zi[ityp][no+4])!=5)
                error->abort("wrong line in eam file: %s",line);
            no+=5;
            
        }
        
        nlines=nrs[ityp]/5;
        no=1;
        for(int i=0;i<nlines;i++)
        {
            err=line_read(fp,line);
            if(err==-1)
                error->abort("eam potential file ended immaturely");
            
            if(sscanf(line,"%lf %lf %lf %lf %lf"
                      ,&tmp_rho[ityp][no],&tmp_rho[ityp][no+1]
                      ,&tmp_rho[ityp][no+2],&tmp_rho[ityp][no+3]
                      ,&tmp_rho[ityp][no+4])!=5)
                error->abort("wrong line in eam file: %s",line);
            no+=5;
            narg=mapp->parse_line(line,arg);
            
        }
        
        if(atoms->my_p_no==0)
            fclose(fp);
    }
    delete [] line;
    
    
    TYPE0 maxr=0.0;
    TYPE0 maxrho=0.0;
    TYPE0 maxdr=0.0;
    TYPE0 maxdrho=0.0;
    
    for(int i=0;i<no_types;i++)
    {
        maxr=MAX(maxr,static_cast<TYPE0>(nrs[i]-1)*drs[i]);
        maxrho=MAX(maxrho,static_cast<TYPE0>(nrhos[i]-1)*drhos[i]);
        maxdr=MAX(maxdr,drs[i]);
        maxdrho=MAX(maxdrho,drhos[i]);
    }
    
    nr=static_cast<int>(maxr/maxdr+0.5);
    nrho=static_cast<int>(maxr/maxdr+0.5);
    dr=maxdr;
    drho=maxdrho;
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    allocate();
    
    TYPE0 r,p,coef1,coef2,coef3,coef4;
    int k;
    TYPE0 sixth=1.0/6.0;
    
    for(int i=0;i<no_types;i++)
    {
        for(int j=0;j<nrho;j++)
        {
            r=static_cast<TYPE0>(j)*drho;
            p=r/drhos[i]+1.0;
            k=static_cast<int> (p);
            k=MIN(k,nrhos[i]-2);
            k=MAX(k,2);
            p-=k;
            p=MIN(p,2.0);
            coef1=-sixth*p*(p-1.0)*(p-2.0);
            coef2=0.5*(p*p-1.0)*(p-2.0);
            coef3=-0.5*p*(p+1.0)*(p-2.0);
            coef4=sixth*p*(p*p-1.0);
            F_arr[i][j][0]=coef1*tmp_F[i][k-1]+
            coef2*tmp_F[i][k]+coef3*tmp_F[i][k+1]+
            coef4*tmp_F[i][k+2];
        }
    }
    
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        no=ityp*(ityp+1);
        for(int i=0;i<nr;i++)
        {
            r=static_cast<TYPE0>(i)*dr;
            p=r/drs[ityp]+1.0;
            k=static_cast<int> (p);
            k=MIN(k,nrs[ityp]-2);
            k=MAX(k,2);
            p-=k;
            p=MIN(p,2.0);
            coef1=-sixth*p*(p-1.0)*(p-2.0);
            coef2=0.5*(p*p-1.0)*(p-2.0);
            coef3=-0.5*p*(p+1.0)*(p-2.0);
            coef4=sixth*p*(p*p-1.0);
            
            rho_arr[no][i][0]=coef1*tmp_rho[ityp][k-1]+
            coef2*tmp_rho[ityp][k]+coef3*tmp_rho[ityp][k+1]+
            coef4*tmp_rho[ityp][k+2];
        }
    }
    
    TYPE0 tmp0,tmp1;
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        for(int jtyp=ityp;jtyp<no_types;jtyp++)
        {
            no=COMP(ityp,jtyp);
            
            for(int i=0;i<nr;i++)
            {
                r=static_cast<TYPE0>(i)*dr;
                
                p=r/drs[ityp]+1.0;
                k=static_cast<int>(p);
                k=MIN(k,nrs[ityp]-2);
                k=MAX(k,2);
                p-=k;
                p=MIN(p,2.0);
                coef1=-sixth*p*(p-1.0)*(p-2.0);
                coef2=0.5*(p*p-1.0)*(p-2.0);
                coef3=-0.5*p*(p+1.0)*(p-2.0);
                coef4=sixth*p*(p*p-1.0);
                tmp0=coef1*tmp_zi[ityp][k-1]+
                coef2*tmp_zi[ityp][k]+coef3*tmp_zi[ityp][k+1]+
                coef4*tmp_zi[ityp][k+2];
                
                p=r/drs[jtyp] + 1.0;
                k=static_cast<int>(p);
                k=MIN(k,nrs[jtyp]-2);
                k=MAX(k,2);
                p-=k;
                p=MIN(p,2.0);
                coef1=-sixth*p*(p-1.0)*(p-2.0);
                coef2=0.5*(p*p-1.0)*(p-2.0);
                coef3=-0.5*p*(p+1.0)*(p-2.0);
                coef4=sixth*p*(p*p-1.0);
                tmp1=coef1*tmp_zi[jtyp][k-1]+
                coef2*tmp_zi[jtyp][k]+coef3*tmp_zi[jtyp][k+1]+
                coef4*tmp_zi[jtyp][k+2];
                
                phi_r_arr[no][i][0]=27.2*0.529*tmp0*tmp1;
            }
        }
    }
    
    
    
    for(int ityp=0;ityp<no_types;ityp++)
    {
        if(nrhos[ityp])
            delete [] tmp_F[ityp];
        
        if(nrs[ityp])
        {
            delete [] tmp_rho[ityp];
            delete [] tmp_zi[ityp];
        }
    }
    if(no_types)
    {
        delete [] tmp_F;
        delete [] tmp_rho;
        delete [] tmp_zi;
        delete [] drhos;
        delete [] drs;
        delete [] nrhos;
        delete [] nrs;
    }
    
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam_dmd_full::set_funcfl(int no_files
                                    ,char** file_names)
{
    if(no_files!=1)
        error->abort("one file is needed for eam");
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    fp=NULL;
    if(atoms->my_p_no==0)
    {
        fp=fopen(file_names[0],"r");
        if(fp==NULL)
            error->abort("file %s not found",file_names[0]);
    }
    
    for(int i=0;i<4;i++)
        if(line_read(fp,line)==-1)
            error->abort("eam potential file ended immaturely");
    
    char** arg;
    int narg;
    
    narg=mapp->parse_line(line,arg);
    if(narg<2 || narg-1 < no_types)
        error->abort("wrong format in eam potential file");
    
    int tot_no_types=narg-1;
    
    int** type_ref;
    CREATE1D(type_ref,no_types);
    for(int i=0;i<no_types;i++)
        CREATE1D(type_ref[i],2);
    
    for(int i=0;i<no_types;i++)
    {
        int ityp=0;
        while(strcmp(arg[ityp+1],atom_types->atom_names[i])!=0
              &&ityp<tot_no_types)
            ityp++;
        
        if(ityp==narg-1)
            error->abort("atom type not found in eam file");
        type_ref[i][0]=ityp;
        type_ref[i][1]=i;
    }
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    delete [] arg;
    
    for(int i=0;i<no_types;i++)
    {
        for(int j=i+1;j<no_types;j++)
        {
            if(type_ref[i][0]>type_ref[j][0])
            {
                int tmp0,tmp1;
                
                tmp0=type_ref[i][0];
                tmp1=type_ref[i][1];
                
                type_ref[i][0]=type_ref[j][0];
                type_ref[i][1]=type_ref[j][1];
                
                type_ref[j][0]=tmp0;
                type_ref[j][1]=tmp1;
            }
        }
    }
    
    if(line_read(fp,line)==-1)
        error->abort("eam potential file ended immaturely");
    narg=mapp->parse_line(line,arg);
    
    
    if(narg!=5)
        error->abort("wrong format in eam potential file");
    
    nrho=atoi(arg[0]);
    nr=atoi(arg[2]);
    drho=atof(arg[1]);
    dr=atof(arg[3]);
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    
    if(nrho%5!=0 || nr%5!=0)
        error->abort("nro and nrho remainder by 5 should be zero");
    
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    delete [] arg;
    
    
    if(line_read(fp,line)==-1)
        error->abort("eam potential file ended immaturely");
    
    
    
    allocate();
    
    
    int no_lines;
    int icur_pos,jcur_pos;
    int err;
    int no;
    int component;
    
    
    icur_pos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ityp==type_ref[icur_pos][0])
        {
            no=0;
            component=type_ref[icur_pos][1];
            no_lines=nrho/5;
            for(int i=0;i<no_lines;i++)
            {
                err=line_read(fp,line);
                if(err==-1)
                    error->abort("eam potential file ended immaturely");
                
                if(sscanf(line,"%lf %lf %lf %lf %lf"
                          ,&F_arr[component][no][0],&F_arr[component][no+1][0]
                          ,&F_arr[component][no+2][0],&F_arr[component][no+3][0]
                          ,&F_arr[component][no+4][0])!=5)
                    error->abort("wrong line in eam file: %s",line);
                no+=5;
                
            }
            
            no=0;
            component=type2rho[type_ref[icur_pos][1]][0];
            no_lines=nr/5;
            for(int i=0;i<no_lines;i++)
            {
                err=line_read(fp,line);
                if(err==-1)
                    error->abort("eam potential file ended immaturely");
                
                if(sscanf(line,"%lf %lf %lf %lf %lf"
                          ,&rho_arr[component][no][0],&rho_arr[component][no+1][0]
                          ,&rho_arr[component][no+2][0],&rho_arr[component][no+3][0]
                          ,&rho_arr[component][no+4][0])!=5)
                    error->abort("wrong line in eam file: %s",line);
                no+=5;
                
            }
            
            icur_pos++;
        }
        else
        {
            no_lines=nr/5+nrho/5;
            for(int i=0;i<no_lines;i++)
                if(line_read(fp,line)==-1)
                    error->abort("eam potential file ended immaturely");
        }
    }
    
    
    no_lines=nr/5;
    icur_pos=0;
    jcur_pos=0;
    
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ityp==type_ref[icur_pos][0])
        {
            jcur_pos=0;
            for(int jtyp=0;jtyp<ityp+1;jtyp++)
            {
                if(jtyp==type_ref[jcur_pos][0])
                {
                    no=0;
                    component=type2phi[type_ref[icur_pos][1]][type_ref[jcur_pos][1]];
                    no_lines=nr/5;
                    for(int i=0;i<no_lines;i++)
                    {
                        err=line_read(fp,line);
                        if(err==-1)
                            error->abort("eam potential file ended immaturely");
                        
                        if(sscanf(line,"%lf %lf %lf %lf %lf"
                                  ,&phi_r_arr[component][no][0],&phi_r_arr[component][no+1][0]
                                  ,&phi_r_arr[component][no+2][0],&phi_r_arr[component][no+3][0]
                                  ,&phi_r_arr[component][no+4][0])!=5)
                            error->abort("wrong line in eam file: %s",line);
                        no+=5;
                        
                    }
                    
                    jcur_pos++;
                }
                else
                {
                    no_lines=nr/5;
                    for(int i=0;i<no_lines;i++)
                        if(line_read(fp,line)==-1)
                            error->abort("eam potential file ended immaturely");
                }
            }
            
            
            icur_pos++;
        }
        else
        {
            no_lines=nr*(ityp+1)/5+nrho/5;
            for(int i=0;i<no_lines;i++)
                if(line_read(fp,line)==-1)
                    error->abort("eam potential file ended immaturely");
        }
        
    }
    
    
    delete [] line;
    
    for(int i=0;i<no_types;i++)
        delete [] type_ref[i];
    if(no_types)
        delete [] type_ref;
    if(atoms->my_p_no==0)
        fclose(fp);
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam_dmd_full::set_fs(int no_files
                                ,char** file_names)
{
    if(no_files!=1)
        error->abort("one file is needed for eam");
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    
    fp=NULL;
    if(atoms->my_p_no==0)
    {
        fp=fopen(file_names[0],"r");
        if(fp==NULL)
            error->abort("file %s not found",file_names[0]);
    }
    
    for(int i=0;i<4;i++)
        if(line_read(fp,line)==-1)
            error->abort("eam potential file ended immaturely");
    
    char** arg;
    int narg;
    
    narg=mapp->parse_line(line,arg);
    if(narg<2 || narg-1 < no_types)
        error->abort("wrong format in eam potential file");
    
    int tot_no_types=narg-1;
    
    int** type_ref;
    CREATE1D(type_ref,no_types);
    for(int i=0;i<no_types;i++)
        CREATE1D(type_ref[i],2);
    
    for(int i=0;i<no_types;i++)
    {
        int ityp=0;
        while(strcmp(arg[ityp+1],atom_types->atom_names[i])!=0
              &&ityp<tot_no_types)
            ityp++;
        
        if(ityp==narg-1)
            error->abort("atom type not found in eam file");
        type_ref[i][0]=ityp;
        type_ref[i][1]=i;
    }
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    delete [] arg;
    
    for(int i=0;i<no_types;i++)
    {
        for(int j=i+1;j<no_types;j++)
        {
            if(type_ref[i][0]>type_ref[j][0])
            {
                int tmp0,tmp1;
                
                tmp0=type_ref[i][0];
                tmp1=type_ref[i][1];
                
                type_ref[i][0]=type_ref[j][0];
                type_ref[i][1]=type_ref[j][1];
                
                type_ref[j][0]=tmp0;
                type_ref[j][1]=tmp1;
            }
        }
    }
    
    if(line_read(fp,line)==-1)
        error->abort("eam potential file ended immaturely");
    narg=mapp->parse_line(line,arg);
    
    
    if(narg!=5)
        error->abort("wrong format in eam potential file");
    
    nrho=atoi(arg[0]);
    nr=atoi(arg[2]);
    drho=atof(arg[1]);
    dr=atof(arg[3]);
    dr_inv=1.0/dr;
    drho_inv=1.0/drho;
    
    if(nrho%5!=0 || nr%5!=0)
        error->abort("nro and nrho remainder by 5 should be zero");
    
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    delete [] arg;
    
    
    if(line_read(fp,line)==-1)
        error->abort("eam potential file ended immaturely");
    
    
    
    allocate();
    
    
    int no_lines;
    int icur_pos,jcur_pos;
    int err;
    int no;
    int component;
    
    icur_pos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        
        
        if(ityp==type_ref[icur_pos][0])
        {
            no=0;
            component=type_ref[icur_pos][1];
            no_lines=nrho/5;
            for(int i=0;i<no_lines;i++)
            {
                err=line_read(fp,line);
                if(err==-1)
                    error->abort("eam potential file ended immaturely");
                
                if(sscanf(line,"%lf %lf %lf %lf %lf"
                          ,&F_arr[component][no][0],&F_arr[component][no+1][0]
                          ,&F_arr[component][no+2][0],&F_arr[component][no+3][0]
                          ,&F_arr[component][no+4][0])!=5)
                    error->abort("wrong line in eam file: %s",line);
                no+=5;
                
            }
            
            jcur_pos=0;
            for(int jtyp=0;jtyp<tot_no_types;jtyp++)
            {
                if(jtyp==type_ref[jcur_pos][0])
                {
                    no=0;
                    component=type2rho[type_ref[icur_pos][1]][type_ref[jcur_pos][1]];
                    no_lines=nr/5;
                    for(int i=0;i<no_lines;i++)
                    {
                        err=line_read(fp,line);
                        if(err==-1)
                            error->abort("eam potential file ended immaturely");
                        
                        if(sscanf(line,"%lf %lf %lf %lf %lf"
                                  ,&rho_arr[component][no][0],&rho_arr[component][no+1][0]
                                  ,&rho_arr[component][no+2][0],&rho_arr[component][no+3][0]
                                  ,&rho_arr[component][no+4][0])!=5)
                            error->abort("wrong line in eam file: %s",line);
                        no+=5;
                        
                    }
                    
                    jcur_pos++;
                }
                else
                {
                    no_lines=nr/5;
                    for(int i=0;i<no_lines;i++)
                        if(line_read(fp,line)==-1)
                            error->abort("eam potential file ended immaturely");
                    
                }
            }
            
            icur_pos++;
        }
        else
        {
            no_lines=nr*tot_no_types/5+nrho/5;
            for(int i=0;i<no_lines;i++)
                if(line_read(fp,line)==-1)
                    error->abort("eam potential file ended immaturely");
        }
    }
    
    
    
    no_lines=nr/5;
    icur_pos=0;
    jcur_pos=0;
    
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ityp==type_ref[icur_pos][0])
        {
            jcur_pos=0;
            for(int jtyp=0;jtyp<ityp+1;jtyp++)
            {
                if(jtyp==type_ref[jcur_pos][0])
                {
                    no=0;
                    component=type2phi[type_ref[icur_pos][1]][type_ref[jcur_pos][1]];
                    no_lines=nr/5;
                    for(int i=0;i<no_lines;i++)
                    {
                        err=line_read(fp,line);
                        if(err==-1)
                            error->abort("eam potential file ended immaturely");
                        
                        if(sscanf(line,"%lf %lf %lf %lf %lf"
                                  ,&phi_r_arr[component][no][0],&phi_r_arr[component][no+1][0]
                                  ,&phi_r_arr[component][no+2][0],&phi_r_arr[component][no+3][0]
                                  ,&phi_r_arr[component][no+4][0])!=5)
                            error->abort("wrong line in eam file: %s",line);
                        no+=5;
                        
                    }
                    
                    jcur_pos++;
                }
                else
                {
                    no_lines=nr/5;
                    for(int i=0;i<no_lines;i++)
                        if(line_read(fp,line)==-1)
                            error->abort("eam potential file ended immaturely");
                    
                }
            }
            
            icur_pos++;
        }
        else
        {
            no_lines=nr*(ityp+1)/5+nrho/5;
            for(int i=0;i<no_lines;i++)
                if(line_read(fp,line)==-1)
                    error->abort("eam potential file ended immaturely");
        }
        
    }
    
    
    delete [] line;
    
    for(int i=0;i<no_types;i++)
        delete [] type_ref[i];
    if(no_types)
        delete [] type_ref;
    if(atoms->my_p_no==0)
        fclose(fp);
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
int ForceField_eam_dmd_full::line_read(FILE* file
                                  ,char*& line)
{
    int lenght;
    int eof=0;
    if(atoms->my_p_no==0)
        if(feof(file))
            eof=-1;
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if(eof==-1)
        return -1;
    
    if(atoms->my_p_no==0)
    {
        fgets(line,MAXCHAR,file);
        lenght=static_cast<int>(strlen(line))+1;
    }
    MPI_Bcast(&lenght,1,MPI_INT,0,world);
    MPI_Bcast(line,lenght,MPI_CHAR,0,world);
    
    return lenght;
}
/*--------------------------------------------
 clean up the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_full::clean_up()
{
    if(allocated==0)
        return;
    
    
    for(int i=0;i<no_types;i++)
    {
        delete [] type2phi[i];
        delete [] type2rho[i];
    }
    if(no_types)
    {
        delete [] type2phi;
        delete [] type2rho;
    }
    
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int i=0;i<no_types*(no_types+1)/2;i++)
        {
            for(int j=0;j<nr;j++)
                delete [] phi_r_arr[i][j];
            delete [] phi_r_arr[i];
        }
        delete [] phi_r_arr;
        
        for(int i=0;i<no_types;i++)
        {
            for(int j=0;j<nr;j++)
                delete [] rho_arr[i][j];
            delete [] rho_arr[i];
        }
        delete [] rho_arr;
        
        
        for(int i=0;i<no_types;i++)
        {
            for(int j=0;j<nrho;j++)
                delete [] F_arr[i][j];
            delete [] F_arr[i];
        }
        delete [] F_arr;
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int i=0;i<no_types*(no_types+1)/2;i++)
        {
            for(int j=0;j<nr;j++)
                delete [] phi_r_arr[i][j];
            delete [] phi_r_arr[i];
        }
        delete [] phi_r_arr;
        
        for(int i=0;i<no_types*no_types;i++)
        {
            for(int j=0;j<nr;j++)
                delete [] rho_arr[i][j];
            delete [] rho_arr[i];
        }
        delete [] rho_arr;
        
        for(int i=0;i<no_types;i++)
        {
            for(int j=0;j<nrho;j++)
                delete [] F_arr[i][j];
            delete [] F_arr[i];
        }
        delete [] F_arr;
    }
    
    allocated=0;
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_full::allocate()
{
    
    CREATE2D(type2phi,no_types,no_types);
    CREATE2D(type2rho,no_types,no_types);
    
    for(int i=0;i<no_types;i++)
        for(int j=0;j<no_types;j++)
        {
            type2phi[i][j]=COMP(i,j);
            type2rho[i][j]=j*no_types+i;
        }
    
    
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        for(int i=0;i<no_types*(no_types+1)/2;i++)
            CREATE1D(phi_r_arr[i],nr);
        for(int i=0;i<no_types*(no_types+1)/2;i++)
            for(int j=0;j<nr;j++)
                CREATE1D(phi_r_arr[i][j],7);
        
        CREATE1D(rho_arr,no_types*no_types);
        for(int i=0;i<no_types;i++)
            CREATE1D(rho_arr[i],nr);
        for(int i=0;i<no_types;i++)
            for(int j=0;j<nr;j++)
                CREATE1D(rho_arr[i][j],7);
        
        for(int i=1;i<no_types;i++)
            for(int j=0;j<no_types;j++)
                rho_arr[type2rho[i][j]]
                =rho_arr[type2rho[i][0]];
        
        CREATE1D(F_arr,no_types);
        for(int i=0;i<no_types;i++)
            CREATE1D(F_arr[i],nrho+1);
        for(int i=0;i<no_types;i++)
            for(int j=0;j<nrho+1;j++)
                CREATE1D(F_arr[i][j],9);
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        for(int i=0;i<no_types*(no_types+1)/2;i++)
            CREATE1D(phi_r_arr[i],nr);
        for(int i=0;i<no_types*(no_types+1)/2;i++)
            for(int j=0;j<nr;j++)
                CREATE1D(phi_r_arr[i][j],7);
        
        CREATE1D(rho_arr,no_types*no_types);
        for(int i=0;i<no_types*no_types;i++)
            CREATE1D(rho_arr[i],nr);
        for(int i=0;i<no_types*no_types;i++)
            for(int k=0;k<nr;k++)
                CREATE1D(rho_arr[i][k],7);
        
        CREATE1D(F_arr,no_types);
        for(int i=0;i<no_types;i++)
            CREATE1D(F_arr[i],nrho+1);
        for(int i=0;i<no_types;i++)
            for(int j=0;j<nrho+1;j++)
                CREATE1D(F_arr[i][j],9);
    }
    
    allocated=1;
}
/*--------------------------------------------
 set the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_full::set_arrays()
{
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate_m(nrho,drho,F_arr[i]);
        
        for(int i=0;i<no_types;i++)
            interpolate(nr,dr,rho_arr[i]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<i+1;j++)
                interpolate(nr,dr,phi_r_arr[type2phi[i][j]]);
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate_m(nrho,drho,F_arr[i]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<no_types;j++)
                interpolate(nr,dr,rho_arr[type2rho[i][j]]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<i+1;j++)
                interpolate(nr,dr,phi_r_arr[type2phi[i][j]]);
    }
}
/*--------------------------------------------
 allocate the arrays
 --------------------------------------------*/
void ForceField_eam_dmd_full::interpolate(int n,TYPE0 delta
,TYPE0** spline)
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for (int i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
                      8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for (int i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;
    for(int i=0;i<n;i++)
    {
        spline[i][4]=spline[i][1]/delta;
        spline[i][5]=2.0*spline[i][2]/delta;
        spline[i][6]=3.0*spline[i][3]/delta;
    }
    
}
/*--------------------------------------------
  allocate the arrays
--------------------------------------------*/
void ForceField_eam_dmd_full::interpolate_m(int n,TYPE0 delta
,TYPE0** spline)
{
    spline[0][1]=spline[1][0]-spline[0][0];
    spline[1][1]=0.5*(spline[2][0]-spline[0][0]);
    spline[n-2][1]=0.5*(spline[n-1][0]-spline[n-3][0]);
    spline[n-1][1]=spline[n-1][0]-spline[n-2][0];
    
    for (int i=2;i<n-2;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
                      8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for (int i=0;i<n-1;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    
    spline[n-1][2]=0.0;
    spline[n-1][3]=0.0;
    
    
    
    for (int i=0;i<n-1;i++)
    {
        spline[i][4]=(spline[i+1][2]-spline[i][2])/6.0-0.5*spline[i][3];
    }
    spline[n-1][4]=0.0;
    
    
    
    for(int i=0;i<n;i++)
    {
        spline[i][5]=spline[i][1]/delta;
        spline[i][6]=2.0*spline[i][2]/delta;
        spline[i][7]=3.0*spline[i][3]/delta;
        spline[i][8]=4.0*spline[i][4]/delta;
    }
    
}
/*--------------------------------------------
 Gaussian-Hermite quadrature weights and
 abscissas for 1 to 14 points
 --------------------------------------------*/
void ForceField_eam_dmd_full::set_weight_abs(int n)
{
    
    if(n<1 || n>14)
        error->abort("n should be between 1 and 14");
    
    if(no_i)
    {
        delete [] xi;
        delete [] wi_0;
        delete [] wi_1;
        delete [] wi_2;
    }
    
    no_i=n;
    
    CREATE1D(xi,n);
    CREATE1D(wi_0,n);
    CREATE1D(wi_1,n);
    CREATE1D(wi_2,n);
    
    if(no_i==1)
    {
        wi_0[0]=1.7724538509055160273;
        
        xi[0]=0.0;
    }
    else if(no_i==2)
    {
        wi_0[0]=0.8862269254527580136;
        wi_0[1]=0.8862269254527580136;
        
        xi[0]=-0.70710678118654752440;
        xi[1]=0.70710678118654752440;
    }
    else if(no_i==3)
    {
        wi_0[0]=0.2954089751509193379;
        wi_0[1]=1.1816359006036773515;
        wi_0[2]=0.2954089751509193379;
        
        xi[0]=-1.2247448713915890491;
        xi[1]=0.0;
        xi[2]=1.2247448713915890491;
    }
    else if(no_i==4)
    {
        wi_0[0]=0.0813128354472451771;
        wi_0[1]=0.8049140900055128365;
        wi_0[2]=0.8049140900055128365;
        wi_0[3]=0.0813128354472451771;
        
        xi[0]=-1.6506801238857845559;
        xi[1]=-0.52464762327529031788;
        xi[2]=0.52464762327529031788;
        xi[3]=1.6506801238857845559;
    }
    else if(no_i==5)
    {
        wi_0[0]=0.01995324205904591321;
        wi_0[1]=0.3936193231522411598;
        wi_0[2]=0.94530872048294188123;
        wi_0[3]=0.3936193231522411598;
        wi_0[4]=0.01995324205904591321;
        
        xi[0]=-2.0201828704560856329;
        xi[1]=-0.95857246461381850711;
        xi[2]=0.0;
        xi[3]=0.95857246461381850711;
        xi[4]=2.0201828704560856329;
    }
    else if(no_i==6)
    {
        wi_0[0]=0.00453000990550884564;
        wi_0[1]=0.1570673203228566439;
        wi_0[2]=0.7246295952243925241;
        wi_0[3]=0.7246295952243925241;
        wi_0[4]=0.1570673203228566439;
        wi_0[5]=0.00453000990550884564;
        
        xi[0]=-2.3506049736744922228;
        xi[1]=-1.3358490740136969497;
        xi[2]=-0.43607741192761650868;
        xi[3]=0.43607741192761650868;
        xi[4]=1.3358490740136969497;
        xi[5]=2.3506049736744922228;
    }
    else if(no_i==7)
    {
        wi_0[0]=0.00097178124509951915;
        wi_0[1]=0.0545155828191270306;
        wi_0[2]=0.4256072526101278005;
        wi_0[3]=0.81026461755680732676;
        wi_0[4]=0.4256072526101278005;
        wi_0[5]=0.0545155828191270306;
        wi_0[6]=0.00097178124509951915;
        
        xi[0]=-2.6519613568352334924;
        xi[1]=-1.6735516287674714450;
        xi[2]=-0.81628788285896466304;
        xi[3]=0.0;
        xi[4]=0.81628788285896466304;
        xi[5]=1.6735516287674714450;
        xi[6]=2.6519613568352334924;
    }
    else if(no_i==8)
    {
        wi_0[0]=0.000199604072211367619;
        wi_0[1]=0.0170779830074134755;
        wi_0[2]=0.207802325814891880;
        wi_0[3]=0.6611470125582412910;
        wi_0[4]=0.6611470125582412910;
        wi_0[5]=0.207802325814891880;
        wi_0[6]=0.0170779830074134755;
        wi_0[7]=0.000199604072211367619;
        
        xi[0]=-2.9306374202572440192;
        xi[1]=-1.9816567566958429259;
        xi[2]=-1.1571937124467801947;
        xi[3]=-0.38118699020732211685;
        xi[4]=0.38118699020732211685;
        xi[5]=1.1571937124467801947;
        xi[6]=1.9816567566958429259;
        xi[7]=2.9306374202572440192;
    }
    else if(no_i==9)
    {
        wi_0[0]=0.000039606977263264382;
        wi_0[1]=0.0049436242755369472;
        wi_0[2]=0.088474527394376573;
        wi_0[3]=0.432651559002555750;
        wi_0[4]=0.72023521560605095712;
        wi_0[5]=0.432651559002555750;
        wi_0[6]=0.088474527394376573;
        wi_0[7]=0.0049436242755369472;
        wi_0[8]=0.000039606977263264382;
        
        xi[0]=-3.1909932017815276072;
        xi[1]=-2.2665805845318431118;
        xi[2]=-1.4685532892166679317;
        xi[3]=-0.72355101875283757332;
        xi[4]=0.0;
        xi[5]=0.72355101875283757332;
        xi[6]=1.4685532892166679317;
        xi[7]=2.2665805845318431118;
        xi[8]=3.1909932017815276072;
    }
    else if(no_i==10)
    {
        wi_0[0]=7.6404328552326206e-6;
        wi_0[1]=0.0013436457467812327;
        wi_0[2]=0.033874394455481063;
        wi_0[3]=0.240138611082314686;
        wi_0[4]=0.6108626337353257988;
        wi_0[5]=0.6108626337353257988;
        wi_0[6]=0.240138611082314686;
        wi_0[7]=0.033874394455481063;
        wi_0[8]=0.0013436457467812327;
        wi_0[9]=7.6404328552326206e-6;
        
        xi[0]=-3.4361591188377376033;
        xi[1]=-2.5327316742327897964;
        xi[2]=-1.7566836492998817735;
        xi[3]=-1.0366108297895136542;
        xi[4]=-0.34290132722370460879;
        xi[5]=0.34290132722370460879;
        xi[6]=1.0366108297895136542;
        xi[7]=1.7566836492998817735;
        xi[8]=2.5327316742327897964;
        xi[9]=3.4361591188377376033;
    }
    else if(no_i==11)
    {
        wi_0[0]=1.4395603937142582e-6;
        wi_0[1]=0.00034681946632334551;
        wi_0[2]=0.011911395444911532;
        wi_0[3]=0.117227875167708503;
        wi_0[4]=0.429359752356125028;
        wi_0[5]=0.65475928691459177920;
        wi_0[6]=0.429359752356125028;
        wi_0[7]=0.117227875167708503;
        wi_0[8]=0.011911395444911532;
        wi_0[9]=0.00034681946632334551;
        wi_0[10]=1.4395603937142582e-6;
        
        xi[0]=-3.6684708465595825185;
        xi[1]=-2.7832900997816517708;
        xi[2]=-2.0259480158257553352;
        xi[3]=-1.3265570844949328559;
        xi[4]=-0.65680956688209976502;
        xi[5]=0.0;
        xi[6]=0.65680956688209976502;
        xi[7]=1.3265570844949328559;
        xi[8]=2.0259480158257553352;
        xi[9]=2.7832900997816517708;
        xi[10]=3.6684708465595825185;
    }
    else if(no_i==12)
    {
        wi_0[0]=2.6585516843563016e-7;
        wi_0[1]=0.00008573687043587859;
        wi_0[2]=0.003905390584629062;
        wi_0[3]=0.051607985615883930;
        wi_0[4]=0.260492310264161129;
        wi_0[5]=0.5701352362624795783;
        wi_0[6]=0.5701352362624795783;
        wi_0[7]=0.260492310264161129;
        wi_0[8]=0.051607985615883930;
        wi_0[9]=0.003905390584629062;
        wi_0[10]=0.00008573687043587859;
        wi_0[11]=2.6585516843563016e-7;
        
        xi[0]=-3.8897248978697819193;
        xi[1]=-3.0206370251208897717;
        xi[2]=-2.2795070805010599002;
        xi[3]=-1.5976826351526047967;
        xi[4]=-0.94778839124016374370;
        xi[5]=-0.31424037625435911128;
        xi[6]=0.31424037625435911128;
        xi[7]=0.94778839124016374370;
        xi[8]=1.5976826351526047967;
        xi[9]=2.2795070805010599002;
        xi[10]=3.0206370251208897717;
        xi[11]=3.8897248978697819193;
    }
    else if(no_i==13)
    {
        wi_0[0]=4.825731850073131e-8;
        wi_0[1]=0.00002043036040270707;
        wi_0[2]=0.0012074599927193859;
        wi_0[3]=0.020862775296169939;
        wi_0[4]=0.140323320687023438;
        wi_0[5]=0.421616296898543222;
        wi_0[6]=0.60439318792116164234;
        wi_0[7]=0.421616296898543222;
        wi_0[8]=0.140323320687023438;
        wi_0[9]=0.020862775296169939;
        wi_0[10]=0.0012074599927193859;
        wi_0[11]=0.00002043036040270707;
        wi_0[12]=4.825731850073131e-8;
        
        xi[0]=-4.1013375961786396412;
        xi[1]=-3.2466089783724099881;
        xi[2]=-2.5197356856782378834;
        xi[3]=-1.8531076516015121420;
        xi[4]=-1.2200550365907484262;
        xi[5]=-0.60576387917106011308;
        xi[6]=0.0;
        xi[7]=0.60576387917106011308;
        xi[8]=1.2200550365907484262;
        xi[9]=1.8531076516015121420;
        xi[10]=2.5197356856782378834;
        xi[11]=3.2466089783724099881;
        xi[12]=4.1013375961786396412;
    }
    else if(no_i==14)
    {
        wi_0[0]=8.628591168125158e-9;
        wi_0[1]=4.716484355018917e-6;
        wi_0[2]=0.0003550926135519236;
        wi_0[3]=0.007850054726457944;
        wi_0[4]=0.06850553422346521;
        wi_0[5]=0.273105609064246603;
        wi_0[6]=0.5364059097120901498;
        wi_0[7]=0.5364059097120901498;
        wi_0[8]=0.273105609064246603;
        wi_0[9]=0.06850553422346521;
        wi_0[10]=0.007850054726457944;
        wi_0[11]=0.0003550926135519236;
        wi_0[12]=4.716484355018917e-6;
        wi_0[13]=8.628591168125158e-9;
        
        xi[0]=-4.3044485704736318126;
        xi[1]=-3.4626569336022705502;
        xi[2]=-2.7484707249854025686;
        xi[3]=-2.0951832585077168157;
        xi[4]=-1.4766827311411408706;
        xi[5]=-0.87871378732939941611;
        xi[6]=-0.29174551067256207845;
        xi[7]=0.29174551067256207845;
        xi[8]=0.87871378732939941611;
        xi[9]=1.4766827311411408706;
        xi[10]=2.0951832585077168157;
        xi[11]=2.7484707249854025686;
        xi[12]=3.4626569336022705502;
        xi[13]=4.3044485704736318126;
    }
    
    for(int i=0;i<no_i;i++)
    {
        wi_1[i]=wi_0[i]*xi[i];
        wi_2[i]=wi_0[i]*xi[i]*xi[i];
    }
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_full::create_2nd_neigh_lst()
{
    
    TYPE0 dx0,dx1,dx2,rsq;
    int iatm,jatm,icomp,jcomp;
    int natms=atoms->natms;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    int* crd;
    atoms->vectors[crd_n].ret(crd);
    
    int* tmp_neigh_list;
    int tmp_neigh_list_size=1024;
    CREATE1D(tmp_neigh_list,1024);
    
    
    if(neigh_lst_sz_sz)
    {
        for(int i=0;i<neigh_lst_sz_sz;i++)
            if(neigh_lst_sz[i])
                delete [] neigh_lst[i];
        
        delete [] neigh_lst_sz;
        delete [] neigh_lst;
    }
    neigh_lst_sz_sz=natms;
    
    CREATE1D(neigh_lst,natms);
    CREATE1D(neigh_lst_sz,natms);
    
    for(int i=0;i<natms;i++)
        crd[i]=neigh_lst_sz[i]=0;
    
    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=(3+no_types)*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=(3+no_types)*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            if(rsq<rsq_crd)
            {
                
                if(neigh_lst_sz[iatm]+1>tmp_neigh_list_size)
                {
                    GROW(tmp_neigh_list,tmp_neigh_list_size,tmp_neigh_list_size+1);
                    tmp_neigh_list_size++;
                }
                tmp_neigh_list[neigh_lst_sz[iatm]]=jatm;
                neigh_lst_sz[iatm]++;
                crd[iatm]++;
                
                if(jatm<natms)
                    crd[jatm]++;
            }
            
            CREATE1D(neigh_lst[iatm],neigh_lst_sz[iatm]);
            memcpy(neigh_lst[iatm],tmp_neigh_list,neigh_lst_sz[iatm]*sizeof(int));
            
        }
    }
    
    
    if(tmp_neigh_list_size)
        delete [] tmp_neigh_list;
    
    atoms->update(crd_n);

}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
TYPE0 ForceField_eam_dmd_full::
Mat(int crd_i,int crd_j,int itype)
{
    return 5.0*beta;
}
/*--------------------------------------------
 return M_{ij}^{\alpha}
 --------------------------------------------*/
TYPE0 ForceField_eam_dmd_full::calc_ent(TYPE0 x)
{
    TYPE0 tmp=(x-0.5)*(x-0.5);
    TYPE0 ans=-0.679898676683507
    +2.14126441923694*tmp
    -1.512638507830134*tmp*tmp
    +15.303838631273946*tmp*tmp*tmp;
    
    return ans;
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
void ForceField_eam_dmd_full::calc_y()
{
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* E;
    atoms->vectors[E_n].ret(E);
    TYPE0* dE;
    atoms->vectors[dE_n].ret(dE);
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* n;
    atoms->vectors[n_n].ret(n);
    int* crd;
    atoms->vectors[crd_n].ret(crd);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    TYPE0* coef;
    TYPE0 p,p2,p3,p4,tmp0,tmp1,s_ij;
    
    int m;
    int iatm,jatm;
    int icomp,jcomp,phi_start,rho_start;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    int natms=atoms->natms;
    for(iatm=0;iatm<natms*no_types;iatm++)
        c_d[iatm]=E[iatm]=dE[iatm]=0.0;
    
    phi_start=0;
    rho_start=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    E[icomp+itype]+=c[jcomp+jtype]*rho[rho_start+type2rho[jtype][itype]];
                    if(jatm<natms)
                        E[jcomp+jtype]+=c[icomp+itype]*rho[rho_start+type2rho[itype][jtype]];
                }
            }
            
            rho_start+=no_types*no_types;
        }
        
        for(int itype=0;itype<no_types; itype++)
        {
            p=E[iatm*no_types+itype]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            p4=p3*p;
            coef=F_arr[itype][m];
            
            tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            tmp1=coef[8]*p3+coef[7]*p2+coef[6]*p+coef[5];
            
            if(E[iatm*no_types+itype]>rho_max)
                tmp0+=tmp1*(E[icomp+itype]-rho_max);
            
            n[icomp+itype]=E[icomp+itype]=tmp0;
            dE[icomp+itype]=tmp1;
        }
    }
    

    atoms->update(dE_n);
    
    phi_start=0;
    rho_start=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    n[icomp+itype]+=c[jcomp+jtype]*(rho[rho_start+type2rho[itype][jtype]]*dE[jatm+jtype]
                                                    +phi[phi_start+type2phi[itype][jtype]]);
                    if(jatm<natms)
                        n[jcomp+jtype]+=c[icomp+itype]*(rho[rho_start+type2rho[jtype][itype]]*dE[iatm+itype]
                                                        +phi[phi_start+type2phi[itype][jtype]]);
                }
            }
            rho_start+=no_types*no_types;
            phi_start+=no_types*(no_types+1)/2;
        }
        
        for(int itype=0;itype<no_types;itype++)
        {
            n[icomp+itype]+=c_0[itype]+1.5*kbT*log(x[(3+no_types)*iatm+3+itype]);
            n[icomp+itype]*=beta;
        }
    }
    
    atoms->update(n_n);

    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neigh_lst_sz[iatm];j++)
        {
            jatm=neigh_lst[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                tmp0=exp(n[icomp+itype]+Mat(crd[iatm],crd[jatm],itype));
                tmp1=exp(n[jcomp+itype]+Mat(crd[iatm],crd[jatm],itype));
                
                s_ij=((1.0-c[icomp+itype])*c[jcomp+itype]*tmp1
                  -(1.0-c[jcomp+itype])*c[icomp+itype]*tmp0);

                c_d[icomp+itype]+=s_ij;
                if(jatm<natms)
                    c_d[jcomp+itype]-=s_ij;
            }
        }
    }
}
/*--------------------------------------------
 claculate F and dF and dFF
 --------------------------------------------*/
TYPE0 ForceField_eam_dmd_full::calc_g(int chk,TYPE0 alpha,
TYPE0* a,TYPE0* g)
{
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* E;
    atoms->vectors[E_n].ret(E);
    TYPE0* dE;
    atoms->vectors[dE_n].ret(dE);
    TYPE0* ddE;
    atoms->vectors[ddE_n].ret(ddE);
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* n;
    atoms->vectors[n_n].ret(n);
    TYPE0* s;
    atoms->vectors[s_n].ret(s);
    int* crd;
    atoms->vectors[crd_n].ret(crd);
    TYPE0* t;
    atoms->vectors[t_n].ret(t);
    TYPE0* v;
    atoms->vectors[v_n].ret(v);
    TYPE0 inner,ans;
    
    
    TYPE0* coef;
    
    TYPE0 p,p2,p3,p4,tmp0,tmp1,s_ij;
    int m;
    
    
    int iatm,jatm;
    
    int icomp,jcomp,phi_start,rho_start;
    
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    int natms=atoms->natms;
    for(int i=0;i<natms*no_types;i++)
        g[i]=E[i]=0.0;
    
    phi_start=0;
    rho_start=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    E[icomp+itype]+=c[jcomp+jtype]*rho[rho_start+type2rho[jtype][itype]];
                    if(jatm<natms)
                        E[jcomp+jtype]+=c[icomp+itype]*rho[rho_start+type2rho[itype][jtype]];
                }
            }
            
            rho_start+=no_types*no_types;
            
        }
        
        for(int itype=0;itype<no_types; itype++)
        {
            p=E[iatm*no_types+itype]*drho_inv;
            m=static_cast<int>(p);
            m=MIN(m,nr-2);
            p-=m;
            p=MIN(p,1.0);
            p2=p*p;
            p3=p2*p;
            p4=p3*p;
            coef=F_arr[itype][m];

            
            tmp0=coef[4]*p4+coef[3]*p3+coef[2]*p2+coef[1]*p+coef[0];
            tmp1=coef[8]*p3+coef[7]*p2+coef[6]*p+coef[5];
            
            if(E[iatm*no_types+itype]>rho_max)
                tmp0+=tmp1*(E[icomp+itype]-rho_max);
            
            n[icomp+itype]=E[icomp+itype]=tmp0;
            dE[icomp+itype]=tmp1;
            ddE[icomp+itype]=(3.0*coef[8]*p2+coef[7]*p+coef[6])*drho_inv;
        }
    }
    
    atoms->update(dE_n);
    
    phi_start=0;
    rho_start=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    n[icomp+itype]+=c[jcomp+jtype]*(rho[rho_start+type2rho[itype][jtype]]*dE[jatm+jtype]
                                                    +phi[phi_start+type2phi[itype][jtype]]);
                    if(jatm<natms)
                        n[jcomp+jtype]+=c[icomp+itype]*(rho[rho_start+type2rho[jtype][itype]]*dE[iatm+itype]
                                                        +phi[phi_start+type2phi[itype][jtype]]);
                }
            }
            rho_start+=no_types*no_types;
            phi_start+=no_types*(no_types+1)/2;
        }
        
        for(int itype=0;itype<no_types;itype++)
        {
            n[icomp+itype]+=c_0[itype];
            n[icomp+itype]*=beta;
            n[icomp+itype]+=1.5*log(x[(3+no_types)*iatm+3+itype]);
        }
        
    }
    
    atoms->update(n_n);
    
    for(int i=0;i<natms*no_types;i++) s[i]=2.0*(c[i]+a[i]);
    
    inner=0.0;
    for(iatm=0;iatm<natms;iatm++)
    {
        
        icomp=no_types*iatm;
        for(int j=0;j<neigh_lst_sz[iatm];j++)
        {
            jatm=neigh_lst[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                tmp0=exp(n[icomp+itype]+Mat(crd[iatm],crd[jatm],itype));
                tmp1=exp(n[jcomp+itype]+Mat(crd[iatm],crd[jatm],itype));
                
                s_ij=-2.0*alpha*((1.0-c[icomp+itype])*c[jcomp+itype]*tmp1-(1.0-c[jcomp+itype])*c[icomp+itype]*tmp0);
                s[icomp+itype]+=s_ij;
                if(jatm<natms)
                    s[jcomp+itype]-=s_ij;
            }
        }
        for(int itype=0;itype<no_types;itype++)
        {
            inner+=s[icomp+itype]*s[icomp+itype];
        }
    }
    ans=0.0;
    MPI_Allreduce(&inner,&ans,1,MPI_TYPE0,MPI_SUM,world);
    ans*=0.25;
    
    
    if(chk)
        return ans;
    
    
    atoms->update(s_n);
    
    for(int i=0;i<natms*no_types;i++)
    {
        v[i]=t[i]=0.0;
        g[i]=-s[i];
    }
    
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neigh_lst_sz[iatm];j++)
        {
            jatm=neigh_lst[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                tmp0=(s[icomp+itype]-s[jcomp+itype])*exp(Mat(crd[iatm],crd[jatm],itype)+n[icomp+itype]);
                tmp1=(s[icomp+itype]-s[jcomp+itype])*exp(Mat(crd[iatm],crd[jatm],itype)+n[jcomp+itype]);
                
                
                t[icomp+itype]+=c[icomp+itype]*tmp0*(1.0-c[jcomp+itype]);
                g[icomp+itype]-=alpha*(tmp0*(1.0-c[jcomp+itype])+tmp1*c[jcomp+itype]);
                if(jatm<natms)
                {
                    t[jcomp+itype]-=c[jcomp+itype]*tmp1*(1.0-c[icomp+itype]);
                    g[jcomp+itype]+=alpha*(tmp1*(1.0-c[icomp+itype])+c[icomp+itype]*tmp0);
                }
            }
        }
    }
    
    atoms->update(t_n);
    
    phi_start=0;
    rho_start=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    
                    tmp0=phi[phi_start+type2phi[jtype][itype]]
                    +rho[rho_start+type2rho[jtype][itype]]*dE[icomp+itype]
                    +rho[rho_start+type2rho[itype][jtype]]*dE[jcomp+jtype];
                    tmp0*=beta*alpha;
                    
                    v[icomp+itype]+=rho[rho_start+type2rho[jtype][itype]]*t[jcomp+jtype];
                    g[icomp+itype]+=tmp0*t[jcomp+jtype];
                    if(jatm<natms)
                    {
                        v[jcomp+jtype]+=rho[rho_start+type2rho[itype][jtype]]*t[icomp+itype];
                        g[jcomp+jtype]+=tmp0*t[icomp+itype];
                    }
                }
            }
            rho_start+=no_types*no_types;
            phi_start+=no_types*(no_types+1)/2;
        }
    }
    
    atoms->update(v_n);
    atoms->update(ddE_n);
    
    rho_start=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        icomp=no_types*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jcomp=no_types*jatm;
            for(int itype=0;itype<no_types;itype++)
            {
                for(int jtype=0;jtype<no_types;jtype++)
                {
                    g[icomp+itype]+=alpha*beta*rho[rho_start+type2rho[itype][jtype]]
                    *ddE[jcomp+jtype]*c[jcomp+jtype]*v[jcomp+jtype];
                    if(jatm<natms)
                        g[jcomp+jtype]+=alpha*beta*rho[rho_start+type2rho[jtype][itype]]
                        *ddE[icomp+itype]*c[icomp+itype]*v[icomp+itype];
                }
            }
            rho_start+=no_types*no_types;
        }
    }

    return ans;
    
}



