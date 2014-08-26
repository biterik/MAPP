
#include "neighbor.h"
#include "ff_dmd.h"
#include "atom_types.h"
using namespace MAPP_NS;
#define PI_IN_SQ 0.564189583547756286948079451561
enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_DMD::
ForceField_DMD(MAPP* mapp) : ForceField(mapp)
{
    no_i=0;
    allocated=0;
    eam_mode=NOT_SET;
    max_pairs=0;
    
    if(atoms->vectors[0].dim!=4)
        error->abort("the dimension of x"
                     " vector should be 4");
    
    
    int no_types=atom_types->no_types;
    CREATE1D(nrgy_strss,7);
    CREATE1D(cut_sk_sq,no_types*(no_types+1)/2);
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_DMD::~ForceField_DMD()
{
    int no_types=atom_types->no_types;
    
    if(no_types) delete [] cut_sk_sq;
    
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
        delete [] drhoi_dr;
        delete [] drhoj_dr;
        delete [] drhoi_dbeta;
        delete [] drhoj_dbeta;
    }
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_DMD::
force_calc(int st_clc,TYPE0* en_st)
{
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    TYPE0* rho;
    atoms->vectors[rho_n].ret(rho);
    TYPE0* dF;
    atoms->vectors[dF_n].ret(dF);
    int* type;
    atoms->vectors[type_n].ret(type);
    
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp,istart;
    TYPE0 dx0,dx1,dx2,rsq,z2p,z2;
    TYPE0 r,p,r_inv=0.0;
    int m;
    TYPE0* coef;
    TYPE0 tmp0,tmp1,tmp2,tmp3;
    TYPE0 upper,lower;
    TYPE0 rho_i,drho_i_dbeta,drho_i_dr
    ,rho_j,drho_j_dbeta,drho_j_dr
    ,phi,dphi_dbeta,dphi_dr;
    TYPE0 coef0,coef1,coef2,coef3,coef4;
    TYPE0 beta_sq,beta,beta_inv;
    TYPE0 g_0_phi,g_1_phi,g_2_phi;
    TYPE0 g_0_rho_i,g_1_rho_i,g_2_rho_i;
    TYPE0 g_0_rho_j,g_1_rho_j,g_2_rho_j;
    TYPE0 fpair,bpair;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    
    if(max_pairs<neighbor->no_pairs)
    {
        if(max_pairs)
        {
            delete [] drhoi_dr;
            delete [] drhoj_dr;
            delete [] drhoi_dbeta;
            delete [] drhoj_dbeta;
        }
        max_pairs=neighbor->no_pairs;
        CREATE1D(drhoi_dr,max_pairs);
        CREATE1D(drhoj_dr,max_pairs);
        CREATE1D(drhoi_dbeta,max_pairs);
        CREATE1D(drhoj_dbeta,max_pairs);
    }
    
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    
    int natms=atoms->natms;
    for(iatm=0;iatm<natms;iatm++) rho[iatm]=0.0;
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=4*iatm;
        nrgy_strss[0]-=3.0*kbT*log(x[icomp+3]);
        f[icomp+3]+=3.0*kbT/x[icomp+3];
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            
            jcomp=4*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            beta_sq=(x[icomp+3]*x[icomp+3]+x[jcomp+3]*x[jcomp+3]);

            
            drhoi_dr[istart]=0.0;
            drhoi_dbeta[istart]=0.0;
            drhoj_dr[istart]=0.0;
            drhoj_dbeta[istart]=0.0;
            if(rsq<cut_sq_mod)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                beta=sqrt(beta_sq);
                beta_inv=1.0/beta;
                phi=rho_i=rho_j=0.0;
                dphi_dbeta=drho_i_dbeta=drho_j_dbeta=0.0;
                dphi_dr=drho_i_dr=drho_j_dr=0.0;
                
                if(beta_min<beta && beta<beta_max)
                {
                    
                    upper=(r+rc)*beta_inv;
                    lower=(r-rc)*beta_inv;
                    g_0_phi=g_1_phi=g_2_phi=
                    g_0_rho_i=g_1_rho_i=g_2_rho_i=
                    g_0_rho_j=g_1_rho_j=g_2_rho_j=0.0;
                    
                    if(lower<xi[no_i-1])
                    {
                        for(int i=0;i<no_i;i++)
                        {
                            if(lower<xi[i] && xi[i]<upper)
                            {
                                tmp0=r-xi[i]*beta;
                                
                                p=fabs(tmp0)*dr_inv+1.0;
                                m=static_cast<int>(p);
                                m=MIN(m,nr-1);
                                p-=m;
                                p=MIN(p,1.0);
                                coef=phi_r_arr[type2rho[itype][jtype]][m];
                                tmp1=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                                if(tmp0<0.0)
                                    tmp1*=-1.0;
                                
                                g_0_phi+=wi_0[i]*tmp1;
                                g_1_phi+=wi_1[i]*tmp1;
                                g_2_phi+=wi_2[i]*tmp1;

                                
                                coef=rho_arr[type2rho[jtype][itype]][m];
                                tmp2=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                                coef=rho_arr[type2rho[itype][jtype]][m];
                                tmp3=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                                tmp2*=tmp0;
                                tmp3*=tmp0;
                                
                                g_0_rho_i+=wi_0[i]*tmp2;
                                g_1_rho_i+=wi_1[i]*tmp2;
                                g_2_rho_i+=wi_2[i]*tmp2;
                                
                                g_0_rho_j+=wi_0[i]*tmp3;
                                g_1_rho_j+=wi_1[i]*tmp3;
                                g_2_rho_j+=wi_2[i]*tmp3;
                            }
                        }
                        
                        coef0=PI_IN_SQ*r_inv;
                        coef1=-coef0*r_inv;
                        coef3=-coef0*beta_inv;
                        coef2=2.0*coef3;
                        coef4=-coef2;

                        
                        phi=coef0*g_0_phi;
                        dphi_dr=coef1*g_0_phi+coef2*g_1_phi;
                        dphi_dbeta=coef3*g_0_phi+coef4*g_2_phi;
                        
                        rho_i=coef0*g_0_rho_i;
                        drho_i_dr=coef1*g_0_rho_i+coef2*g_1_rho_i;
                        drho_i_dbeta=coef3*g_0_rho_i+coef4*g_2_rho_i;
                        
                        rho_j=coef0*g_0_rho_j;
                        drho_j_dr=coef1*g_0_rho_j+coef2*g_1_rho_j;
                        drho_j_dbeta=coef3*g_0_rho_j+coef4*g_2_rho_j;
                    }
                    
                }
                else if (beta<=beta_min)
                {
                    if(rsq < cut_sq)
                    {
                        r=sqrt(rsq);
                        p=r*dr_inv+1.0;
                        m=static_cast<int>(p);
                        m=MIN(m,nr-1);
                        p-=m;
                        p=MIN(p,1.0);
                        
                        coef=rho_arr[type2rho[jtype][itype]][m];
                        rho_i=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        drho_i_dr=(coef[6]*p+coef[5])*p+coef[4];
                        coef=rho_arr[type2rho[itype][jtype]][m];
                        rho_j=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        drho_j_dr=(coef[6]*p+coef[5])*p+coef[4];
                        
                        coef=phi_r_arr[type2phi[itype][jtype]][m];
                        z2p=(coef[6]*p + coef[5])*p+coef[4];
                        z2=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        
                        phi=z2*r_inv;
                        dphi_dr=z2p*r_inv-phi*r_inv;
                        
                    }
                }
                
                
                fpair=-dphi_dr*r_inv;
                bpair=-2.0*dphi_dbeta;
                rho[iatm]+=rho_i;
                f[icomp]+=dx0*fpair;
                f[icomp+1]+=dx1*fpair;
                f[icomp+2]+=dx2*fpair;
                f[icomp+3]+=bpair*x[icomp+3];
                
                if(jatm<natms)
                {
                    rho[jatm]+=rho_j;
                    
                    f[jcomp]-=dx0*fpair;
                    f[jcomp+1]-=dx1*fpair;
                    f[jcomp+2]-=dx2*fpair;
                    f[jcomp+3]+=bpair*x[jcomp+3];
                    
                    nrgy_strss[0]+=phi;
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
                    nrgy_strss[0]+=0.5*phi;
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
                
                drhoi_dr[istart]=-drho_i_dr*r_inv;
                drhoi_dbeta[istart]=-2.0*drho_i_dbeta;
                
                drhoj_dr[istart]=-drho_j_dr*r_inv;
                drhoj_dbeta[istart]=-2.0*drho_j_dbeta;
            }
            
            
            istart++;
        }
    }
    
    atoms->update(rho_n);
    
    
    int tot_natms=atoms->natms+atoms->natms_ph;
    
    for(iatm=0;iatm<tot_natms;iatm++)
    {
        
        p=rho[iatm]*drho_inv+1.0;
        m=static_cast<int> (p);
        m=MIN(m,nr-1);
        p-=m;
        p=MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        dF[iatm]=(coef[6]*p+coef[5])*p+coef[4];
        
        if(iatm<natms)
        {
            tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
            if(rho[iatm]>rho_max)
                tmp0+=dF[iatm]*(rho[iatm]-rho_max);
            nrgy_strss[0]+=tmp0;
        }
    }
    
    
    istart=0;
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=4*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            jcomp=4*jatm;
            
            fpair=dF[iatm]*drhoi_dr[istart]+dF[jatm]*drhoj_dr[istart];
            bpair=dF[iatm]*drhoi_dbeta[istart]+dF[jatm]*drhoj_dbeta[istart];
            
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            
            f[icomp]+=dx0*fpair;
            f[icomp+1]+=dx1*fpair;
            f[icomp+2]+=dx2*fpair;
            f[icomp+3]+=bpair*x[icomp+3];
            
            if(jatm<natms)
            {
                f[jcomp]-=dx0*fpair;
                f[jcomp+1]-=dx1*fpair;
                f[jcomp+2]-=dx2*fpair;
                f[jcomp+3]+=bpair*x[jcomp+3];
                
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
    en_st[0]+=c_0;
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
TYPE0 ForceField_DMD::energy_calc()
{
    
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* rho;
    atoms->vectors[rho_n].ret(rho);
    TYPE0* dF;
    atoms->vectors[dF_n].ret(dF);
    int* type;
    atoms->vectors[type_n].ret(type);
    
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,z2;
    TYPE0 r,p,r_inv=0.0;
    int m;
    TYPE0* coef;
    TYPE0 tmp0,tmp1,tmp2,tmp3;
    TYPE0 upper,lower;
    TYPE0 rho_i,rho_j,phi;
    TYPE0 coef0;
    TYPE0 beta_sq,beta,beta_inv;
    TYPE0 g_0_phi;
    TYPE0 g_0_rho_i;
    TYPE0 g_0_rho_j;
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    

    int natms=atoms->natms;
    for(iatm=0;iatm<natms;iatm++) rho[iatm]=0.0;
    TYPE0 en=0.0;
    TYPE0 en_tot=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=4*iatm;
        
        en-=3.0*kbT*log(x[icomp+3]);
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            
            jcomp=4*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            beta_sq=(x[icomp+3]*x[icomp+3]+x[jcomp+3]*x[jcomp+3]);
            
            if(rsq<cut_sq_mod)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                beta=sqrt(beta_sq);
                phi=rho_i=rho_j=0.0;
                
                if(beta_min<beta && beta<beta_max)
                {
                    beta_inv=1.0/beta;
                    
                    upper=(r+rc)*beta_inv;
                    lower=(r-rc)*beta_inv;
                    g_0_phi=0.0;
                    g_0_rho_i=0.0;
                    g_0_rho_j=0.0;
                    
                    if(lower<xi[no_i-1])
                    {
                        for(int i=0;i<no_i;i++)
                        {
                            if(lower<xi[i] && xi[i]<upper)
                            {
                                tmp0=r-xi[i]*beta;
                                
                                p=fabs(tmp0)*dr_inv+1.0;
                                m=static_cast<int>(p);
                                m=MIN(m,nr-1);
                                p-=m;
                                p=MIN(p,1.0);
                                coef=phi_r_arr[type2rho[itype][jtype]][m];
                                tmp1=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                                if(tmp0<0.0)
                                    tmp1*=-1.0;
                                
                                g_0_phi+=wi_0[i]*tmp1;
                                
                                coef=rho_arr[type2rho[jtype][itype]][m];
                                tmp2=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                                coef=rho_arr[type2rho[itype][jtype]][m];
                                tmp3=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                                tmp2*=tmp0;
                                tmp3*=tmp0;
                                
                                g_0_rho_i+=wi_0[i]*tmp2;
                                
                                g_0_rho_j+=wi_0[i]*tmp3;
                            }
                        }
                        
                        
                        coef0=PI_IN_SQ*r_inv;

                        
                        phi=coef0*g_0_phi;
                        
                        rho_i=coef0*g_0_rho_i;
                        rho_j=coef0*g_0_rho_j;

                    }
                    
                }
                else if (beta<=beta_min)
                {
                    if(rsq < cut_sq)
                    {
                        r=sqrt(rsq);
                        p=r*dr_inv+1.0;
                        m=static_cast<int>(p);
                        m=MIN(m,nr-1);
                        p-=m;
                        p=MIN(p,1.0);
                        
                        coef=rho_arr[type2rho[jtype][itype]][m];
                        rho_i=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        coef=rho_arr[type2rho[itype][jtype]][m];
                        rho_j=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        
                        coef=phi_r_arr[type2phi[itype][jtype]][m];
                        z2=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                        
                        phi=z2*r_inv;
                        
                    }
                }
                
                rho[iatm]+=rho_i;
                
                if(jatm<natms)
                {
                    en+=phi;
                    rho[jatm]+=rho_j;
                }
                else
                    en+=0.5*phi;
                
            }
        }
        
        p=rho[iatm]*drho_inv+1.0;
        m=static_cast<int> (p);
        m=MIN(m,nr-1);
        p-=m;
        p=MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[iatm]>rho_max)
            tmp0+=((coef[6]*p+coef[5])*p+coef[4])*(rho[iatm]-rho_max);
        en+=tmp0;
        
    }
    

    MPI_Allreduce(&en,&en_tot,1,MPI_TYPE0,MPI_SUM,world);
    en_tot+=c_0;
    return en_tot;
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_DMD::init()
{
    TYPE0 skin=atoms->skin;
    TYPE0 ph_cut=0.0;
    int no_types=atom_types->no_types;
    for (int i=0;i<no_types*(no_types+1)/2;i++)
        cut_sk_sq[i]=cut_sq_mod+(skin)*(skin)
        +2*sqrt(cut_sq_mod)*(skin);
    
    ph_cut=sqrt(cut_sq_mod);
    
    atoms->set_ph(ph_cut);
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    type_n=atoms->find("type");
    
    neighbor->pair_wise=1;
    
    rho_n=atoms->add<TYPE0>(1,1,"rho");
    dF_n=atoms->add<TYPE0>(1,1,"dF");
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_DMD::fin()
{
    atoms->del(dF_n);
    atoms->del(rho_n);
    
    if(max_pairs)
    {
        delete [] drhoi_dr;
        delete [] drhoj_dr;
        delete [] drhoi_dbeta;
        delete [] drhoj_dbeta;
        max_pairs=0;
    }
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_DMD::coef(int narg,char** arg)
{
    TYPE0 alpha_min,alpha_max;
    TYPE0 kb,T,hbar;
    
    set_weight_abs(atoi(arg[1]));
    alpha_min=atof(arg[2]);
    alpha_max=atof(arg[3]);
    T=atof(arg[4]);
    kb=atof(arg[5]);
    hbar=atof(arg[6]);
    
    clean_up();
    if(strcmp(arg[7],"FS")==0)
    {
        eam_mode=FINNIS_FL;
        set_funcfl(narg-8,&arg[8]);
    }
    else if(strcmp(arg[7],"SetFL")==0)
    {
        eam_mode=SET_FL;
        set_setfl(narg-8,&arg[8]);
    }
    else if(strcmp(arg[7],"FuncFL")==0)
    {
        eam_mode=FUNC_FL;
        set_funcfl(narg-8,&arg[8]);
    }
    else
        error->abort("wrong coeff command "
                     "for eam Force Field");
    set_arrays();
    
    int no_types=atom_types->no_types;
    TYPE0 mass;
    TYPE0 deb_l;
    TYPE0* vec_0;
    CREATE1D(vec_0,no_types);
    for(int i=0;i<no_types;i++)
    {
        mass=atom_types->mass[i];
        mass*=1.0364269184093291236e-28;
        deb_l=hbar*2.506628274631/sqrt(mass*kb*T);
        vec_0[i]=1.5*kb*T*log(0.1170996630*deb_l*deb_l);
    }
    
    int natms=atoms->natms;
    int* type;
    atoms->vectors[atoms->find("type")].ret(type);
    int no_itype;
    int tot_no_itype;
    c_0=0.0;
    for(int itype=0;itype<no_types;itype++)
    {
        no_itype=0;
        for(int i=0;i<natms;i++)
            if(type[i]==itype)
                no_itype++;
        
        tot_no_itype=0;
        MPI_Allreduce(&no_itype,&tot_no_itype,1,MPI_INT,MPI_SUM,world);
        c_0+=static_cast<TYPE0>(tot_no_itype)*vec_0[itype];
    }
    kbT=kb*T;
    delete [] vec_0;
    
    if(alpha_min==0.0)
        error->abort("minimum alpha cannot be zero");
    
    beta_min=1.0/sqrt(alpha_max);
    beta_max=1.0/sqrt(alpha_min);
    
    rc=(static_cast<TYPE0>(nr)-1.0)*dr;
    rho_max=(static_cast<TYPE0>(nrho)-1.0)*drho;
    cut_sq=rc*rc;
    mod_rc=rc+xi[no_i-1]*beta_max;
    cut_sq_mod=mod_rc*mod_rc;
    
    TYPE0* x;
    atoms->vectors[0].ret(x);
    int tot_natms=atoms->natms+atoms->natms_ph;
    
    TYPE0 alpha_ave=0.5*(alpha_min+alpha_max);
    TYPE0 beta_ave=0.5*(beta_min+beta_max);
    alpha_ave=20.0;
    beta_ave=1.0/sqrt(alpha_ave);
    for(int i=0;i<tot_natms;i++)
        x[i*4+3]=beta_ave;
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_DMD::set_setfl(int no_files
,char** file_names)
{
    int no_types=atom_types->no_types;
    if(no_files!=no_types)
        error->abort("no of types and number files should be equal in eam");
    
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
        rho_arr[i][0][0]=0.0;
        for(int j=1;j<nrho+1;j++)
        {
            r=(j-1)*drho;
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
        rho_arr[no][0][0]=0.0;
        for(int i=1;i<nr+1;i++)
        {
            r=(i-1)*dr;
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
            
            phi_r_arr[no][0][0]=0.0;
            for(int i=1;i<nr+1;i++)
            {
                r=(i-1)*dr;
                
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
void ForceField_DMD::set_funcfl(int no_files
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
    int no_types=atom_types->no_types;
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
            no=1;
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
            
            no=1;
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
                    no=1;
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
void ForceField_DMD::set_fs(int no_files
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
    int no_types=atom_types->no_types;
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
            no=1;
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
                    no=1;
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
                    no=1;
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
int ForceField_DMD::line_read(FILE* file
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
void ForceField_DMD::clean_up()
{
    if(allocated==0)
        return;
    
    int no_types=atom_types->no_types;
    
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
void ForceField_DMD::allocate()
{
    int no_types=atom_types->no_types;
    
    CREATE2D(type2phi,no_types,no_types);
    CREATE2D(type2rho,no_types,no_types);
    
    for(int i=0;i<no_types;i++)
        for(int j=0;j<no_types;j++)
        {
            type2phi[i][j]=COMP(i,j);
            type2rho[i][j]=i*no_types+j;
        }
    
    
    
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        for(int i=0;i<no_types*(no_types+1)/2;i++)
            CREATE1D(phi_r_arr[i],nr+1);
        for(int i=0;i<no_types*(no_types+1)/2;i++)
            for(int j=0;j<nr+1;j++)
                CREATE1D(phi_r_arr[i][j],7);
        
        CREATE1D(rho_arr,no_types*no_types);
        for(int i=0;i<no_types;i++)
            CREATE1D(rho_arr[i],nr+1);
        for(int i=0;i<no_types;i++)
            for(int j=0;j<nr+1;j++)
                CREATE1D(rho_arr[i][j],7);
        
        for(int i=1;i<no_types;i++)
            for(int j=0;j<no_types;j++)
                rho_arr[type2rho[i][j]]
                =rho_arr[type2rho[0][j]];
        
        CREATE1D(F_arr,no_types);
        for(int i=0;i<no_types;i++)
            CREATE1D(F_arr[i],nrho+1);
        for(int i=0;i<no_types;i++)
            for(int j=0;j<nrho+1;j++)
                CREATE1D(F_arr[i][j],7);
        
    }
    else if(eam_mode==FINNIS_FL)
    {
        
        CREATE1D(phi_r_arr,no_types*(no_types+1)/2);
        for(int i=0;i<no_types*(no_types+1)/2;i++)
            CREATE1D(phi_r_arr[i],nr+1);
        for(int i=0;i<no_types*(no_types+1)/2;i++)
            for(int j=0;j<nr+1;j++)
                CREATE1D(phi_r_arr[i][j],7);
        
        CREATE1D(rho_arr,no_types*no_types);
        for(int i=0;i<no_types*no_types;i++)
            CREATE1D(rho_arr[i],nr+1);
        for(int i=0;i<no_types*no_types;i++)
            for(int k=0;k<nr+1;k++)
                CREATE1D(rho_arr[i][k],7);
        
        CREATE1D(F_arr,no_types);
        for(int i=0;i<no_types;i++)
            CREATE1D(F_arr[i],nrho+1);
        for(int i=0;i<no_types;i++)
            for(int j=0;j<nrho+1;j++)
                CREATE1D(F_arr[i][j],7);
    }
    
    allocated=1;
}
/*--------------------------------------------
 set the arrays
 --------------------------------------------*/
void ForceField_DMD::set_arrays()
{
    int no_types=atom_types->no_types;
    if(eam_mode==FUNC_FL || eam_mode==SET_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate(nrho,drho,F_arr[i]);
        
        for(int i=0;i<no_types;i++)
            interpolate(nr,dr,rho_arr[i]);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<i+1;j++)
                interpolate(nr,dr,phi_r_arr[type2phi[i][j]]);
    }
    else if(eam_mode==FINNIS_FL)
    {
        for(int i=0;i<no_types;i++)
            interpolate(nrho,drho,F_arr[i]);
        
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
void ForceField_DMD::interpolate(int n,TYPE0 delta
,TYPE0** spline)
{
    spline[1][1]=spline[2][0]-spline[1][0];
    spline[2][1]=0.5*(spline[3][0]-spline[1][0]);
    spline[n-1][1]=0.5*(spline[n][0]-spline[n-2][0]);
    spline[n][1]=spline[n][0]-spline[n-1][0];
    
    for (int i=3;i<n-1;i++)
        spline[i][1]=((spline[i-2][0]-spline[i+2][0])+
                      8.0*(spline[i+1][0]-spline[i-1][0]))/12.0;
    
    for (int i=1;i<n;i++)
    {
        spline[i][2]=3.0*(spline[i+1][0]-spline[i][0])-
        2.0*spline[i][1]-spline[i+1][1];
        spline[i][3]=spline[i][1]+spline[i+1][1]-
        2.0*(spline[i+1][0]-spline[i][0]);
    }
    
    spline[n][2]=0.0;
    spline[n][3]=0.0;
    for(int i=1;i<n+1;i++)
    {
        spline[i][4]=spline[i][1]/delta;
        spline[i][5]=2.0*spline[i][2]/delta;
        spline[i][6]=3.0*spline[i][3]/delta;
    }
    
}
/*--------------------------------------------
 Gaussian-Hermite quadrature weights and
 abscissas for 1 to 14 points
 --------------------------------------------*/
void ForceField_DMD::set_weight_abs(int n)
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

