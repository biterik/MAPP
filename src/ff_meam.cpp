//
//  ff_meam.cpp
//  MAPP
//
//  Created by Sina on 7/23/14.
//  Copyright (c) 2014 Li Group/Sina. All rights reserved.
//
#include "neighbor.h"
#include "ff_meam.h"
#include "atom_types.h"
using namespace MAPP_NS;

enum{FCC,BCC,HCP,DIM,DIAMOND,B1,C11,L12,B2};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_meam::
ForceField_meam(MAPP* mapp):ForceField(mapp)
{
    int no_types=atom_types->no_types;
    
    CREATE1D(c_min,no_types);
    CREATE1D(c_max,no_types);
    for(int i=0;i<no_types;i++)
    {
        CREATE1D(c_min[i],no_types);
        CREATE1D(c_max[i],no_types);
    }
    
    for(int j=0;j<no_types;j++)
    {
        for(int i=0;i<no_types;i++)
        {
            CREATE1D(c_min[i][j],no_types);
            CREATE1D(c_max[i][j],no_types);
        }
    }
    
    rho_dim=38;
    
    
    CREATE1D(v2d,6);
    CREATE1D(v3d,10);
    v2d[0]=1;
    v2d[1]=2;
    v2d[2]=2;
    v2d[3]=1;
    v2d[4]=2;
    v2d[5]=1;
    
    v3d[0]=1;
    v3d[1]=3;
    v3d[2]=3;
    v3d[3]=3;
    v3d[4]=6;
    v3d[5]=3;
    v3d[6]=1;
    v3d[7]=3;
    v3d[8]=3;
    v3d[9]=1;
    
    CREATE2D(vind2d,3,3);
    CREATE2D(vind3d,3,3);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            CREATE1D(vind3d[i][j],3);
    
    int comp0=0;
    int comp1=0;
    for(int i=0;i<3;i++)
    {
        for(int j=i;j<3;j++)
        {
            vind2d[i][j]=comp0;
            vind2d[j][i]=comp0;
            comp0++;
            for(int k=j;k<3;k++)
            {
                vind3d[i][j][k]=comp1;
                vind3d[i][k][j]=comp1;
                vind3d[j][i][k]=comp1;
                vind3d[j][k][i]=comp1;
                vind3d[k][i][j]=comp1;
                vind3d[k][j][i]=comp1;
                comp1++;
            }
        }
    }
    
    
    
    
    /*
     rho dimension is 38!!!!
**** rho functions (4)
     00: rho0
     01: rho1
     02: rho2
     03: rho3
     
**** rho_1 terms (3)
     04: Arho1_1
     05: Arho1_2
     06: Arho1_3
     
**** rho_2 terms (7)
******** 1st terms (6)
     07: Arho2_1
     08: Arho2_2
     09: Arho2_3
     10: Arho2_4
     11: Arho2_5
     12: Arho2_6
******** 2nd term (1)
     13: Arho2b
     
**** rho_3 terms
******** 1st terms (10)
     14: Arho3_1
     15: Arho3_2
     16: Arho3_3
     17: Arho3_4
     18: Arho3_5
     19: Arho3_6
     20: Arho3_7
     21: Arho3_8
     22: Arho3_9
     23: Arho3_10
******** 2nd terms (3)
     24: Arho3b_1
     25: Arho3b_2
     26: Arho3b_3
     
**** t_ave terms (3)
     27: t_ave_1
     28: t_ave_2
     29: t_ave_3

**** tsq_ave terms (3)
     30: tsq_ave_1
     31: tsq_ave_2
     32: tsq_ave_3
     
**** Gamma function (1)
     33: gamma

**** dGamma functions (3)
     34: dgamma1
     35: dgamma2
     36: dgamma3
     
**** fhop term  (1)
     37: fhop
     
     */
    
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_meam::~ForceField_meam()
{
    int no_types=atom_types->no_types;
    
    for(int j=0;j<no_types;j++)
    {
        for(int i=0;i<no_types;i++)
        {
            delete [] c_min[i][j];
            delete [] c_max[i][j];
        }
    }
    
    for(int i=0;i<no_types;i++)
    {
        delete [] c_min[i];
        delete [] c_max[i];
    }
    if(no_types)
    {
        delete [] c_min;
        delete [] c_max;
    }
    
    
    delete [] v2d;
    delete [] v3d;
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            delete [] vind3d[i][j];
    
    for(int i=0;i<3;i++)
        delete [] vind3d[i];
    delete [] vind3d;
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_meam::screen()
{
    int* my_list;
    int my_list_size;
    
    int jatm,katm;
    
    TYPE0 xij,xjk,xki,yij,yjk,yki,zij,zjk,zki;
    TYPE0 rij2,rjk2,rki2;
    TYPE0 rij;
    int icomp,jcomp,kcomp,itype,jtype,ktype,coef1,coef2;
    int icomp_rho,jcomp_rho;
    
    
    TYPE0 fcij,dfcij,sij,dsij;
    TYPE0 rnorm,rbound,x_ki,x_jk,a;
    TYPE0 cikj,sikj=0.0,cmax,cmin,dfikj,dcikj;
    TYPE0 ai,aj,ro0i,ro0j,rhoa0i,rhoa0j,rhoa1i,rhoa1j
    ,rhoa2i,rhoa2j,rhoa3i,rhoa3j;
    TYPE0 Z;
    TYPE0 G,dG,Gbar,dGbar,gam,rho_bkgd,rhob,denom,B;
    int iArho1_comp,jArho1_comp,iArho2_comp
    ,jArho2_comp,iArho3_comp,jArho3_comp
    ,iArhob3_comp,jArhob3_comp;
    TYPE0 A1i,A2i,A3i,A1j,A2j,A3j,B1i,B1j;
    
    TYPE0* delij;
    CREATE1D(delij,3);
    TYPE0* s;
    CREATE1D(s,3);
    
    
    int istart=0;
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        my_list_size=neighbor->neighbor_list_size[iatm];
        my_list=neighbor->neighbor_list[iatm];
        icomp=iatm*3;
        icomp_rho=iatm*rho_dim;
        itype=type[iatm];
        for(int jn=0;jn<my_list_size;jn++)
        {
            jatm=my_list[jn];
            if(jatm>iatm)
            {
                jtype=type[jatm];
                jcomp=jatm*3;
                jcomp_rho=jatm*rho_dim;
                
                xij=x[icomp]-x[jcomp];
                yij=x[icomp+1]-x[jcomp+1];
                zij=x[icomp+2]-x[jcomp+2];
                delij[0]=-xij;
                delij[1]=-yij;
                delij[2]=-zij;
                rij2=xij*xij+yij*yij+zij*zij;
                rij=sqrt(rij2);
                
                if(rij<rc_meam)
                {
                    rnorm=(rc_meam-rij)*delr_meam_inv;
                    sij=1.0;
                    dsij=0.0;
                    rbound = ebound_meam[itype][jtype]*rij2;
                    
                    int kn=0;
                    dfcut(rnorm,fcij,dfcij);
                    dfcij*=delr_meam_inv;
                    
                    while(kn<my_list_size && sij!=0.0)
                    {
                        katm=my_list[kn];
                        
                        if(katm!=jatm)
                        {
                            ktype=type[katm];
                            kcomp=3*katm;
                            
                            xjk=x[jcomp]-x[kcomp];
                            yjk=x[jcomp+1]-x[kcomp+1];
                            zjk=x[jcomp+2]-x[kcomp+2];
                            rjk2=xjk*xjk+yjk*yjk+zjk*zjk;
                            
                            if(rjk2<rbound)
                            {
                                xki=x[kcomp]-x[icomp];
                                yki=x[kcomp+1]-x[icomp+1];
                                zki=x[kcomp+2]-x[icomp+2];
                                rki2=xki*xki+yki*yki+zki*zki;
                                
                                if(rki2<rbound)
                                {
                                    x_ki=rki2/rij2;
                                    x_jk=rjk2/rij2;
                                    a=1.0-(x_ki-x_jk)*(x_ki-x_jk);
                                    if(a>0)
                                    {
                                        cikj=(2.0*(x_ki+x_jk)+a-2.0)/a;
                                        cmax=c_max[itype][jtype][ktype];
                                        cmin=c_min[itype][jtype][ktype];
                                        if(cikj<cmax&&cikj>cmin)
                                        {
                                            cikj=(cikj-cmin)/(cmax-cmin);
                                            dfcut(cikj,sikj,dfikj);
                                            sij*=sikj;
                                            coef1=dfikj/((cmax-cmin)*sikj);
                                            dCfunc(rij2,rki2,rjk2,dcikj);
                                            dsij+=coef1*dcikj;
                                        }
                                        else if(cikj<cmin)
                                        {
                                            sij=0.0;
                                            dsij=0.0;
                                        }
                                    }
                                }
                            }
                        }
                        kn++;
                    }
                    coef1=sij*fcij;
                    coef2=sij*dfcij/rij;
                    dsij*=coef1;
                    dsij-=coef2;
                    if(sij==1.0 || sij==0.0)
                        dsij=0.0;
                }
                else
                {
                    fcij=0.0;
                    dfcij=0.0;
                    sij=0.0;
                    dsij=0.0;
                }
                
                fcpair[istart]=fcij;
                scrfcn[istart]=sij;
                dscrfcn[istart]=dsij;
                
                if(rij<rc_meam && sij*fcij!=0.0)
                {
                    ai=rij/re_meam[itype][itype]-1.0;
                    aj=rij/re_meam[jtype][jtype]-1.0;
                    ro0i=rho0_meam[itype];
                    ro0j=rho0_meam[jtype];
                    
                    rhoa0i=sij*ro0i*exp(-beta0_meam[itype]*ai);
                    rhoa1i=sij*ro0i*exp(-beta1_meam[itype]*ai);
                    rhoa2i=sij*ro0i*exp(-beta2_meam[itype]*ai);
                    rhoa3i=sij*ro0i*exp(-beta3_meam[itype]*ai);
                    
                    rhoa0j=sij*ro0j*exp(-beta0_meam[jtype]*aj);
                    rhoa1j=sij*ro0j*exp(-beta1_meam[jtype]*aj);
                    rhoa2j=sij*ro0j*exp(-beta2_meam[jtype]*aj);
                    rhoa3j=sij*ro0j*exp(-beta3_meam[jtype]*aj);
                    
                    
                    if(ialloy==1)
                    {
                        rhoa1i*=t1_meam[itype];
                        rhoa2i*=t2_meam[itype];
                        rhoa3i*=t3_meam[itype];
                        
                        rhoa1j*=t1_meam[jtype];
                        rhoa2j*=t2_meam[jtype];
                        rhoa3j*=t3_meam[jtype];
                    }
                    
                    rho[icomp_rho]+=rhoa0j;
                    rho[jcomp_rho]+=rhoa0i;
                    
                    
                    if(ialloy!=2)
                    {
                        
                        rho[icomp_rho+27]+=t1_meam[jtype]*rhoa0j;
                        rho[icomp_rho+28]+=t2_meam[jtype]*rhoa0j;
                        rho[icomp_rho+29]+=t3_meam[jtype]*rhoa0j;
                        rho[jcomp_rho+27]+=t1_meam[itype]*rhoa0i;
                        rho[jcomp_rho+28]+=t2_meam[itype]*rhoa0i;
                        rho[jcomp_rho+29]+=t3_meam[itype]*rhoa0i;
                    }
                    if(ialloy==1)
                    {
                        
                        rho[icomp_rho+30]+=t1_meam[jtype]*t1_meam[jtype]*rhoa0j;
                        rho[icomp_rho+31]+=t2_meam[jtype]*t2_meam[jtype]*rhoa0j;
                        rho[icomp_rho+32]+=t3_meam[jtype]*t3_meam[jtype]*rhoa0j;
                        rho[jcomp_rho+30]+=t1_meam[itype]*t1_meam[itype]*rhoa0i;
                        rho[jcomp_rho+31]+=t2_meam[itype]*t2_meam[itype]*rhoa0i;
                        rho[jcomp_rho+32]+=t3_meam[itype]*t3_meam[itype]*rhoa0i;
                    }
                    
                    rho[icomp_rho+13]+=rhoa0i;
                    rho[jcomp_rho+13]+=rhoa0j;
                    
                    
                    
                    
                    
                    iArho1_comp=icomp_rho+4;
                    jArho1_comp=jcomp_rho+4;
                    
                    iArho2_comp=icomp_rho+7;
                    jArho2_comp=jcomp_rho+7;
                    
                    iArho3_comp=icomp_rho+14;
                    jArho3_comp=jcomp_rho+14;
                    
                    iArhob3_comp=icomp_rho+24;
                    jArhob3_comp=jcomp_rho+24;
                    
                    A1j=rhoa1j/rij;
                    B1j=rhoa3j/rij;
                    A2j=rhoa2j/rij2;
                    A3j=rhoa3j/(rij2*rij);
                    
                    A1i=rhoa1i/rij;
                    B1i=rhoa3i/rij;
                    A2i=rhoa2i/rij2;
                    A3i=rhoa3i/(rij2*rij);
                    
                    
                    for(int i=0;i<3;i++)
                    {
                        rho[iArho1_comp++]+=delij[i]*A1j;
                        rho[jArho1_comp++]-=delij[i]*A1i;
                        rho[iArho3_comp++]+=delij[i]*B1j;
                        rho[jArho3_comp++]-=delij[i]*B1i;
                        
                        for(int j=i;j<3;j++)
                        {
                            rho[iArho2_comp++]+=delij[i]*delij[j]*A2j;
                            rho[jArho2_comp++]+=delij[i]*delij[j]*A2i;
                            for(int k=j;k<3;k++)
                            {
                                rho[iArhob3_comp++]+=delij[i]*delij[j]*delij[k]*A3j;
                                rho[jArhob3_comp++]-=delij[i]*delij[j]*delij[k]*A3i;
                            }
                        }
                    }
                }
                
                istart++;
            }
        }
        
        rho[icomp_rho+1]=0.0;
        rho[icomp_rho+2]=-rho[icomp_rho+13]*rho[icomp_rho+13]/3.0;
        rho[icomp_rho+3]=0.0;
        for(int i=0;i<3;i++)
        {
            rho[icomp_rho+1]+=rho[icomp_rho+4+i]*rho[icomp_rho+4+i];
            rho[icomp_rho+3]-=0.6*rho[icomp_rho+24+i]*rho[icomp_rho+24+i];
        }
        
        for(int i=0;i<6;i++)
        {
            rho[icomp_rho+2]+=v2d[i]*rho[icomp_rho+7+i]*rho[icomp_rho+7+i];
        }
        
        for(int i=0;i<10;i++)
        {
            rho[icomp_rho+3]+=v3d[i]*rho[icomp_rho+14+i]*rho[icomp_rho+14+i];
        }
        
        if(rho[icomp_rho]>0.0)
        {
            if(ialloy==1)
            {
                rho[icomp_rho+27]=rho[icomp_rho+27]/rho[icomp_rho+30];
                rho[icomp_rho+28]=rho[icomp_rho+28]/rho[icomp_rho+31];
                rho[icomp_rho+29]=rho[icomp_rho+29]/rho[icomp_rho+32];
            }
            else if(ialloy==2)
            {
                rho[icomp_rho+27]=t1_meam[itype];
                rho[icomp_rho+28]=t2_meam[itype];
                rho[icomp_rho+29]=t3_meam[itype];
            }
            else
            {
                rho[icomp_rho+27]=rho[icomp_rho+27]/rho[icomp_rho];
                rho[icomp_rho+28]=rho[icomp_rho+28]/rho[icomp_rho];
                rho[icomp_rho+29]=rho[icomp_rho+29]/rho[icomp_rho];
            }
        }
        
        
        rho[icomp_rho+33]=rho[icomp_rho+27]*rho[icomp_rho+1]
        +rho[icomp_rho+28]*rho[icomp_rho+3]+rho[icomp_rho+29]*rho[icomp_rho+3];
        
        if(rho[icomp_rho]>0.0)
            rho[icomp_rho+33]=rho[icomp_rho+33]/(rho[icomp_rho]*rho[icomp_rho]);
        
        
        
        Z=Z_meam[itype];
        
        G_gam(rho[icomp_rho+33],ibar_meam[itype],gsmooth_factor,G);
        
        if(ibar_meam[itype]<=0)
            Gbar=1.0;
        else
        {
            get_shpfcn(s,lattice[itype][itype]);
            if(mix_ref_t==1)
                gam=(rho[icomp_rho+27]*s[0]+rho[icomp_rho+28]*s[1]+rho[icomp_rho+29]*s[2])/(Z*Z);
            else
                gam=(t1_meam[itype]*s[0]+t2_meam[itype]*s[1]+t3_meam[itype]*s[2])/(Z*Z);
            
            G_gam(gam,ibar_meam[itype],gsmooth_factor,Gbar);
        }
        rho_vec[iatm]=rho[icomp_rho]*G;
        if(mix_ref_t==1)
        {
            if(ibar_meam[itype]<=0.0)
            {
                Gbar=1.0;
                dGbar=0.0;
            }
            else
            {
                get_shpfcn(s,lattice[itype][itype]);
                gam=(rho[icomp_rho+27]*s[0]+rho[icomp_rho+28]*s[1]+rho[icomp_rho+29]*s[2])/(Z*Z);
                dG_gam(gam,ibar_meam[itype],gsmooth_factor,Gbar,dGbar);
            }
            rho_bkgd=rho0_meam[itype]*Z*Gbar;
        }
        else
        {
            if(bkgd_dyn==1)
            {
                rho_bkgd=rho0_meam[itype]*Z;
            }
            else
                rho_bkgd=rho_ref_meam[itype];
        }
        
        rhob=rho_vec[iatm]/rho_bkgd;
        denom=1.0/rho_bkgd;
        dG_gam(rho[icomp_rho+33],ibar_meam[itype],gsmooth_factor,G,dG);
        rho[icomp_rho+34]=(G-2.0*dG*rho[icomp_rho+33])*denom;
        
        if(rho[icomp_rho]!=0)
            rho[icomp_rho+35]=(dG/rho[icomp_rho])*denom;
        else
            rho[icomp_rho+35]=0.0;
        
        if(mix_ref_t==1)
            rho[icomp_rho+36]=dG/rho[icomp_rho]*denom;
        else
            rho[icomp_rho+36]=0.0;
        
        B=A_meam[itype]*Ec_meam[itype][itype];
        
        if(rhob!=0.0)
        {
            if(emb_lin_neg==1 && rhob<=0.0)
                rho[icomp_rho+37]=-B;
            else
                rho[icomp_rho+37]=B*(log(rhob)+1.0);
        }
        else
        {
            if(emb_lin_neg==1)
                rho[icomp_rho+37]=-B;
            else
                rho[icomp_rho+37]=B;
        }
        
    }
    
    delete [] delij;
    delete [] s;

}

/*--------------------------------------------
 claculate the forces
 --------------------------------------------*/
void ForceField_meam::ff()
{
    int* my_list;
    int my_list_size;

    int jatm,katm;
    
    TYPE0* f;
    TYPE0 xij,xjk,xki,yij,yjk,yki,zij,zjk,zki;
    TYPE0 rij2,rjk2,rij3,rki2;
    TYPE0 rij;
    int icomp,jcomp,kcomp,itype,jtype,ktype;
    int icomp_rho,jcomp_rho;
    
    
    TYPE0 sij;
    TYPE0 rbound,x_ki,x_jk,a;
    TYPE0 cikj,sikj=0.0,cmax,cmin;
    TYPE0 ai,aj,ro0i,ro0j,rhoa0i,rhoa0j,rhoa1i,rhoa1j
    ,rhoa2i,rhoa2j,rhoa3i,rhoa3j;
    TYPE0 drhoa0i,drhoa0j,drhoa1i,drhoa1j
    ,drhoa2i,drhoa2j,drhoa3i,drhoa3j,recip;
    
    TYPE0 arg1i1,arg1j1,arg1i2,arg1j2,arg1i3,arg1j3,arg3i3,arg3j3,arg;
    int nv2_comp,nv3_comp;
    TYPE0 drho0dr1,drho0dr2,a1,a2,a3,a3a,drho1dr1,drho1dr2,drho2dr1,drho2dr2,drho3dr1,drho3dr2;
    TYPE0 drho0ds1,drho0ds2,drho1ds1,drho1ds2,drho2ds1,drho2ds2,drho3ds1,drho3ds2;
    TYPE0 dt1ds1,dt1ds2,dt2ds1,dt2ds2,dt3ds1,dt3ds2;
    TYPE0 dsij1,dsij2,dcikj1,dcikj2,dfc,force1,force2;
    TYPE0 drhods1=0.0,drhods2=0.0;
    TYPE0 re_meam_inv_i,re_meam_inv_j;
    
    
    TYPE0* delij;
    CREATE1D(delij,3);
    TYPE0* delik;
    CREATE1D(delik,3);
    TYPE0* deljk;
    CREATE1D(deljk,3);
    TYPE0* si;
    CREATE1D(si,3);
    TYPE0* sj;
    CREATE1D(sj,3);
    TYPE0* drho1drm1;
    CREATE1D(drho1drm1,3);
    TYPE0* drho1drm2;
    CREATE1D(drho1drm2,3);
    TYPE0* drho2drm1;
    CREATE1D(drho2drm1,3);
    TYPE0* drho2drm2;
    CREATE1D(drho2drm2,3);
    TYPE0* drho3drm1;
    CREATE1D(drho3drm1,3);
    TYPE0* drho3drm2;
    CREATE1D(drho3drm2,3);
    TYPE0* drhodrm1;
    CREATE1D(drhodrm1,3);
    TYPE0* drhodrm2;
    CREATE1D(drhodrm2,3);
    TYPE0* dUdrijm;
    CREATE1D(dUdrijm,3);
    
    int istart=0;
    for(int iatm=0;iatm<atoms->natms;iatm++)
    {
        my_list=neighbor->neighbor_list[iatm];
        my_list_size=neighbor->neighbor_list_size[iatm];
        icomp=iatm*3;
        icomp_rho=iatm*rho_dim;
        itype=type[iatm];
        
        for(int jn=0;jn<my_list_size;jn++)
        {
            jatm=my_list[jn];
            if(jatm>iatm)
            {
                if(scrfcn[istart]!=0.0)
                {
                    jtype=type[jatm];
                    jcomp=jatm*3;
                    jcomp_rho=jatm*rho_dim;
                    
                    xij=x[icomp]-x[jcomp];
                    yij=x[icomp+1]-x[jcomp+1];
                    zij=x[icomp+2]-x[jcomp+2];
                    delij[0]=-xij;
                    delij[1]=-yij;
                    delij[2]=-zij;
                    rij2=xij*xij+yij*yij+zij*zij;
                    rij=sqrt(rij2);
                    if(rij<rc_meam)
                    {
                        /*
                         ind = eltind(elti,eltj)
                         pp = rij*rdrar + 1.0D0
                         kk = pp
                         kk = min(kk,nrar-1)
                         pp = pp - kk
                         pp = min(pp,1.0D0)
                         phi = ((phirar3(kk,ind)*pp + phirar2(kk,ind))*pp+ phirar1(kk,ind))*pp + phirar(kk,ind)
                         phip = (phirar6(kk,ind)*pp + phirar5(kk,ind))*pp+phirar4(kk,ind)
                         recip = 1.0d0/r
                         
                         if (eflag_either.ne.0) then
                         if (eflag_global.ne.0) eng_vdwl = eng_vdwl + phi*sij
                         if (eflag_atom.ne.0) then
                         eatom(i) = eatom(i) + 0.5*phi*sij
                         eatom(j) = eatom(j) + 0.5*phi*sij
                         endif
                         endif
                         */
                        
                        re_meam_inv_i=1.0/re_meam[itype][itype];
                        ai=rij*re_meam_inv_i-1.0;
                        ro0i=rho0_meam[itype];
                        rhoa0i=ro0i*exp(-beta0_meam[itype]*ai);
                        drhoa0i=-beta0_meam[itype]*re_meam_inv_i*rhoa0i;
                        
                        rhoa1i=ro0i*exp(-beta1_meam[itype]*ai);
                        drhoa1i=-beta1_meam[itype]*re_meam_inv_i*rhoa1i;
                        
                        rhoa2i=ro0i*exp(-beta2_meam[itype]*ai);
                        drhoa2i=-beta2_meam[itype]*re_meam_inv_i*rhoa2i;
                        
                        rhoa3i=ro0i*exp(-beta3_meam[itype]*ai);
                        drhoa3i=-beta3_meam[itype]*re_meam_inv_i*rhoa3i;
                        
                        if(itype!=jtype)
                        {
                            
                            re_meam_inv_j=1.0/re_meam[jtype][jtype];
                            aj=rij*re_meam_inv_j-1.0;
                            ro0j=rho0_meam[jtype];
                            rhoa0j=ro0j*exp(-beta0_meam[jtype]*aj);
                            drhoa0j=-beta0_meam[jtype]*re_meam_inv_j*rhoa0j;
                            
                            rhoa1j=ro0j*exp(-beta1_meam[jtype]*aj);
                            drhoa1j=-beta1_meam[jtype]*re_meam_inv_j*rhoa1j;
                            
                            rhoa2j=ro0j*exp(-beta2_meam[jtype]*aj);
                            drhoa2j=-beta2_meam[jtype]*re_meam_inv_j*rhoa2j;
                            
                            rhoa3j=ro0j*exp(-beta3_meam[jtype]*aj);
                            drhoa3j=-beta3_meam[jtype]*re_meam_inv_j*rhoa3j;
                        }
                        else
                        {
                            rhoa0j=rhoa0i;
                            drhoa0j=drhoa0i;
                            rhoa1j=rhoa1i;
                            drhoa1j=drhoa1i;
                            rhoa2j=rhoa2i;
                            drhoa2j=drhoa2i;
                            rhoa3j=rhoa3i;
                            drhoa3j=drhoa3i;
                        }
                        
                        if(ialloy==1)
                        {
                            rhoa1j*=t1_meam[jtype];
                            rhoa2j*=t2_meam[jtype];
                            rhoa3j*=t3_meam[jtype];
                            rhoa1i*=t1_meam[itype];
                            rhoa2i*=t2_meam[itype];
                            rhoa3i*=t3_meam[itype];
                            drhoa1j*=t1_meam[jtype];
                            drhoa2j*=t2_meam[jtype];
                            drhoa3j*=t3_meam[jtype];
                            drhoa1i*=t1_meam[itype];
                            drhoa2i*=t2_meam[itype];
                            drhoa3i*=t3_meam[itype];
                        }
                        
                        
                        
                        
                        arg1i1=0.0;
                        arg1j1=0.0;
                        arg1i2=0.0;
                        arg1j2=0.0;
                        arg1i3=0.0;
                        arg1j3=0.0;
                        arg3i3=0.0;
                        arg3j3=0.0;
                        
                        nv2_comp=0;
                        nv3_comp=0;
                        for(int i=0;i<3;i++)
                        {
                            for(int j=i;j<3;j++)
                            {
                                for(int k=j;k<3;k++)
                                {
                                    arg=delij[i]*delij[j]*delij[k]*v3d[nv3_comp];
                                    arg1i3+=rho[icomp_rho+14+nv3_comp]*arg;
                                    arg1j3-=rho[jcomp_rho+14+nv3_comp]*arg;
                                    nv3_comp++;
                                }
                                arg=delij[i]*delij[j]*v2d[nv2_comp];
                                arg1i2+=rho[icomp_rho+7+nv2_comp];
                                arg1j2+=rho[jcomp_rho+7+nv2_comp];
                                nv2_comp++;
                            }
                            arg1i1+=rho[icomp_rho+4+i]*delij[i];
                            arg1j1-=rho[jcomp_rho+4+i]*delij[i];
                            arg3i3+=rho[icomp_rho+24+i]*delij[i];
                            arg3j3-=rho[jcomp_rho+24+i]*delij[i];
                        }
                        
                        
                        sij = scrfcn[istart]*fcpair[istart];
                        
                        //rho0 terms
                        drho0dr1=drhoa0j*sij;
                        drho0dr2=drhoa0i*sij;
                        
                        //rho1 terms
                        a1=2.0*sij/rij;
                        drho1dr1=a1*(drhoa1j-rhoa1j/rij)*arg1i1;
                        drho1dr2=a1*(drhoa1i-rhoa1i/rij)*arg1j1;
                        
                        
                        for(int i=0;i<3;i++)
                        {
                            drho1drm1[i]=a1*rhoa1j*rho[icomp_rho+4+i];
                            drho1drm2[i]=-a1*rhoa1i*rho[jcomp_rho+4+i];
                        }
                        
                        //rho2 terms
                        a2=2.0*sij/rij2;
                        drho2dr1=a2*(drhoa2j-2.0*rhoa2j/rij)*arg1i2-2.0/3.0*rho[icomp_rho+13]*drhoa2j*sij;
                        drho2dr2=a2*(drhoa2i-2.0*rhoa2i/rij)*arg1j2-2.0/3.0*rho[jcomp_rho+13]*drhoa2i*sij;
                        a2=4.0*sij/rij2;
                        
                        
                        
                        for(int i=0;i<3;i++)
                        {
                            drho2drm1[i]=0.0;
                            drho2drm2[i]=0.0;
                            for(int j=0;j<3;j++)
                            {
                                drho2drm1[i]+=rho[icomp_rho+7+vind2d[i][j]]*delij[j];
                                drho2drm1[i]-=rho[jcomp_rho+7+vind2d[i][j]]*delij[j];
                            }
                            drho2drm1[i]*=a2*rhoa2j;
                            drho2drm2[i]*=-a2*rhoa2i;
                        }
                        
                        rij3=rij*rij2;
                        a3=2.0*sij/rij3;
                        a3a=1.2*sij/rij;
                        drho3dr1 = a3*(drhoa3j - 3*rhoa3j/rij)*arg1i3- a3a*(drhoa3j - rhoa3j/rij)*arg3i3;
                        drho3dr2 = a3*(drhoa3i - 3*rhoa3i/rij)*arg1j3- a3a*(drhoa3i - rhoa3i/rij)*arg3j3;
                        a3 = 6.0*sij/rij3;
                        a3a = 1.2*sij/rij;
                        
                        
                        for(int i=0;i<3;i++)
                        {
                            drho3drm1[i]=0.0;
                            drho3drm2[i]=0.0;
                            nv2_comp=0;
                            for(int j=0;j<3;j++)
                            {
                                for(int k=j;k<3;k++)
                                {
                                    arg = delij[j]*delij[k]*v2d[nv2_comp];
                                    drho3drm1[i]+=rho[icomp_rho+14+vind3d[i][j][k]]*arg;
                                    drho3drm2[i]+=rho[jcomp_rho+14+vind3d[i][j][k]]*arg;
                                    nv2_comp++;
                                }
                            }
                            drho3drm1[i]+=(a3*drho3drm1[i]-a3a*rho[icomp_rho+24+i])*rhoa3j;
                            drho3drm2[i]-=(a3*drho3drm1[i]-a3a*rho[jcomp_rho+24+i])*rhoa3i;
                        }
                        
                        TYPE0 t1i,t2i,t3i,t1j,t2j,t3j;
                        t1i=rho[icomp_rho+27];
                        t2i=rho[icomp_rho+28];
                        t3i=rho[icomp_rho+29];
                        t1j=rho[jcomp_rho+27];
                        t2j=rho[jcomp_rho+28];
                        t3j=rho[jcomp_rho+29];
                        
                        TYPE0 a1i,a2i,a3i,a1j,a2j,a3j,dt1dr1,dt1dr2,dt2dr1,dt2dr2,dt3dr1,dt3dr2;
                        
                        if(ialloy==1)
                        {
                            a1i=0.0;
                            a2i=0.0;
                            a3i=0.0;
                            a1j=0.0;
                            a2j=0.0;
                            a3j=0.0;
                            
                            if(rho[icomp_rho+30]!=0.0)
                                a1i=drhoa0j*sij/rho[icomp_rho+30];
                            if(rho[jcomp_rho+30]!=0.0)
                                a1j=drhoa0i*sij/rho[jcomp_rho+30];
                            if(rho[icomp_rho+31]!=0.0)
                                a2i=drhoa0j*sij/rho[icomp_rho+31];
                            if(rho[jcomp_rho+31]!=0.0)
                                a2j=drhoa0i*sij/rho[jcomp_rho+31];
                            if(rho[icomp_rho+32]!=0.0)
                                a3i=drhoa0j*sij/rho[icomp_rho+32];
                            if(rho[jcomp_rho+32]!=0.0)
                                a3j=drhoa0i*sij/rho[jcomp_rho+32];
                            
                            dt1dr1=a1i*(t1_meam[jtype]-t1i*t1_meam[jtype]*t1_meam[jtype]);
                            dt1dr2=a1j*(t1_meam[itype]-t1j*t1_meam[itype]*t1_meam[itype]);
                            dt2dr1=a2i*(t2_meam[jtype]-t1i*t2_meam[jtype]*t2_meam[jtype]);
                            dt2dr2=a2j*(t2_meam[itype]-t1j*t2_meam[itype]*t2_meam[itype]);
                            dt3dr1=a3i*(t3_meam[jtype]-t1i*t3_meam[jtype]*t3_meam[jtype]);
                            dt3dr2=a3j*(t3_meam[itype]-t1j*t3_meam[itype]*t3_meam[itype]);
                        }
                        else if(ialloy==2)
                        {
                            dt1dr1=0.0;
                            dt1dr2=0.0;
                            dt2dr1=0.0;
                            dt2dr2=0.0;
                            dt3dr1=0.0;
                            dt3dr2=0.0;
                        }
                        else
                        {
                            ai=0.0;
                            if(rho[icomp_rho]!=0.0)
                                ai=drhoa0j*sij/rho[icomp_rho];
                            aj=0.0;
                            if(rho[jcomp_rho]!=0.0)
                                aj=drhoa0i*sij/rho[jcomp_rho];
                            
                            dt1dr1=ai*(t1_meam[jtype]-t1i);
                            dt1dr2 = aj*(t1_meam[itype]-t1j);
                            dt2dr1 = ai*(t2_meam[jtype]-t2i);
                            dt2dr2 = aj*(t2_meam[itype]-t2j);
                            dt3dr1 = ai*(t3_meam[jtype]-t3i);
                            dt3dr2 = aj*(t3_meam[itype]-t3j);
                            
                        }
                        
                        get_shpfcn(si,lattice[itype][jtype]);
                        get_shpfcn(sj,lattice[jtype][itype]);
                        
                        TYPE0 drhodr1,drhodr2;
                        drhodr1=rho[icomp_rho+34]*drho0dr1+rho[icomp_rho+35]
                        *(dt1dr1*rho[icomp_rho+1]+t1i*drho1dr1
                          +dt2dr1*rho[icomp_rho+2]+t2i*drho2dr1
                          +dt3dr1*rho[icomp_rho+3]+t3i*drho3dr1)
                        -rho[icomp_rho+36]*
                        (si[0]*dt1dr1+si[1]*dt2dr1+si[2]*dt3dr1);
                        
                        drhodr2=rho[jcomp_rho+34]*drho0dr2+rho[jcomp_rho+35]
                        *(dt1dr2*rho[jcomp_rho+1]+t1j*drho1dr2
                          +dt2dr2*rho[jcomp_rho+2]+t2j*drho2dr2
                          +dt3dr2*rho[jcomp_rho+3]+t3j*drho3dr2)
                        -rho[jcomp_rho+36]*
                        (sj[0]*dt1dr2+sj[1]*dt2dr2+sj[2]*dt3dr2);
                        
                        for(int i=0;i<3;i++)
                        {
                            drhodrm1[i]=rho[icomp_rho+35]*(t1i*drho1drm1[i]
                                                           +t2i*drho2drm1[i]+t3i*drho3drm1[i]);
                            drhodrm2[i]=rho[jcomp_rho+35]*(t1j*drho1drm2[i]
                                                           +t2j*drho2drm2[i]+t3j*drho3drm2[i]);
                        }
                        
                        
                        if(dscrfcn[istart]!=0.0)
                        {
                            drho0ds1=rhoa0j;
                            drho0ds2=rhoa0i;
                            a1=2.0/rij;
                            drho1ds1=a1*rhoa1j*arg1i1;
                            drho1ds2=a1*rhoa1i*arg1j1;
                            a2=2.0/rij2;
                            
                            drho2ds1=a2*rhoa2j*arg1i2
                            -2.0/3.0*rho[icomp_rho+13]*rhoa2j;
                            drho2ds2=a2*rhoa2i*arg1j2
                            -2.0/3.0*rho[jcomp_rho+13]*rhoa2i;
                            
                            a3=2.0/rij3;
                            a3a=1.2/rij;
                            
                            drho3ds1=a3*rhoa3j*arg1i3-a3a*rhoa3j*arg3i3;
                            drho3ds2=a3*rhoa3i*arg1j3-a3a*rhoa3i*arg3j3;
                            
                            if(ialloy==1)
                            {
                                a1i=0.0;
                                a1j=0.0;
                                a2i=0.0;
                                a2j=0.0;
                                a3i=0.0;
                                a3j=0.0;
                                if(rho[icomp_rho+30]!=0.0)
                                    a1i=drhoa0j*sij/rho[icomp_rho+30];
                                if(rho[jcomp_rho+30]!=0.0)
                                    a1j=drhoa0i*sij/rho[jcomp_rho+30];
                                if(rho[icomp_rho+31]!=0.0)
                                    a2i=drhoa0j*sij/rho[icomp_rho+31];
                                if(rho[jcomp_rho+31]!=0.0)
                                    a2j=drhoa0i*sij/rho[jcomp_rho+31];
                                if(rho[icomp_rho+32]!=0.0)
                                    a3i=drhoa0j*sij/rho[icomp_rho+32];
                                if(rho[jcomp_rho+32]!=0.0)
                                    a3j=drhoa0i*sij/rho[jcomp_rho+32];
                                
                                dt1ds1=a1i*(t1_meam[jtype]-t1i*t1_meam[jtype]*t1_meam[jtype]);
                                dt1ds2=a1j*(t1_meam[itype]-t1j*t1_meam[itype]*t1_meam[itype]);
                                dt2ds1=a2i*(t2_meam[jtype]-t2i*t2_meam[jtype]*t1_meam[jtype]);
                                dt2ds2=a2j*(t2_meam[itype]-t2j*t2_meam[itype]*t1_meam[itype]);
                                dt3ds1=a3i*(t3_meam[jtype]-t3i*t3_meam[jtype]*t1_meam[jtype]);
                                dt3ds2=a3j*(t3_meam[itype]-t3j*t3_meam[itype]*t1_meam[itype]);
                                
                            }
                            else if(ialloy==2)
                            {
                                dt1ds1=0.0;
                                dt1ds2=0.0;
                                dt2ds1=0.0;
                                dt2ds2=0.0;
                                dt3ds1=0.0;
                                dt3ds2=0.0;
                            }
                            else
                            {
                                ai=0.0;
                                if(rho[icomp_rho]!=0.0)
                                    ai=drhoa0j/rho[icomp_rho];
                                aj=0.0;
                                if(rho[jcomp_rho]!=0.0)
                                    aj=drhoa0i/rho[jcomp_rho];
                                
                                dt1ds1=ai*(t1_meam[jtype]-t1i);
                                dt1ds2=aj*(t1_meam[itype]-t1j);
                                dt2ds1=ai*(t2_meam[jtype]-t2i);
                                dt2ds2=aj*(t2_meam[itype]-t2j);
                                dt3ds1=ai*(t3_meam[jtype]-t3i);
                                dt3ds2=aj*(t3_meam[itype]-t3j);
                            }
                            
                            
                            drhods1=rho[icomp_rho+34]*drho0ds1+rho[icomp_rho+35]
                            *(dt1ds1*rho[icomp_rho+1]+t1i*drho1ds1
                              +dt2ds1*rho[icomp_rho+2]+t2i*drho2ds1
                              +dt3ds1*rho[icomp_rho+3]+t3i*drho3ds1)
                            -rho[icomp_rho+36]*
                            (si[0]*dt1ds1+si[1]*dt2ds1+si[2]*dt3ds1);
                            
                            drhods2=rho[jcomp_rho+34]*drho0ds2+rho[jcomp_rho+35]
                            *(dt1ds2*rho[jcomp_rho+1]+t1j*drho1ds2
                              +dt2ds2*rho[jcomp_rho+2]+t2j*drho2ds2
                              +dt3ds2*rho[jcomp_rho+3]+t3j*drho3ds2)
                            -rho[jcomp_rho+36]*
                            (sj[0]*dt1ds2+sj[1]*dt2ds2+sj[2]*dt3ds2);
                        }
                        
                        TYPE0 dUdrij,dUdsij,phi,phip,force,forcem;
                        dUdrij=phip*sij+rho[icomp_rho+37]*drhodr1
                        +rho[jcomp_rho+37]*drhodr2;
                        
                        dUdsij=0.0;
                        if(dscrfcn[istart]!=0.0)
                            dUdsij=phi+rho[icomp_rho+37]*drhods1
                            +rho[jcomp_rho+37]*drhods2;
                        
                        for(int i=0;i<3;i++)
                            dUdrijm[i]=rho[icomp_rho+37]*drhodrm1[i]
                            +rho[jcomp_rho+37]*drhodrm2[i];
                        
                        force=dUdrij*recip+dUdsij*dscrfcn[istart];
                        
                        for(int i=0;i<3;i++)
                        {
                            forcem=delij[i]*force+dUdrijm[i];
                            f[icomp+i]+=forcem;
                            if(jatm<atoms->natms)
                                f[jcomp+i]-=forcem;
                        }
                        
                        /*
                         stress and all the other shit
                         
                         remember to fix it
                         */
                        if(sij!=0.0 && sij!=1.0)
                        {
                            rbound=rij2*ebound_meam[itype][jtype];
                            
                            for(int kn=0;kn<my_list_size;kn++)
                            {
                                katm=my_list[kn];
                                
                                if(katm!=jatm)
                                {
                                    dsij1=dsij2=0.0;
                                    ktype=type[katm];
                                    kcomp=3*katm;
                                    
                                    xjk=x[jcomp]-x[kcomp];
                                    yjk=x[jcomp+1]-x[kcomp+1];
                                    zjk=x[jcomp+2]-x[kcomp+2];
                                    rjk2=xjk*xjk+yjk*yjk+zjk*zjk;
                                    
                                    if(rjk2<rbound)
                                    {
                                        xki=x[kcomp]-x[icomp];
                                        yki=x[kcomp+1]-x[icomp+1];
                                        zki=x[kcomp+2]-x[icomp+2];
                                        rki2=xki*xki+yki*yki+zki*zki;
                                        
                                        if(rki2<rbound)
                                        {
                                            x_ki=rki2/rij2;
                                            x_jk=rjk2/rij2;
                                            a=1.0-(x_ki-x_jk)*(x_ki-x_jk);
                                            if(a!=0.0)
                                            {
                                                cikj=(2.0*(x_ki+x_jk)+a-2.0)/a;
                                                cmax=c_max[itype][jtype][ktype];
                                                cmin=c_min[itype][jtype][ktype];
                                                if(cikj<cmax&&cikj>cmin)
                                                {
                                                    cikj=(cikj-cmin)/(cmax-cmin);
                                                    dfcut(cikj,sikj,dfc);
                                                    dCfunc2(rij2,rki2,rjk2,dcikj1,dcikj2);
                                                    a=(sij/(cmax-cmin))*(dfc/sikj);
                                                    dsij1=a*dcikj1;
                                                    dsij2=a*dcikj2;
                                                }
                                            }
                                        }
                                        
                                        if(dsij1!=0.0 || dsij2!=0.0)
                                        {
                                            force1=dUdsij*dsij1;
                                            force2=dUdsij*dsij2;
                                            delik[0]=xki;
                                            deljk[0]=-xjk;
                                            delik[1]=yki;
                                            deljk[1]=-yjk;
                                            delik[2]=zki;
                                            deljk[2]=-zjk;
                                            for(int i=0;i<3;i++)
                                                f[icomp+i]+=force1*delik[i];
                                            if(jatm<atoms->natms)
                                                for(int i=0;i<3;i++)
                                                    f[jcomp+i]+=force2*deljk[i];
                                            
                                            if(katm<atoms->natms)
                                                for(int i=0;i<3;i++)
                                                    f[kcomp+i]-=force1*delik[i]+force2*deljk[i];
                                            
                                            /*
                                             the stress contrbution crap
                                             remember to set it
                                             */
                                            
                                        }
                                    }
                                    
                                }
                                
                            }
                        }
                        
                    }
                }
                
                istart++;
            }
        }
    }
    
    
    delete [] delik;
    delete [] deljk;
    delete [] dUdrijm;
    delete [] drho3drm1;
    delete [] drho3drm2;
    delete [] drho1drm1;
    delete [] drho1drm2;
    delete [] drho2drm1;
    delete [] drho2drm2;
    delete [] delij;
    delete [] si;
    delete [] sj;
    
}
/*--------------------------------------------
 G_gam
 --------------------------------------------*/
void ForceField_meam::G_gam(TYPE0 Gamma,int ibar
,TYPE0 gsmooth_factor,TYPE0& G)
{
    TYPE0 gsmooth_switchpoint;
    if(ibar==0||ibar==4)
    {
        gsmooth_switchpoint=-gsmooth_factor/(gsmooth_factor+1.0);
        if(Gamma<gsmooth_switchpoint)
        {
            G=pow(gsmooth_switchpoint/Gamma,gsmooth_factor)
            /(gsmooth_factor+1.0);
            G=sqrt(G);
        }
        else
            G=sqrt(1.0+Gamma);
    }
    else if (ibar==1)
    {
        G=exp(Gamma*0.5);
    }
    else if (ibar==3)
    {
        G=2.0/(1.0+exp(-Gamma));
    }
    else if (ibar==-5)
    {
        if ((1.0+Gamma)>=0.0)
        {
            G = sqrt(1.0+Gamma);
        }
        else
        {
            G=-sqrt(-1.0-Gamma);
        }
    }
    else
        error->abort("ibar error in MEAM");
}
/*--------------------------------------------
 dG_gam
 --------------------------------------------*/
void ForceField_meam::dG_gam(TYPE0 Gamma,int ibar
,TYPE0 gsmooth_factor,TYPE0& G,TYPE0& dG)
{
    TYPE0 gsmooth_switchpoint;
    if(ibar==0||ibar==4)
    {
        gsmooth_switchpoint=-gsmooth_factor/(gsmooth_factor+1.0);
        if(Gamma<gsmooth_switchpoint)
        {
            G=pow(gsmooth_switchpoint/Gamma,gsmooth_factor)
            /(gsmooth_factor+1.0);
            G=sqrt(G);
            dG=0.5/G;
        }
        else
        {
            G=sqrt(1.0+Gamma);
            dG=0.5/G;
        }
    }
    else if (ibar==1)
    {
        G=exp(Gamma*0.5);
        dG=0.5*G;
    }
    else if (ibar==3)
    {
        G=2.0/(1.0+exp(-Gamma));
        dG=G*(1.0-0.5*G);
    }
    else if (ibar==-5)
    {
        if ((1.0+Gamma)>=0.0)
        {
            G=sqrt(1.0+Gamma);
            dG=0.5/G;
        }
        else
        {
            G=-sqrt(-1.0-Gamma);
            dG=-0.5/G;
        }
    }
    else
        error->abort("ibar error in MEAM");
}
/*--------------------------------------------
 fcut
 --------------------------------------------*/
void ForceField_meam::fcut(TYPE0 xi,TYPE0& fc)
{
    if(xi>1.0)
    {
        fc=1.0;
    }
    else if (xi<0.0)
    {
        fc=0.0;
    }
    else
    {
        TYPE0 a=1.0-xi;
        a*=a;
        a*=a;
        a=1.0-a;
        fc=a*a;
    }
}
/*--------------------------------------------
 dfcut
 --------------------------------------------*/
void ForceField_meam::dfcut(TYPE0 xi,TYPE0& fc,TYPE0& dfc)
{
    if(xi>1.0)
    {
        fc=1.0;
        dfc=0.0;
    }
    else if (xi<0.0)
    {
        fc=0.0;
        dfc=0.0;
    }
    else
    {
        TYPE0 a=1.0-xi;
        TYPE0 a3=a*a*a;
        TYPE0 a4=a*a*a*a;
        fc=(1.0-a4)*(1.0-a4);
        dfc=8*(1.0-a4)*a3;
    }
}
/*--------------------------------------------
 dCfunc
 --------------------------------------------*/
void ForceField_meam::dCfunc(TYPE0 rij2
,TYPE0 rjk2,TYPE0 rki2,TYPE0& dCikj)
{
    TYPE0 rij4=rij2*rij2;
    TYPE0 a=rki2-rjk2;
    TYPE0 b=rki2+rjk2;
    TYPE0 denom=rij4-a*a;
    denom *= denom;
    dCikj=-4.0*(-2.0*rij2*a*a +rij4*b+a*a*b)/denom;
}
/*--------------------------------------------
 dCfunc2
 --------------------------------------------*/
void ForceField_meam::dCfunc2(TYPE0 rij2
,TYPE0 rjk2,TYPE0 rki2,TYPE0& dCikj1,TYPE0& dCikj2)
{
    TYPE0 rij4=rij2*rij2;
    TYPE0 rjk4=rjk2*rjk2;
    TYPE0 rki4=rki2*rki2;
    TYPE0 a=rki2-rjk2;
    TYPE0 denom=rij4-a*a;
    denom *= denom;
    
    dCikj1=4.0*rij2*(rij4+rki4+2.0*rki2*rjk2-3.0*rjk4-2.0*rij2*a)/denom;
    dCikj2=4.0*rij2*(rij4-3.0*rki4+2.0*rki2*rjk2+rjk4+2.0*rij2*a)/denom;
}
/*--------------------------------------------
 dCfunc2
 --------------------------------------------*/
void ForceField_meam::get_shpfcn(TYPE0* s,int latt)
{
    if(latt==FCC ||
       latt==BCC ||
       latt==B1  ||
       latt==B2)
    {
        s[0]=0.0;
        s[1]=0.0;
        s[2]=0.0;
    }
    else if(latt==HCP)
    {
        s[0]=0.0;
        s[1]=0.0;
        s[2]=1.0/3.0;
    }
    else if(latt==DIAMOND)
    {
        s[0]=0.0;
        s[1]=0.0;
        s[2]=32.0/9.0;
    }
    else if(latt==DIM)
    {
        s[0]=1.0;
        s[1]=2.0/3.0;
        s[2]=0.4;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void ForceField_meam::init()
{
}
/*--------------------------------------------
 run
 --------------------------------------------*/
void ForceField_meam::force_calc
(int st_clc,TYPE0* en_st)
{
}
/*--------------------------------------------
 fin
 --------------------------------------------*/
void ForceField_meam::fin()
{
}
/*--------------------------------------------
 only energy calculation this is useful for
 minimization/linesearch methods that do not
 use derivatives of energy
 --------------------------------------------*/
TYPE0 ForceField_meam::energy_calc()
{
    return 0.0;
}


