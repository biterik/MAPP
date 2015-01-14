
#include "neighbor.h"
#include "ff_eam.h"
#include "atom_types.h"
using namespace MAPP_NS;

enum{NOT_SET,FUNC_FL,SET_FL,FINNIS_FL};

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam::
ForceField_eam(MAPP* mapp) : ForceField(mapp)
{
    if(mapp->mode!=MD)
        error->abort("this forcefield works only with md mode");

    no_types=atom_types->no_types;
    
    allocated=0;
    eam_mode=NOT_SET;
    max_pairs=0;

    int no_types=atom_types->no_types;
    CREATE1D(nrgy_strss,7);
    CREATE1D(cut_sk_sq,no_types*(no_types+1)/2);
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam::~ForceField_eam()
{
    delete [] nrgy_strss;
    int no_types=atom_types->no_types;
    if(no_types)
        delete [] cut_sk_sq;
    
    if(allocated) clean_up();
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam::
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
    /*
    TYPE0* dF;
    atoms->vectors[dF_n].ret(dF);
     */
    int* type;
    atoms->vectors[type_n].ret(type);
    
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,z2p,z2;
    TYPE0 r,p,r_inv,fpair,tmp0,tmp1;
    TYPE0 drho_i_dr,drho_j_dr,dphi_dr;
    TYPE0 rho_i,rho_j,phi;
    int m,istart;
    TYPE0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
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
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            drhoi_dr[istart]=drhoj_dr[istart]=0.0;
            if(rsq < cut_sq)
            {
                r=sqrt(rsq);
                r_inv=1.0/r;
                p=r*dr_inv;
                m=static_cast<int>(p);
                m=MIN(m,nr-2);
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
            
                rho[iatm]+=rho_i;
                fpair=-dphi_dr*r_inv;
                f[icomp]+=fpair*dx0;
                f[icomp+1]+=fpair*dx1;
                f[icomp+2]+=fpair*dx2;
                if(jatm<natms)
                {
                    rho[jatm]+=rho_j;
                    f[jcomp]-=fpair*dx0;
                    f[jcomp+1]-=fpair*dx1;
                    f[jcomp+2]-=fpair*dx2;
                    
                    nrgy_strss[0]+=phi;
                    if (st_clc)
                    {
                        nrgy_strss[1]-=fpair*dx0*dx0;
                        nrgy_strss[2]-=fpair*dx1*dx1;
                        nrgy_strss[3]-=fpair*dx2*dx2;
                        nrgy_strss[4]-=fpair*dx1*dx2;
                        nrgy_strss[5]-=fpair*dx2*dx0;
                        nrgy_strss[6]-=fpair*dx0*dx1;
                    }
                }
                else
                {
                    nrgy_strss[0]+=0.5*phi;
                    if (st_clc)
                    {
                        nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                        nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                        nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                        nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                        nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                        nrgy_strss[6]-=0.5*fpair*dx0*dx1;
                    }
                }
                
                drhoi_dr[istart]=-drho_i_dr*r_inv;
                drhoj_dr[istart]=-drho_j_dr*r_inv;
            }
            
            istart++;
        }
        p=rho[iatm]*drho_inv;
        m=static_cast<int> (p);
        m=MIN(m,nr-2);
        p-=m;
        p=MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        tmp1=(coef[6]*p+coef[5])*p+coef[4];
        tmp0=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
        if(rho[iatm]>rho_max)
            tmp0+=tmp1*(rho[iatm]-rho_max);
        nrgy_strss[0]+=tmp0;
        rho[iatm]=tmp1;
        
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
                jtype=type[jatm];
                jcomp=3*jatm;
                
                fpair=rho[iatm]*drhoi_dr[istart]+rho[jatm]*drhoj_dr[istart];
                
                dx0=x[icomp]-x[jcomp];
                dx1=x[icomp+1]-x[jcomp+1];
                dx2=x[icomp+2]-x[jcomp+2];
                
                f[icomp]+=dx0*fpair;
                f[icomp+1]+=dx1*fpair;
                f[icomp+2]+=dx2*fpair;
                
                if(jatm<natms)
                {
                    f[jcomp]-=dx0*fpair;
                    f[jcomp+1]-=dx1*fpair;
                    f[jcomp+2]-=dx2*fpair;
                    
                    if (st_clc)
                    {
                        nrgy_strss[1]-=fpair*dx0*dx0;
                        nrgy_strss[2]-=fpair*dx1*dx1;
                        nrgy_strss[3]-=fpair*dx2*dx2;
                        nrgy_strss[4]-=fpair*dx1*dx2;
                        nrgy_strss[5]-=fpair*dx2*dx0;
                        nrgy_strss[6]-=fpair*dx0*dx1;
                    }
                }
                else
                {
                    if (st_clc)
                    {
                        nrgy_strss[1]-=0.5*fpair*dx0*dx0;
                        nrgy_strss[2]-=0.5*fpair*dx1*dx1;
                        nrgy_strss[3]-=0.5*fpair*dx2*dx2;
                        nrgy_strss[4]-=0.5*fpair*dx1*dx2;
                        nrgy_strss[5]-=0.5*fpair*dx2*dx0;
                        nrgy_strss[6]-=0.5*fpair*dx0*dx1;
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
 energy calculation
 --------------------------------------------*/
TYPE0 ForceField_eam::energy_calc()
{
    TYPE0* x;
    atoms->vectors[x_n].ret(x);
    TYPE0* f;
    atoms->vectors[f_n].ret(f);
    TYPE0* rho;
    atoms->vectors[rho_n].ret(rho);
    int* type;
    atoms->vectors[type_n].ret(type);
    
    int iatm,jatm;
    
    int itype,jtype,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq;
    TYPE0 r,p,phi,tmp0;
    int m;
    TYPE0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    TYPE0 en=0.0;
    TYPE0 en_tot=0.0;
    
    int natms=atoms->natms;
    
    for(iatm=0;iatm<natms;iatm++) rho[iatm]=0.0;
    
    for(iatm=0;iatm<natms;iatm++)
    {
        itype=type[iatm];
        icomp=3*iatm;
        for(int j=0;j<neighbor_list_size[iatm];j++)
        {
            jatm=neighbor_list[iatm][j];
            jtype=type[jatm];
            
            jcomp=3*jatm;
            dx0=x[icomp]-x[jcomp];
            dx1=x[icomp+1]-x[jcomp+1];
            dx2=x[icomp+2]-x[jcomp+2];
            rsq=dx0*dx0+dx1*dx1+dx2*dx2;
            
            if(rsq<cut_sq)
            {
                r=sqrt(rsq);
                
                p=r*dr_inv;
                m=static_cast<int>(p);
                m=MIN(m,nr-2);
                p-=m;
                p=MIN(p,1.0);

                coef=phi_r_arr[type2phi[itype][jtype]][m];
                phi=(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                
                coef=rho_arr[type2rho[jtype][itype]][m];
                rho[iatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(jatm<natms)
                {
                    coef=rho_arr[type2rho[itype][jtype]][m];
                    rho[jatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                    en+=phi;
                }
                else
                    en+=0.5*phi;
                

            }
        }
        
        p=rho[iatm]*drho_inv;
        m=static_cast<int>(p);
        m=MIN(m,nr-2);
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
    return en_tot;
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_eam::init()
{
    TYPE0 skin=atoms->skin;
    TYPE0 ph_cut=0.0;
    int no_types=atom_types->no_types;
    for (int i=0;i<no_types*(no_types+1)/2;i++)
        cut_sk_sq[i]=cut_sq+(skin)*(skin)
            +2*sqrt(cut_sq)*(skin);
    
    ph_cut=sqrt(cut_sq);
    
    atoms->set_ph(ph_cut);
    
    x_n=atoms->find("x");
    f_n=atoms->find("f");
    type_n=atoms->find("type");
    
    neighbor->pair_wise=1;
    
    rho_n=atoms->add<TYPE0>(1,1,"rho");
    //dF_n=atoms->add<TYPE0>(1,1,"dF");
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_eam::fin()
{
    if(max_pairs)
    {
        delete [] drhoi_dr;
        delete [] drhoj_dr;
        max_pairs=0;
    }
    
    //atoms->del(dF_n);
    atoms->del(rho_n);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_eam::coef(int narg,char** arg)
{
    if (narg<3)
        error->abort("wrong coeff command "
            "for eam Force Field");
    clean_up();
    if(strcmp(arg[1],"FS")==0)
    {
        eam_mode=FINNIS_FL;
        set_fs(narg-2,&arg[2]);
    }
    else if(strcmp(arg[1],"SetFL")==0)
    {
        eam_mode=SET_FL;
        set_setfl(narg-2,&arg[2]);
    }
    else if(strcmp(arg[1],"FuncFL")==0)
    {
        eam_mode=FUNC_FL;
        set_funcfl(narg-2,&arg[2]);
    }
    else
        error->abort("wrong coeff command "
            "for eam Force Field");
    cut_sq=(nr-1.0)*(nr-1.0)*dr*dr;
    rho_max=(nrho-1.0)*drho;
    
    set_arrays();

}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam::set_funcfl(int no_files
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
    TYPE0* tmp;
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    int narg;
    char** arg;
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
        
        int tot=nrhos[ityp]+2*nrs[ityp];
        
        CREATE1D(tmp,tot);
        
        int ipos=0;
        while (ipos<tot)
        {
            if(line_read(fp,line)==-1)
                error->abort("eam potential file ended immaturely");
            
            narg=mapp->parse_line(line,arg);
            
            if(ipos+narg>tot)
                error->abort("eam potential file ended immaturely");
            
            for(int i=0;i<narg;i++)
            {
                tmp[ipos]=atof(arg[i]);
                ipos++;
                delete [] arg[i];
            }
            if(narg)
                delete [] arg;
            
        }
        
        if(atoms->my_p_no==0)
            fclose(fp);
        
        CREATE1D(tmp_F[ityp],nrhos[ityp]+1);
        CREATE1D(tmp_zi[ityp],nrs[ityp]+1);
        CREATE1D(tmp_rho[ityp],nrs[ityp]+1);
        
        tmp_rho[ityp][0]=tmp_zi[ityp][0]=tmp_F[ityp][0]=0.0;
        memcpy(&tmp_F[ityp][1],&tmp[0],nrhos[ityp]*sizeof(TYPE0));
        memcpy(&tmp_zi[ityp][1],&tmp[nrhos[ityp]],nrs[ityp]*sizeof(TYPE0));
        memcpy(&tmp_rho[ityp][1],&tmp[nrhos[ityp]+nrs[ityp]],nrs[ityp]*sizeof(TYPE0));
        
        delete [] tmp;
        
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
void ForceField_eam::set_setfl(int no_files
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
    

    
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    delete [] arg;
    
    
    if(line_read(fp,line)==-1)
        error->abort("eam potential file ended immaturely");
    
    
    int ipos=0;
    int tot=tot_no_types*(nrho+nr)+tot_no_types*(tot_no_types+1)/2*nr;
    
    TYPE0* tmp;
    CREATE1D(tmp,tot);
    while (ipos<tot)
    {
        if(line_read(fp,line)==-1)
            error->abort("eam potential file ended immaturely");
        
        narg=mapp->parse_line(line,arg);
        
        if(ipos+narg>tot)
            error->abort("eam potential file ended immaturely");
        
        for(int i=0;i<narg;i++)
        {
            tmp[ipos]=atof(arg[i]);
            ipos++;
            delete [] arg[i];
        }
        if(narg)
            delete [] arg;
        
    }
    
    delete [] line;
    if(atoms->my_p_no==0)
        fclose(fp);
    
    
    allocate();
    
    
    int icur_pos,jcur_pos;
    int component;
    icur_pos=0;
    ipos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ityp==type_ref[icur_pos][0])
        {
            component=type_ref[icur_pos][1];
            
            for(int i=0;i<nrho;i++)
                F_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nrho;
            
            for(int i=0;i<nr;i++)
                rho_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nr;
            
            icur_pos++;
        }
        else
        {
            ipos+=nrho+nr;
        }
        
    }
    
    icur_pos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ityp==type_ref[icur_pos][0])
        {
            jcur_pos=0;
            for(int jtyp=0;jtyp<ityp+1;jtyp++)
            {
                if(jtyp==type_ref[jcur_pos][0])
                {
                    component=type2phi[type_ref[icur_pos][1]][type_ref[jcur_pos][1]];
                    
                    for(int i=0;i<nr;i++)
                        phi_r_arr[component][i][0]=tmp[ipos+i];
                    
                    ipos+=nr;
                    
                    jcur_pos++;
                }
                else
                {
                    ipos+=nr;
                }
            }
            
            
            icur_pos++;
        }
        else
        {
            ipos+=nr*(ityp+1);
        }
    }
    
    
    delete [] tmp;
    
    for(int i=0;i<no_types;i++)
        delete [] type_ref[i];
    if(no_types)
        delete [] type_ref;

}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam::set_fs(int no_files
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
    if(narg)
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
    

    
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    if(narg)
        delete [] arg;
    
    
    if(line_read(fp,line)==-1)
        error->abort("eam potential file ended immaturely");
    
    
    
    int ipos=0;
    int tot=tot_no_types*nrho+tot_no_types*tot_no_types*nr
    +tot_no_types*(tot_no_types+1)/2*nr;
    TYPE0* tmp;
    CREATE1D(tmp,tot);
    while (ipos<tot)
    {
        if(line_read(fp,line)==-1)
            error->abort("eam potential file ended immaturely");
        
        narg=mapp->parse_line(line,arg);
        
        if(ipos+narg>tot)
            error->abort("eam potential file ended immaturely");
        
        for(int i=0;i<narg;i++)
        {
            tmp[ipos]=atof(arg[i]);
            ipos++;
            delete [] arg[i];
        }
        if(narg)
            delete [] arg;
        
    }
    
    delete [] line;
    if(atoms->my_p_no==0)
        fclose(fp);
    
    allocate();
    
    
    int icur_pos,jcur_pos;
    int component;
    
    
    
    ipos=0;
    icur_pos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ityp==type_ref[icur_pos][0])
        {
            component=type_ref[icur_pos][1];
            
            for(int i=0;i<nrho;i++)
                F_arr[component][i][0]=tmp[ipos+i];
            
            ipos+=nrho;
            
            
            jcur_pos=0;
            for(int jtyp=0;jtyp<tot_no_types;jtyp++)
            {
                if(jtyp==type_ref[jcur_pos][0])
                {
                    
                    component=type2rho[type_ref[icur_pos][1]][type_ref[jcur_pos][1]];
                    
                    for(int i=0;i<nr;i++)
                        rho_arr[component][i][0]=tmp[ipos+i];
                    
                    ipos+=nr;
                    jcur_pos++;
                }
                else
                {
                    ipos+=nr;
                }
            }
            
            icur_pos++;
        }
        else
        {
            ipos+=nr*tot_no_types+nrho;
        }
    }
    
    
    icur_pos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ityp==type_ref[icur_pos][0])
        {
            jcur_pos=0;
            for(int jtyp=0;jtyp<ityp+1;jtyp++)
            {
                if(jtyp==type_ref[jcur_pos][0])
                {
                    component=type2phi[type_ref[icur_pos][1]][type_ref[jcur_pos][1]];
                    
                    for(int i=0;i<nr;i++)
                        phi_r_arr[component][i][0]=tmp[ipos+i];
                    
                    ipos+=nr;
                    
                    jcur_pos++;
                }
                else
                {
                    ipos+=nr;
                }
            }
            
            
            icur_pos++;
        }
        else
        {
            ipos+=nr*(ityp+1);
        }
    }
    
    
    delete [] tmp;
    
    for(int i=0;i<no_types;i++)
        delete [] type_ref[i];
    if(no_types)
        delete [] type_ref;
}
/*--------------------------------------------
 read line and broadcast
 --------------------------------------------*/
int ForceField_eam::line_read(FILE* file
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
void ForceField_eam::clean_up()
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
void ForceField_eam::allocate()
{
    int no_types=atom_types->no_types;
    
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
                CREATE1D(F_arr[i][j],7);
        
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
                CREATE1D(F_arr[i][j],7);
    }
    
    allocated=1;
}
/*--------------------------------------------
 set the arrays
 --------------------------------------------*/
void ForceField_eam::set_arrays()
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
void ForceField_eam::interpolate(int n,TYPE0 delta
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

