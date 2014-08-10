
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
    allocated=0;
    eam_mode=NOT_SET;
    
    

    CREATE1D(nrgy_strss,7);
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam::~ForceField_eam()
{
    delete nrgy_strss;
    if(allocated) clean_up();
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam::
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
    
    int itype,jtype,icomp,jcomp;
    TYPE0 dx0,dx1,dx2,rsq,rhoip,rhojp,z2p,z2;
    TYPE0 r,p,r_inv,psip,phi,phip,fpair;
    int m;
    TYPE0* coef;
    
    int** neighbor_list=neighbor->neighbor_list;
    int* neighbor_list_size=neighbor->neighbor_list_size;
    
    nrgy_strss[0]=0.0;
    if (st_clc)
        for (int i=1;i<7;i++)
            nrgy_strss[i]=0.0;
    
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
            
            if(rsq < cut_sq)
            {
                r=sqrt(rsq);
                
                p = r*dr_inv + 1.0;
                m = static_cast<int> (p);
                m = MIN(m,nr-1);
                p -= m;
                p = MIN(p,1.0);
                
                coef=rho_arr[type2rho[itype][jtype]][m];
                rho[iatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(jatm<natms)
                {
                    coef=rho_arr[type2rho[jtype][itype]][m];
                    rho[jatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                }
            }
        }
    }
    
    atoms->update(rho_n);

    
    int tot_natms=atoms->natms+atoms->natms_ph;
    
    for(iatm=0;iatm<tot_natms;iatm++)
    {
        
        p = rho[iatm]*drho_inv + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        dF[iatm]=(coef[6]*p+coef[5])*p+coef[4];
        
        if(iatm<natms)
            nrgy_strss[0]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
    }
    
    
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
            
            if(rsq < cut_sq)
            {
                r=sqrt(rsq);
                
                p=r*dr_inv+1.0;
                m=static_cast<int>(p);
                m=MIN(m,nr-1);
                p-=m;
                p=MIN(p,1.0);
    
                
                coef=rho_arr[type2rho[itype][jtype]][m];
                rhoip=(coef[6]*p + coef[5])*p+coef[4];
                coef=rho_arr[type2rho[jtype][itype]][m];
                rhojp=(coef[6]*p + coef[5])*p+coef[4];
                coef= phi_r_arr[type2phi[itype][jtype]][m];
                z2p =(coef[6]*p + coef[5])*p+coef[4];
                z2 =((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                
                r_inv=1.0/r;
                phi=z2*r_inv;
                phip=z2p*r_inv-phi*r_inv;
                psip=dF[iatm]*rhojp+dF[jatm]*rhoip+phip;
                fpair=-psip*r_inv;
                
                f[icomp]+=dx0*fpair;
                f[icomp+1]+=dx1*fpair;
                f[icomp+2]+=dx2*fpair;
                
                

                if(jatm<natms)
                {
                    nrgy_strss[0]+=phi;
                    f[jcomp]-=dx0*fpair;
                    f[jcomp+1]-=dx1*fpair;
                    f[jcomp+2]-=dx2*fpair;
                    
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
                    nrgy_strss[0]+=phi;
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
    TYPE0 r,p;
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
                
                p=r*dr_inv+1.0;
                m=static_cast<int>(p);
                m=MIN(m,nr-1);
                p-=m;
                p=MIN(p,1.0);

                coef=rho_arr[type2phi[itype][jtype]][m];
                en+=(((coef[3]*p+coef[2])*p+coef[1])*p+coef[0])/r;
                
                coef=rho_arr[type2rho[itype][jtype]][m];
                rho[iatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                if(jatm<natms)
                {
                    coef=rho_arr[type2rho[jtype][itype]][m];
                    rho[jatm]+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
                }
                

            }
        }
        
        p=rho[iatm]*drho_inv+1.0;
        m=static_cast<int> (p);
        m=MIN(m,nr-1);
        p-=m;
        p=MIN(p,1.0);
        itype=type[iatm];
        coef=F_arr[itype][m];
        en+=((coef[3]*p+coef[2])*p+coef[1])*p+coef[0];
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
    dF_n=atoms->add<TYPE0>(1,1,"dF");
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_eam::fin()
{
    atoms->del(dF_n);
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
        set_funcfl(narg-2,&arg[2]);
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
    set_arrays();
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam::set_setfl(int no_files
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
        
        for(int i=0;i<2;i++)
            if(line_read(fp,line)==-1)
                error->abort("eam potential file ended immaturely");
        
        
        narg=mapp->parse_line(line,arg);
        if(narg!=5)
            error->abort("wrong line in eam file");
        
        nrhos[ityp]=atoi(arg[0]);
        nrs[ityp]=atoi(arg[2]);
        drhos[ityp]=atof(arg[1]);
        drs[ityp]=atof(arg[3]);
        
        
        CREATE1D(tmp_F[ityp],nrhos[ityp]);
        CREATE1D(tmp_rho[ityp],nrs[ityp]);
        CREATE1D(tmp_zi[ityp],nrs[ityp]);
        
        for(int i=0;i<5;i++)
            delete [] arg[i];
        delete [] arg;
        
        if(nrhos[ityp]%5!=0 || nrs[ityp]%5!=0)
            error->abort("nro and nrho remainder by 5 should be zero");
        
        nlines=nrhos[ityp]/5;
        no=0;
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
        no=0;
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
        no=0;
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
                tmp0=coef1*tmp_rho[ityp][k-1]+
                coef2*tmp_rho[ityp][k]+coef3*tmp_rho[ityp][k+1]+
                coef4*tmp_rho[ityp][k+2];
                
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
                tmp1=coef1*tmp_rho[jtyp][k-1]+
                coef2*tmp_rho[jtyp][k]+coef3*tmp_rho[jtyp][k+1]+
                coef4*tmp_rho[jtyp][k+2];
                
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
    
    cut_sq=(nr-1.0)*(nr-1.0)*dr*dr;
}
/*--------------------------------------------
 read setfiles
 --------------------------------------------*/
void ForceField_eam::set_funcfl(int no_files
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
            component=type2rho[0][type_ref[icur_pos][1]];
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
    cut_sq=(nr-1.0)*(nr-1.0)*dr*dr;
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
    cut_sq=(nr-1.0)*(nr-1.0)*dr*dr;
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

