
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
    
    

    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam::~ForceField_eam()
{
    if(allocated) clean_up();
}
/*--------------------------------------------
 force calculation
 --------------------------------------------*/
void ForceField_eam::
force_calc(int st_clc,TYPE0* en_st)
{
    
}
/*--------------------------------------------
 energy calculation
 --------------------------------------------*/
TYPE0 ForceField_eam::energy_calc()
{
    return 0.0;
}
/*--------------------------------------------
 init before running
 --------------------------------------------*/
void ForceField_eam::init()
{
    
}
/*--------------------------------------------
 fin after running
 --------------------------------------------*/
void ForceField_eam::fin()
{
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
void ForceField_eam::coef(int narg,char** arg)
{

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
            
            narg=mapp->parse_line(line,arg);
            if(narg!=5)
                error->abort("wrong line in eam file");
            
            for(int j=0;j<5;j++)
                tmp_F[ityp][no++]=atof(arg[j]);
            
            for(int j=0;j<5;j++)
                delete [] arg[j];
            delete [] arg;
            
        }
        
        nlines=nrs[ityp]/5;
        no=0;
        for(int i=0;i<nlines;i++)
        {
            err=line_read(fp,line);
            if(err==-1)
                error->abort("eam potential file ended immaturely");
            
            narg=mapp->parse_line(line,arg);
            if(narg!=5)
                error->abort("wrong line in eam file");
            
            for(int j=0;j<5;j++)
                tmp_zi[ityp][no++]=atof(arg[j]);
            
            for(int j=0;j<5;j++)
                delete [] arg[j];
            delete [] arg;
            
        }
        
        nlines=nrs[ityp]/5;
        no=0;
        for(int i=0;i<nlines;i++)
        {
            err=line_read(fp,line);
            if(err==-1)
                error->abort("eam potential file ended immaturely");
            
            narg=mapp->parse_line(line,arg);
            if(narg!=5)
                error->abort("wrong line in eam file");
            
            for(int j=0;j<5;j++)
                tmp_rho[ityp][no++]=atof(arg[j]);
            
            for(int j=0;j<5;j++)
                delete [] arg[j];
            delete [] arg;
            
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
    
    if(atoms->my_p_no==0)
        fclose(fp);
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
              &&ityp<narg-1)
            ityp++;
        
        if(ityp==narg-1)
            error->abort("atom type not found in eam file");
        type_ref[i][0]=ityp;
        type_ref[i][1]=i;
    }
    
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
    
    
        
    
    for(int j=0;j<narg;j++)
        delete [] arg[j];
    delete [] arg;
    
    if(line_read(fp,line)==-1)
        error->abort("eam potential file ended immaturely");
    narg=mapp->parse_line(line,arg);
    
    
    if(narg!=5)
        error->abort("wrong format in eam potential file");
    
    nrho=atoi(arg[0]);
    nr=atoi(arg[2]);
    drho=atof(arg[1]);
    dr=atof(arg[3]);
    
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
    
    no_lines=nrho/5;
    icur_pos=0;
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        if(ityp==type_ref[icur_pos][0])
        {
            no=0;
            component=type_ref[icur_pos][1];
            for(int i=0;i<no_lines;i++)
            {
                err=line_read(fp,line);
                if(err==-1)
                    error->abort("eam potential file ended immaturely");
                
                narg=mapp->parse_line(line,arg);
                if(narg!=5)
                    error->abort("wrong line in eam file");
                
                for(int j=0;j<5;j++)
                    F_arr[component][no++][0]=atof(arg[j]);
                
                for(int j=0;j<narg;j++)
                    delete [] arg[j];
                delete [] arg;
            }
            
            icur_pos++;
        }
        else
        {
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
        for(int jtyp=0;jtyp<tot_no_types;jtyp++)
        {
            if(ityp==type_ref[icur_pos][0] && jtyp==type_ref[jcur_pos][0])
            {
                no=0;
                component=type_ref[icur_pos][1]*no_types+type_ref[jcur_pos][1];
                for(int i=0;i<no_lines;i++)
                {
                    err=line_read(fp,line);
                    if(err==-1)
                        error->abort("eam potential file ended immaturely");
                    
                    narg=mapp->parse_line(line,arg);
                    if(narg!=5)
                        error->abort("wrong line in eam file");
                    
                    for(int j=0;j<5;j++)
                        rho_arr[component][no++][0]=atof(arg[j]);
                    
                    for(int j=0;j<narg;j++)
                        delete [] arg[j];
                    delete [] arg;
                }
                
                icur_pos++;
                jcur_pos++;
            }
            else
            {
                for(int i=0;i<no_lines;i++)
                    if(line_read(fp,line)==-1)
                        error->abort("eam potential file ended immaturely");
                
                if(ityp==type_ref[icur_pos][0])
                    icur_pos++;
                if(jtyp==type_ref[jcur_pos][0])
                    jcur_pos++;
            }
        }
        
    }
    
    
    
    no_lines=nr/5;
    icur_pos=0;
    jcur_pos=0;
    
    for(int ityp=0;ityp<tot_no_types;ityp++)
    {
        for(int jtyp=0;jtyp<ityp+1;jtyp++)
        {
            if(ityp==type_ref[icur_pos][0] && jtyp==type_ref[jcur_pos][0])
            {
                no=0;
                component=COMP(type_ref[icur_pos][1],type_ref[jcur_pos][1]);
                for(int i=0;i<no_lines;i++)
                {
                    err=line_read(fp,line);
                    if(err==-1)
                        error->abort("eam potential file ended immaturely");
                    
                    narg=mapp->parse_line(line,arg);
                    if(narg!=5)
                        error->abort("wrong line in eam file");
                    
                    for(int j=0;j<5;j++)
                        phi_r_arr[component][no++][0]=atof(arg[j]);
                    
                    for(int j=0;j<narg;j++)
                        delete [] arg[j];
                    delete [] arg;
                }
                
                icur_pos++;
                jcur_pos++;
            }
            else
            {
                for(int i=0;i<no_lines;i++)
                    if(line_read(fp,line)==-1)
                        error->abort("eam potential file ended immaturely");
                
                if(ityp==type_ref[icur_pos][0])
                    icur_pos++;
                if(jtyp==type_ref[jcur_pos][0])
                    jcur_pos++;
            }
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
                delete [] rho_arr[i*(no_types+1)][j];
            delete [] rho_arr[i*(no_types+1)];
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
            CREATE1D(rho_arr[i*(no_types+1)],nr+1);
        for(int i=0;i<no_types;i++)
            for(int j=0;j<nr+1;j++)
                CREATE1D(rho_arr[i*(no_types+1)][j],7);
        
        for(int i=0;i<no_types;i++)
            for(int j=0;j<no_types;j++)
                if(i!=j)
                    rho_arr[i*no_types+j]
                    =rho_arr[i*(no_types+1)];
        
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






