/*--------------------------------------------
 Created by Sina on 07/23/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include <cmath>
#include <stdlib.h>
#include "read_cfg.h"
#include "error.h"
#include "memory.h"
#include "atom_types.h"
using namespace MAPP_NS;
enum{TYPE_XD_ID,TYPE_ID_XD,ID_TYPE_XD,ID_XD_TYPE,XD_ID_TYPE,XD_TYPE_ID};
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ReadCFG::ReadCFG(MAPP* mapp,int narg,char** args)
:Read(mapp)
{
    curr_id=-1;
    
    if(narg!=3)
        error->abort("wrong command");
    
    if(atoms->dimension!=3)
        error->abort("in order to read a cfg"
                     " file the box dimension should be 3");
    
    // setup the defaults
    basic_length=1.0;
    R=1.0;
    entry_count=6;
    ext_cfg=0;
    header_cmplt=0;
    atom_cmplt=0;
    vel_chk=1;

    CREATE2D(H_x,3,3);
    CREATE2D(H_x_d,3,3);
    CREATE2D(H0,3,3);
    CREATE2D(eta,3,3);
    CREATE2D(eta_sq,3,3);
    CREATE2D(trns,3,3);
    
    M3ZERO(H_x);
    M3ZERO(H_x_d);
    M3ZERO(H0);
    M3ZERO(eta);
    M3ZERO(trns);
    
    trns[0][0]=trns[1][1]=trns[2][2]=1.0;
    H0[0][0]=H0[1][1]=H0[2][2]=1.0;
    // end of defaults
    
    
    CREATE1D(line,MAXCHAR);
    if(atoms->my_p_no==0)
    {
        cfgfile=fopen(args[2],"r");
        if(cfgfile==NULL)
            error->abort("file %s not found",args[2]);
    }
    
    
    fpos_t pos;
    while (!header_cmplt)
    {
        if(atoms->my_p_no==0)
        {
            fgetpos(cfgfile,&pos);
            fgets(line,MAXCHAR,cfgfile);
        }
        MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
        read_header();
    }
    if(atoms->my_p_no==0)
        fsetpos(cfgfile,&pos);
    
    set_box();
    
    atoms->auto_grid_proc();
 
    if(mapp->mode==MD)
    {
        x_n=atoms->find("x");
        x_d_n=atoms->find_exist("x_d");
        type_n=atoms->find("type");
        id_n=atoms->find("id");
        
        if(x_d_n>0)
        {
            vec_list=new VecLst(mapp,4,x_n,x_d_n,type_n,id_n);
            if(x_d_n<type_n && type_n<id_n)
            {
                ch_x_d=atoms->vectors[x_n].byte_size;
                ch_type=ch_x_d+atoms->vectors[x_d_n].byte_size;
                ch_id=ch_type+atoms->vectors[type_n].byte_size;
            }
            else if(x_d_n<id_n && id_n<type_n)
            {
                ch_x_d=atoms->vectors[x_n].byte_size;
                ch_id=ch_x_d+atoms->vectors[x_d_n].byte_size;
                ch_type=ch_id+atoms->vectors[id_n].byte_size;
            }
            else if(id_n<x_d_n && x_d_n<type_n)
            {
                ch_id=atoms->vectors[x_n].byte_size;
                ch_x_d=ch_id+atoms->vectors[id_n].byte_size;
                ch_type=ch_x_d+atoms->vectors[x_d_n].byte_size;
            }
            else if(id_n<type_n && type_n<x_d_n)
            {
                ch_id=atoms->vectors[x_n].byte_size;
                ch_type=ch_id+atoms->vectors[id_n].byte_size;
                ch_x_d=ch_type+atoms->vectors[type_n].byte_size;
            }
            else if(type_n<x_d_n && x_d_n<id_n)
            {
                ch_type=atoms->vectors[x_n].byte_size;
                ch_x_d=ch_type+atoms->vectors[type_n].byte_size;
                ch_id=ch_x_d+atoms->vectors[x_d_n].byte_size;
            }
            else if(type_n<id_n && id_n<x_d_n)
            {
                ch_type=atoms->vectors[x_n].byte_size;
                ch_id=ch_type+atoms->vectors[type_n].byte_size;
                ch_type=ch_id+atoms->vectors[id_n].byte_size;
            }
        }
        else
        {
            vec_list=new VecLst(mapp,3,x_n,type_n,id_n);
            if(type_n<id_n)
            {
                ch_type=atoms->vectors[x_n].byte_size;
                ch_id=ch_type+atoms->vectors[type_n].byte_size;
            }
            else
            {
                ch_id=atoms->vectors[x_n].byte_size;
                ch_type=ch_id+atoms->vectors[id_n].byte_size;
            }
        }

    }
    else if(mapp->mode==DMD)
    {
        if((entry_count-3)%2!=0)
            error->abort("wrong number of entry for DMD in cfg file");
        dmd_no_types=(entry_count-3)/2;
        if(atoms->vectors[0].dim!=dmd_no_types+3)
            atoms->vectors[0].change_dimension(3+dmd_no_types);
        x_n=0;
        c_n=atoms->find_exist("c");
        id_n=atoms->find("id");
        if(c_n==-1)
        {
            c_n=atoms->add<TYPE0>(1,dmd_no_types,"c");
        }
        else
        {
            if(atoms->vectors[c_n].dim!=dmd_no_types)
                atoms->vectors[c_n].change_dimension(dmd_no_types);
        }
        vec_list=new VecLst(mapp,3,x_n,c_n,id_n);
        if(c_n<id_n)
        {
            ch_c=atoms->vectors[x_n].byte_size;
            ch_id=ch_c+atoms->vectors[c_n].byte_size;
        }
        else
        {
            ch_id=atoms->vectors[x_n].byte_size;
            ch_c=ch_id+atoms->vectors[id_n].byte_size;
        }
        
        if(ext_cfg==0)
            error->abort("DMD mode can only read extended cfg");
    }
    
    CREATE1D(ch_buff,vec_list->byte_size);
    
    while (!atom_cmplt)
    {
        if(atoms->my_p_no==0)
            fgets(line,MAXCHAR,cfgfile);
        MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
        if(mapp->mode==MD)
            read_atom();
        else if(mapp->mode==DMD)
            read_atom_dmd();
    }
    int loc_no=atoms->natms;
    int tot_no=0;
    MPI_Allreduce(&loc_no,&tot_no,1,MPI_INT,MPI_SUM,world);
    
    if (tot_no!=atoms->tot_natms)
        error->abort("the number of atoms dont match"
                     " recheck your cfg file: %d %d",
                     atoms->tot_natms,tot_no);
    
    delete vec_list;
    delete [] line;
    delete [] trns;
    
    for(int i=0;i<3;i++)
    {
        delete [] eta_sq[i];
        delete [] eta[i];
        delete [] H0[i];
        delete [] H_x[i];
        delete [] H_x_d[i];
    }
    
    delete [] eta_sq;
    delete [] eta;
    delete [] H0;
    delete [] H_x;
    delete [] H_x_d;
    delete [] ch_buff;
    
    if(atoms->my_p_no==0)
        fclose(cfgfile);

    
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
ReadCFG::~ReadCFG()
{
    
}
/*--------------------------------------------
 reads the header of the cfg file
 --------------------------------------------*/
void ReadCFG::read_header()
{
    char* command;
    int narg = mapp->hash_remover(line,command);
    TYPE0 tmp;
    int icmp,jcmp,tmpno;
    char* strtmp1;
    char* strtmp2;
    char* strtmp3;
    
    CREATE1D(strtmp1,MAXCHAR);
    CREATE1D(strtmp2,MAXCHAR);
    CREATE1D(strtmp3,MAXCHAR);
    
    if (narg==0)
    {
        delete [] strtmp3;
        delete [] strtmp2;
        delete [] strtmp1;
        return;
    }
    else if(narg==1)
    {
        if(!strcmp(command,".NO_VELOCITY."))
        {
            vel_chk=0;
            ext_cfg=1;
            entry_count-=3;
        }
        else if(ext_cfg)
            header_cmplt=1;
        
        else
            error->abort("unknown command: %s",command);
    }
    else if(narg==3)
    {
        if(sscanf(command,"Transform(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if (icmp>2 || icmp<0
                || jcmp>2 || jcmp<0)
                error->abort("unknown component: %s",command);
            trns[icmp][jcmp]=tmp;
        }
        else if(sscanf(command,"eta(%d,%d) = %lf",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if (icmp>2||icmp<0||jcmp>2||jcmp<0)
                error->abort("unknown component %s",command);
            eta[icmp][jcmp]=tmp;
        }
        else if(sscanf(command,
                       "entry_count = %d",&tmpno)==1)
        {
            entry_count=tmpno;
            ext_cfg=1;
            int mincomp=3+(3*vel_chk);
            if (entry_count < mincomp)
                error->abort("wrong number of entry_count");
        }
        else error->abort("unknown command: %s",command);
    }
    else if(narg==4)
    {
        if(sscanf(command,
                  "H0(%d,%d) = %lf A",&icmp,&jcmp,&tmp)==3)
        {
            icmp--;
            jcmp--;
            if (icmp>2||icmp<0||jcmp>2||jcmp<0)
                error->abort("unknown component: %s",command);
            else
                H0[icmp][jcmp]=tmp;
        }
        else if(sscanf(command,"R = %lf %s",&tmp,strtmp1)==2)
        {
            R=tmp;
        }
        else if(sscanf(command,"auxiliary[%d] = %s %s",&icmp,strtmp1,strtmp2)==3)
        {
            int mincomp=3+(3*vel_chk);
            if(icmp+mincomp+1>entry_count)
                error->abort("auxilary component larger than entry_count");
        }
        else
            error->abort("unknown command: %s",command);
    }
    else if(narg==5)
    {
        if(sscanf(command,
                  "Number of particles = %d",&tmpno)==1)
            atoms->tot_natms=tmpno;
        else if(sscanf(command,"auxiliary[%d] = %s %s %s",&icmp,strtmp1,strtmp2,strtmp3)==4)
        {
            int mincomp=3+(3*vel_chk);
            if(icmp+mincomp+1>entry_count)
                error->abort("auxilary component larger than entry_count");
        }

        else
            error->abort("unknown command: %s",command);
    }
    else if(narg==6)
    {
        if(sscanf(command,
                  "A = %lf Angstrom (basic length-scale)",&tmp)==1)
            basic_length=tmp;
        else
            error->abort("unknown command: %s",command);
    }
    else if (narg==8&&ext_cfg==0)
        header_cmplt=1;
    else
        error->abort("unknown command: %s",command);
    
    delete [] command;
    delete [] strtmp3;
    delete [] strtmp2;
    delete [] strtmp1;
}
/*--------------------------------------------
 calculates H from H0, Transform, and eta;
 ** make sure H is zeroed before;
 for now we disregard eta;
 remeber to fix it later;
 --------------------------------------------*/
void ReadCFG::set_box()
{
    TYPE0** Ht;
    TYPE0* sq;
    TYPE0* b;
    CREATE2D(Ht,3,3);
    CREATE1D(sq,3);
    CREATE1D(b,3);
    TYPE0 babs;
    
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            for (int k=0;k<3;k++)
                H_x[i][j]+=H0[i][k]*trns[k][j];
    
    int chk=1;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            if(eta[i][j]!=0.0)
                chk=0;
    if(chk==0)
    {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                eta[i][j]*=2.0;
        for(int i=0;i<3;i++)
            eta[i][i]++;
        
        M3sqroot(eta,eta_sq);
        
        M3EQV(H_x,H0);
        M3ZERO(H_x);
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                for (int k=0;k<3;k++)
                    H_x[i][j]+=H0[i][k]*eta_sq[k][j];
    }
    
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
        {
            H_x[i][j]*=basic_length;
            atoms->H[i][j]=H_x[i][j];
        }
    M3EQV(H_x,H_x_d);
    
    if (M3DET(H_x)==0)
        error->abort("H determinant is zero");
    TYPE0 det;
    
    for (int i=0;i<3;i++)
    {
        sq[i]=0.0;
        for (int j=0;j<3;j++)
        {
            sq[i]+=H_x[i][j]*H_x[i][j];
            Ht[i][j]=0.0;
        }
    }
    
    Ht[0][0]=sqrt(sq[0]);
    
    for (int i=0;i<3;i++)
        Ht[1][0]+=H_x[0][i]*H_x[1][i];
    Ht[1][0]/=Ht[0][0];
    
    Ht[1][1]=sqrt(sq[1]-Ht[1][0]*Ht[1][0]);
    
    for (int i=0;i<3;i++)
        Ht[2][0]+=H_x[0][i]*H_x[2][i];
    
    Ht[2][0]/=Ht[0][0];
    
    
    b[0]=H_x[0][1]*H_x[1][2]-H_x[0][2]*H_x[1][1];
    b[1]=H_x[0][2]*H_x[1][0]-H_x[0][0]*H_x[1][2];
    b[2]=H_x[0][0]*H_x[1][1]-H_x[0][1]*H_x[1][0];
    babs=sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
    b[0]/=babs;
    b[1]/=babs;
    b[2]/=babs;
    
    for (int i=0;i<3;i++)
        Ht[2][2]+=H_x[2][i]*b[i];
    Ht[2][1]=sqrt(sq[2]-Ht[2][2]*Ht[2][2]-Ht[2][0]*Ht[2][0]);
    
    M3EQV(Ht,atoms->H);
    M3EQV(Ht,H_x);
    M3INV(atoms->H,atoms->B,det);
    
    for(int i=0;i<3;i++)
        delete Ht[i];
    delete [] Ht;
    delete [] b;
    delete [] sq;
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void ReadCFG::read_atom()
{
    char** arg;
    TYPE0 mass;
    TYPE0* buff;
    int narg=mapp->parse_line(line,arg);
    
    if(atoms->my_p_no==0)
    {
        if(feof(cfgfile))
            atom_cmplt=1;
    }
    MPI_Bcast(&atom_cmplt,1,MPI_INT,0,world);
    if(atom_cmplt)
    {
        if(narg)
            delete [] arg;
        return;
    }
    
    
    if (narg!=8 && ext_cfg==0)
        error->abort("wrong format: %s",line);
    if (ext_cfg && !(narg==1 || narg==entry_count))
        error->abort("wrong extended format: %s",line);
    
    if(ext_cfg)
    {
        if(narg==1)
        {
            mass=static_cast<TYPE0>(atof(arg[0]));
            if(atoms->my_p_no==0)
                fgets(line,MAXCHAR,cfgfile);
            MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
            
            narg=mapp->parse_line(line,arg);
            if(narg!=1)
                error->abort("error cfg: %s",line);
            if(mass<=0.0)
                error->abort("mass cannot be zero");
            last_type=
            atom_types->add_type(mass,arg[0]);
        }
        else if(narg==entry_count)
        {
            CREATE1D(buff,6);
            for (int i=0;i<6;i++)
                buff[i]=0.0;
            if (vel_chk)
            {
                for (int i=0;i<3;i++)
                    buff[i]=static_cast<TYPE0>(atof(arg[i]));
                for (int i=3;i<6;i++)
                    buff[i]=static_cast<TYPE0>(atof(arg[i]))*R;
                add_atom_read_x(last_type,buff);
            }
            else
            {
                for (int i=0;i<3;i++)
                    buff[i]=static_cast<TYPE0>(atof(arg[i]));
                add_atom_read_x(last_type,buff);
            }
            
            delete [] buff;
        }
        else
            error->abort("unknown line: %s",line);
    }
    else
    {
        if(narg==8)
        {
            mass=static_cast<TYPE0>(atof(arg[0]));
            last_type=
            atom_types->add_type(mass,arg[1]);
            CREATE1D(buff,6);
            for (int i=0;i<6;i++) buff[i]=0.0;
            
            for (int i=0;i<3;i++)
                buff[i]=static_cast<TYPE0>(atof(arg[i+2]));
            for (int i=3;i<6;i++)
                buff[i]=static_cast<TYPE0>(atof(arg[i+2]))*R;
            add_atom_read_x(last_type,buff);
            delete [] buff;
        }
        else
            error->abort("unknown line: %s",line);
    }
    
    
    for(int i=0;i<narg;i++)
        delete [] arg[i];
    if(narg)
        delete [] arg;
}
/*--------------------------------------------
 reads the atom section of the cfg file
 --------------------------------------------*/
void ReadCFG::read_atom_dmd()
{
    char** arg;
    TYPE0 mass;
    TYPE0* buff;
    int narg=mapp->parse_line(line,arg);
    
    if(atoms->my_p_no==0)
    {
        if(feof(cfgfile))
            atom_cmplt=1;
    }
    MPI_Bcast(&atom_cmplt,1,MPI_INT,0,world);
    if(atom_cmplt)
    {
        if(narg)
            delete [] arg;
        return;
    }
    
    
    if (narg!=8 && ext_cfg==0)
        error->abort("wrong format: %s",line);
    if (ext_cfg && !(narg==1 || narg==entry_count))
        error->abort("wrong extended format: %s",line);
    

    if(narg==1)
    {
        mass=static_cast<TYPE0>(atof(arg[0]));
        if(atoms->my_p_no==0)
            fgets(line,MAXCHAR,cfgfile);
        MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
        
        narg=mapp->parse_line(line,arg);
        if(narg!=1)
            error->abort("error cfg: %s",line);
        if(mass<=0.0)
            error->abort("mass cannot be zero");
        last_type=
        atom_types->add_type(mass,arg[0]);
        if(last_type>=dmd_no_types)
            error->abort("exceeded the number of types in DMD");
    }
    else if(narg==entry_count)
    {
        CREATE1D(buff,entry_count);
        for (int i=0;i<entry_count;i++)
            buff[i]=static_cast<TYPE0>(atof(arg[i]));
        
        add_atom_read_x(last_type,buff);
        
        delete [] buff;
    }
    else
        error->abort("unknown line: %s",line);

    
    
    for(int i=0;i<narg;i++)
        delete [] arg[i];
    if(narg)
        delete [] arg;
}
/*--------------------------------------------
 addatom_read_x
 --------------------------------------------*/
void ReadCFG::add_atom_read_x(int t,TYPE0* buff)
{
    curr_id++;
    
    if(mapp->mode==MD)
    {
        TYPE0* x_d;
        CREATE1D(x_d,3);
        
        for(int i=0;i<3;i++)
        {
            while(buff[i]>=1)
                buff[i]--;
            while(buff[i]<0)
                buff[i]++;
        }
        
        V3ZERO(x_d);
        
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                x_d[j]+=buff[i+3]*H_x_d[i][j]*R;
        for (int i=0;i<3;i++)
            buff[i+3]=x_d[i];
        
        delete [] x_d;
        
        for(int i=0;i<3;i++)
            if(!(atoms->s_lo[i]<=buff[i]
                 && buff[i]<atoms->s_hi[i]))
                return;
        
        memcpy(ch_buff,buff,3*sizeof(TYPE0));
        memcpy(&ch_buff[ch_type],&t,sizeof(int));
        memcpy(&ch_buff[ch_id],&curr_id,sizeof(int));
        
        if(x_d_n>-1)
            memcpy(&ch_buff[ch_x_d],&buff[3],3*sizeof(TYPE0));
        
    }
    else
    {

        
        for(int i=0;i<3;i++)
        {
            while(buff[i]>=1)
                buff[i]--;
            while(buff[i]<0)
                buff[i]++;
        }
        
        for(int i=0;i<3;i++)
            if(!(atoms->s_lo[i]<=buff[i]
                 && buff[i]<atoms->s_hi[i]))
                return;
        
        memcpy(ch_buff,buff,(3+dmd_no_types)*sizeof(TYPE0));
        memcpy(&ch_buff[ch_id],&curr_id,sizeof(int));
        memcpy(&ch_buff[ch_c],&buff[3+dmd_no_types],dmd_no_types*sizeof(TYPE0));

    }
    
    atoms->unpack(ch_buff,0,1,vec_list);
    
}
/*--------------------------------------------
 calculates square root of 3x3 matrix
 ref: L. P. Franca
 An Algorithm to Compute The Square Root of
 a 3x3 Positive Definite Matrix
 Computers Math. Applic. Vol. 18, No. 5,
 pp. 459-466, 1989
 --------------------------------------------*/
void ReadCFG::M3sqroot(TYPE0** A,TYPE0** Asq)
{
    TYPE0 IA=0;
    for(int i=0;i<3;i++)
        IA+=A[i][i];
    TYPE0 IIA=0;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            IIA-=A[i][j]*A[j][i];
    
    IIA+=IA*IA;
    IIA*=0.5;
    TYPE0 IIIA=M3DET(A);
    TYPE0 k=IA*IA-3*IIA;
    if(k<0.0)
        error->abort("error in taking the square root of"
                     " the matrix: it is not positive definite");
    if(k<TOLERANCE)
    {
        if(IA<=0.0)
            error->abort("error in taking the square root of"
                         " the matrix: it is not positive definite");
        M3ZERO(Asq);
        for(int i=0;i<3;i++)
            Asq[i][i]=sqrt(IA/3.0);
        return;
    }
    
    TYPE0 l=IA*(IA*IA -4.5*IIA)+13.5*IIIA;
    
    TYPE0 temp=l/(k*sqrt(k));
    if(temp>1.0||temp<-1.0)
        error->abort("error in taking the square root of"
                     " the matrix: it is not positive definite");
    TYPE0 phi=acos(temp);
    TYPE0 lambda=sqrt((1.0/3.0)*(IA+2*sqrt(k)*cos(phi/3.0)));
    
    TYPE0 IIIAsq=sqrt(IIIA);
    TYPE0 y=-lambda*lambda+IA+2*IIIAsq/lambda;
    if(y<0.0)
        error->abort("error in taking the square root of"
                     " the matrix: it is not positive definite");
    TYPE0 IAsq=lambda+sqrt(y);
    TYPE0 IIAsq=0.5*(IAsq*IAsq-IA);
    
    TYPE0 coef0=IAsq*IIAsq-IIIAsq;
    if(coef0==0)
        error->abort("error in taking the square root of"
                     " the matrix: it is not positive definite");
    coef0=1.0/coef0;
    
    M3ZERO(Asq);
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                Asq[i][j]-=coef0*A[i][k]*A[k][j];
    
    TYPE0 coef1=coef0*(IAsq*IAsq-IIAsq);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            Asq[i][j]+=coef1*A[i][j];
    
    TYPE0 coef2=coef0*IAsq*IIIAsq;
    for(int i=0;i<3;i++)
        Asq[i][i]+=coef2;
    
}
