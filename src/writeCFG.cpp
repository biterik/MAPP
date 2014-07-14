/*--------------------------------------------
 Created by Sina on 1/29/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include "writeCFG.h"
#include "error.h"
#include "memory.h"
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
WriteCFG::
WriteCFG(MAPP* mapp,int narg,char** arg) : InitPtrs(mapp)
{
    /*
    
    if (narg<3)
        error->abort("wrong number of arguments in the write cfg command");
    
    
    if(strcmp(arg[2],"x"))
        error->abort("first argument should be x");
    
    CREATE1D(args,narg-1);
    int length;
    for(int i=1;i<narg;i++)
    {
        length=static_cast<int>(strlen(arg[i]))+1;
        CREATE1D(args[i-1],length);
        for(int j=0;j<length;j++)
            args[i-1][j]=arg[i][j];
    }
    ncomp=narg-2;
    
    
    
    CREATE1D(vec_nos,ncomp);
    tot_dims=0;
    for (int i=0;i<ncomp;i++)
    {
        vec_nos[i]=peratom->find<double>((const char*)args[i+1]);
        tot_dims+=peratom->dim<double>(vec_nos[i]);
    }
    */
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
WriteCFG::~WriteCFG()
{
    /*
    delete [] args;
    delete [] vec_nos;
     */
}
/*--------------------------------------------
 write the cfg file
 --------------------------------------------*/
void WriteCFG::write_file(int timestep)
{
    /*
    char filename[MAXCHAR];
    sprintf (filename, "%s.%08d.cfg",args[0],timestep);
    FILE* cfgfile;
    cfgfile=fopen(filename,"w");
    

    
    int ityp_natms;
    int typ_idx=peratom->find<int>("type");
    
    int n,m,vcomp,dim;

    
    
    
    
    peratom->x2s(peratom->natms);
    
    //wirte the header
    if(comm->myno==0)
    {
        fprintf(cfgfile,"Number of particles = %d\n",peratom->natms_tot);
        fprintf(cfgfile,"A = 1 Angstrom (basic length-scale)\n");
        
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                fprintf(cfgfile,"H0(%d,%d) = %lf A\n",i+1,j+1,peratom->H[i][j]);
        
        fprintf(cfgfile,".NO_VELOCITY.\n");
    }

    
    if(comm->myno==0)
    {
        fprintf(cfgfile,"entry_count = %d\n",tot_dims);
        n=0;
        if(ncomp>1)
        {
            for(int icomp=1;icomp<ncomp;icomp++)
                for(int idim=0;idim<peratom->dim<double>(vec_nos[icomp]);idim++)
                    fprintf(cfgfile,
                    "auxiliary[%d] = %s%d [reduced unit]\n"
                    ,n++,args[1+icomp],idim+1);
        }
    }
    
    
 
    double* buff;
    int size,tot_size;
    int* buff_size;
    int* tot_buff_size;
    CREATE1D(buff_size,comm->totp);
    CREATE1D(tot_buff_size,comm->totp);
    MPI_Status status;
    
    
    for(int ityp=0;ityp<peratom->ntyps;ityp++)
    {
        
        ityp_natms=0;
        for(int iatm=0;iatm<peratom->natms;iatm++)
            if(peratom->int_vec[typ_idx][iatm]==ityp)
                ityp_natms++;
        
        size=ityp_natms*tot_dims;
        CREATE1D(buff,size);
        
        
        m=0;
        for(int iatm=0;iatm<peratom->natms;iatm++)
        {
            if(peratom->int_vec[typ_idx][iatm]==ityp)
            {
                
                for(int ivec=0;ivec<ncomp;ivec++)
                {
                    vcomp=vec_nos[ivec];
                    dim=peratom->dim<double>(vcomp);
                    
                    for(int idim=0;idim<dim;idim++)
                        buff[m++]=peratom->double_vec[vcomp][iatm*dim+idim];
                }
            }
        }
        
        MPI_Allreduce(&size,&tot_size,
        1,MPI_INT,MPI_SUM,world);
        
        if(comm->myno==0)
            GROW(buff,size,tot_size);

        
        for(int i=0;i<comm->totp;i++)
        {
            buff_size[i]=0;
            tot_buff_size[i]=0;
        }
        
        buff_size[comm->myno]=size;
        
        MPI_Allreduce(&(buff_size[0]),&(tot_buff_size[0]),
        comm->totp,MPI_INT,MPI_MAX,world);
        
        
        
        
        if(comm->myno!=0)
            MPI_Send(&buff[0],size,MPI_DOUBLE,0,comm->myno,world);
        
        for (int iprc=1;iprc<comm->totp;iprc++)
        {
            if(comm->myno==0)
            {
                MPI_Recv(&(buff[size]),tot_buff_size[iprc],MPI_DOUBLE,iprc,iprc,world,&status);
                size+=tot_buff_size[iprc];
            }
        }
        
        
        
        if(comm->myno==0)
        {
            fprintf(cfgfile,"%lf\n",peratom->m[ityp]);
            fprintf(cfgfile,"%s\n",peratom->atom_names[ityp]);
            
            m=tot_size/tot_dims;
            n=0;
            for(int i=0;i<m;i++)
            {
                for(int j=0;j<tot_dims;j++)
                    fprintf(cfgfile,"%lf ",buff[n++]);
                fprintf(cfgfile,"\n");
                
            }
        }
        
        if(size)
            delete [] buff;
        
        
    }
    
    
    delete [] buff_size;
    delete [] tot_buff_size;
    fclose(cfgfile);
    peratom->s2x(peratom->natms);*/
}

