#include "command_fix.h"
#include "atom_types.h"
#include <stdlib.h>
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
CommandFix::CommandFix(MAPP* mapp
,int narg,char** args):InitPtrs(mapp)
{
    if(mapp->mode==MD)
        md(narg,args);
    else if(mapp->mode==DMD)
        dmd(narg,args);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
CommandFix::~CommandFix()
{
    
}
/*--------------------------------------------
 md
 --------------------------------------------*/
void CommandFix::md(int narg,char** args)
{
    int dof_n,iarg;
    int x_dim=atoms->vectors[0].dim;
    char* dof_ch;
    CREATE1D(dof_ch,x_dim);
    memset(dof_ch,0,x_dim);
    
    int id_n;

    if(strcmp(args[1],"release")==0)
    {
        if(narg!=2)
            error->abort("incorrect fix command");
        
        dof_n=atoms->find_exist("dof");
        if(dof_n>-1)
            atoms->del(dof_n);
        return;
    }
    
    
    int files_started=0;
    iarg=1;
    while(iarg<narg && files_started==0)
    {
        if(strcmp(args[iarg],"x")==0)
            dof_ch[0]=1;
        else if(strcmp(args[iarg],"y")==0)
            dof_ch[1]=1;
        else if(strcmp(args[iarg],"z")==0)
            dof_ch[2]=1;
        else
            files_started=1;
        iarg++;
    }
    iarg--;
    
    if(iarg>=narg)
        error->abort("incorrect fix command");
    
    
    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    int no_atoms;
    
    int* list=NULL;
    int list_size=0;
    
    while (iarg<narg)
    {
        if(atoms->my_p_no==0)
        {
            fp=fopen(args[iarg],"r");
            if(fp==NULL)
                error->abort("file %s not found",args[iarg]);
            
            fgets(line,MAXCHAR,fp);
            sscanf(line,"%d",&no_atoms);
        }
        
        MPI_Bcast(&no_atoms,1,MPI_INT,0,world);
        GROW(list,list_size,list_size+no_atoms);
        
        if(atoms->my_p_no==0)
        {
            for(int iatm=list_size;iatm<list_size+no_atoms;iatm++)
            {
                fgets(line,MAXCHAR,fp);
                sscanf(line,"%d",&list[iatm]);
            }
            
            fclose(fp);
        }
        
        MPI_Bcast(&list[list_size],no_atoms,MPI_INT,0,world);
        
        list_size+=no_atoms;
        
        iarg++;
    }
    delete [] line;
    
    char* dof;
    dof_n=atoms->find_exist("dof");
    
    int natms=atoms->natms;
    if(dof_n==-1)
    {
        dof_n=atoms->add<char>(0,x_dim,"dof");
        atoms->vectors[dof_n].ret(dof);
        memset(dof,0,x_dim*natms);
    }
    
    atoms->vectors[dof_n].ret(dof);
    
    int* id;
    id_n=atoms->find("id");
    atoms->vectors[id_n].ret(id);
    
    for(int i=0;i<natms;i++)
        for(int j=0;j<list_size;j++)
            if(id[i]==list[j])
            {
                memcpy(&dof[i*x_dim],dof_ch,x_dim);
            }
    
    if(list_size)
        delete [] list;
    
    delete [] dof_ch;
}
/*--------------------------------------------
 dmd
 --------------------------------------------*/
void CommandFix::dmd(int narg,char** args)
{
    int dof_n,cdof_n,iarg;
    int x_dim=atoms->vectors[0].dim;
    char* dof_ch;
    CREATE1D(dof_ch,x_dim);
    memset(dof_ch,0,x_dim);
    int no_types=atom_types->no_types;
    char* cdof_ch;
    CREATE1D(cdof_ch,no_types);
    memset(cdof_ch,0,no_types);
    
    
    int id_n;
    if(strcmp(args[1],"release")==0)
    {
        if(narg!=2)
            error->abort("incorrect fix command");
        
        dof_n=atoms->find_exist("dof");
        if(dof_n>-1)
            atoms->del(dof_n);
        
        cdof_n=atoms->find_exist("cdof");
        if(cdof_n>-1)
            atoms->del(cdof_n);
        return;
    }
    
    
    
    iarg=1;
    int files_started=0,c_comp,alpha_comp;
    
    
    while(iarg<narg && files_started==0)
    {
        if(strcmp(args[iarg],"x")==0)
            dof_ch[0]=1;
        else if(strcmp(args[iarg],"y")==0)
            dof_ch[1]=1;
        else if(strcmp(args[iarg],"z")==0)
            dof_ch[2]=1;
        else if(sscanf(args[iarg],"c[%d]",&c_comp)==1)
        {
            if(c_comp>=no_types || c_comp<0)
                error->abort("invalid c component %d",c_comp);
            cdof_ch[c_comp]=1;
        }
        else if(sscanf(args[iarg],"alpha[%d]",&alpha_comp)==1)
        {
            if(alpha_comp>=no_types || alpha_comp<0)
                error->abort("invalid alpha component %d",alpha_comp);
            dof_ch[3+alpha_comp]=1;
        }
        else
            files_started=1;
        iarg++;
    }
    iarg--;
    
    if(iarg>=narg)
        error->abort("incorrect fix command");

    
    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    int no_atoms;
    
    int* list=NULL;
    int list_size=0;
    
    while (iarg<narg)
    {
        if(atoms->my_p_no==0)
        {
            fp=fopen(args[iarg],"r");
            if(fp==NULL)
                error->abort("file %s not found",args[iarg]);
            
            fgets(line,MAXCHAR,fp);
            sscanf(line,"%d",&no_atoms);
        }
        
        MPI_Bcast(&no_atoms,1,MPI_INT,0,world);
        GROW(list,list_size,list_size+no_atoms);
        
        if(atoms->my_p_no==0)
        {
            for(int iatm=list_size;iatm<list_size+no_atoms;iatm++)
            {
                fgets(line,MAXCHAR,fp);
                sscanf(line,"%d",&list[iatm]);
            }
            
            fclose(fp);
        }
        
        MPI_Bcast(&list[list_size],no_atoms,MPI_INT,0,world);
        
        list_size+=no_atoms;
        
        iarg++;
    }
    

    delete [] line;
    
    char* dof;
    dof_n=atoms->find_exist("dof");
    
    int natms=atoms->natms;
    if(dof_n==-1)
    {
        dof_n=atoms->add<char>(0,x_dim,"dof");
        atoms->vectors[dof_n].ret(dof);
        memset(dof,0,x_dim*natms);
    }
    
    
    char* cdof;
    cdof_n=atoms->find_exist("cdof");
    
    if(cdof_n==-1)
    {
        cdof_n=atoms->add<char>(0,no_types,"cdof");
        atoms->vectors[cdof_n].ret(cdof);
        memset(cdof,0,no_types*natms);
    }
    
    atoms->vectors[dof_n].ret(dof);
    atoms->vectors[cdof_n].ret(cdof);
    
    
    int* id;
    id_n=atoms->find("id");
    atoms->vectors[id_n].ret(id);
    
    for(int i=0;i<natms;i++)
        for(int j=0;j<list_size;j++)
            if(id[i]==list[j])
            {
                memcpy(&dof[i*x_dim],dof_ch,x_dim);
                memcpy(&cdof[i*no_types],cdof_ch,no_types);
            }
    
    if(list_size)
        delete [] list;
    
    delete [] dof_ch;
    delete [] cdof_ch;
}


