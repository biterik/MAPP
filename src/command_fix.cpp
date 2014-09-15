#include "command_fix.h"
#include <stdlib.h>
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
CommandFix::CommandFix(MAPP* mapp
,int narg,char** args):InitPtrs(mapp)
{
    int dof_n;
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
    else if(strcmp(args[1],"x")==0)
        dof_ch[0]=1;
    else if(strcmp(args[1],"y")==0)
        dof_ch[1]=1;
    else if(strcmp(args[1],"z")==0)
        dof_ch[2]=1;
    else if(strcmp(args[1],"xy")==0 ||
            strcmp(args[1],"yx")==0)
        dof_ch[0]=dof_ch[1]=1;
    else if(strcmp(args[1],"yz")==0 ||
            strcmp(args[1],"zy")==0)
        dof_ch[2]=dof_ch[3]=1;
    else if(strcmp(args[1],"zx")==0 ||
            strcmp(args[1],"xz")==0)
        dof_ch[2]=dof_ch[0]=1;
    else if(strcmp(args[1],"xyz")==0 ||
            strcmp(args[1],"xzy")==0 ||
            strcmp(args[1],"yzx")==0 ||
            strcmp(args[1],"yxz")==0 ||
            strcmp(args[1],"zxy")==0 ||
            strcmp(args[1],"zyx")==0)
        dof_ch[0]=dof_ch[1]=dof_ch[2]=1;
    else
        error->abort("incorrect fix command");

    FILE* fp=NULL;
    char* line;
    CREATE1D(line,MAXCHAR);
    int no_atoms;
    
    int* list=NULL;
    int list_size=0;
    
    for(int iarg=2;iarg<narg;iarg++)
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

    TYPE0* x;
    atoms->vectors[0].ret(x);
    
    for(int i=0;i<natms;i++)
        for(int j=0;j<list_size;j++)
            if(id[i]==list[j])
            {
                memcpy(&dof[i*x_dim],dof_ch,x_dim);
                //printf("%lf %lf %lf \n",x[i*x_dim],x[i*x_dim+1],x[i*x_dim+2]);
            }
    
    if(list_size)
        delete [] list;
    
    delete [] dof_ch;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
CommandFix::~CommandFix()
{
    
}
