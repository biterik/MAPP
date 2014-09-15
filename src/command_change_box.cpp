#include "command_change_box.h"
#include "atoms.h"
#include <stdlib.h>
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
CommandChangeBox::CommandChangeBox(MAPP* mapp
,int narg,char** args):InitPtrs(mapp)
{
    
    int dim=atoms->dimension;
    TYPE0** A;
    CREATE2D(A,dim,dim);
    
    int iarg=1;
    
    if(strcmp(args[iarg],"strain")==0)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                A[i][j]=0.0;
    }
    else if(strcmp(args[iarg],"dilation")==0)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                A[i][j]=1.0;
    }
    else if(strcmp(args[iarg],"equal")==0)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                A[i][j]=atoms->H[i][j];
    }
    else
        error->abort("wrong change box command");
    
    iarg++;
    
    
    int icmp,jcmp;
    while(iarg<narg)
    {
        if(sscanf(args[iarg],"H[%d][%d]",&icmp,&jcmp)==2)
        {
            
            if(icmp<0 || icmp>=atoms->dimension)
                error->abort("wrong command %s",args[iarg]);
            if(jcmp<0 || jcmp>=atoms->dimension)
                error->abort("wrong command %s",args[iarg]);
            
            iarg++;
            if(iarg==narg)
                error->abort("change box component ratio missing");
            
            if(icmp<=jcmp)
                A[jcmp][icmp]=atof(args[iarg]);
            else
                A[icmp][jcmp]=atof(args[iarg]);
            iarg++;
            
        }
        else
            error->abort("wrong command %s",args[iarg]);
    }
    
    
    
    if(strcmp(args[1],"strain")==0)
    {
        for(int i=0;i<dim;i++)
            A[i][i]+=1.0;
        TYPE0** C;
        CREATE2D(C,dim,dim);
        
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                C[i][j]=0.0;
        
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                for(int k=0;k<dim;k++)
                    C[i][j]+=atoms->H[i][k]*A[k][j];
        
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                atoms->H[i][j]=C[i][j];
        
        for(int i=0;i<dim;i++)
            delete [] C[i];
        if(dim)
            delete [] C;
    }
    else if(strcmp(args[1],"dilation")==0)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                atoms->H[i][j]*=A[i][j];
    }
    else if(strcmp(args[1],"equal")==0)
    {
        for(int i=0;i<dim;i++)
            for(int j=0;j<dim;j++)
                atoms->H[i][j]=A[i][j];
    }
    
    
    for(int i=0;i<dim;i++)
        delete [] A[i];
    if(dim)
        delete [] A;
    
    if(dim==3)
        M3INV_TRI_LOWER(atoms->H,atoms->B);
    else
        atoms->invert_lower_triangle(atoms->H,atoms->B,dim);
    

}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
CommandChangeBox::~CommandChangeBox()
{
    
}
