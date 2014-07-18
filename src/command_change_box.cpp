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
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            A[i][j]=1.0;
    
    int iarg=1;
    
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
    }
    
    TYPE0** H=atoms->H;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            H[i][j]*=A[i][j];
    
    for(int i=0;i<dim;i++)
        delete [] A[i];
    if(dim)
        delete [] A;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
CommandChangeBox::~CommandChangeBox()
{
    
}
