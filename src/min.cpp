/*--------------------------------------------
 Created by Sina on 2/1/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include "min.h"
#include "ff.h"
#include "neighbor.h"
#include "error.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::Min(MAPP* mapp):InitPtrs(mapp)
{
    if(forcefield==NULL)
        error->abort("force field is not initiated");
    
    
    chng_box=0;
    dim=atoms->dimension;
    err=LS_S;
    
    char** args;
    int narg=mapp->parse_line((char*)
    "PE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    pe_idx=0;
    stress_idx=1;
    thermo=new ThermoDynamics(mapp,narg,args);
    for(int i=0;i<narg;i++)
        delete [] args[i];
    delete [] args;
    
    // default values 
    max_iter=10000;
    energy_tolerance=1.0e-8;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Min::~Min()
{
    
}
/*--------------------------------------------
 error messages
 --------------------------------------------*/
void Min::errors()
{
    
    if(err==LS_F_DOWNHILL)
    {
        if(atoms->my_p_no==0)
            printf("line search failed: not downhill direction\n\n");
    }
    else if(err==LS_F_GRAD0)
    {
        if(atoms->my_p_no==0)
            printf("line search failed: gradient is zero\n\n");
    }
    else if(err==LS_F_ALPHAMIN)
    {
        if(atoms->my_p_no==0)
            printf("line search failed: reached minimum alpha\n\n");
    }
    else if (err==MIN_F_MAX_ITER)
    {
        if(atoms->my_p_no==0)
            printf("minimization finished: maximum iteration reached\n\n");
    }
    else if (err==MIN_F_TOLERANCE)
    {
        if(atoms->my_p_no==0)
            printf("minimization finished: energy tolerance reached\n\n");
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Min::rectify_f(TYPE0* f)
{
    if(dof_n==-1)
        return;
    atoms->vectors[dof_n].ret(dof);
    int tot=(atoms->natms)*(atoms->vectors[0].dim);
    for(int i=0;i<tot;i++) if(dof[i]==1) f[i]=0.0;
}

