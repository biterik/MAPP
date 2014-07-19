
#include "error.h"
#include "md.h"
#include <stdlib.h>
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
MD::MD(MAPP* mapp):InitPtrs(mapp)
{
    if(forcefield==NULL)
        error->abort("force field is not initiated");

}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
MD::~MD()
{
    
}
/*--------------------------------------------
 add time step
 --------------------------------------------*/
void MD::add_dt(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong command");
    dt=atof(args[1]);
    if(dt<=0.0)
        error->abort("time step cannot be"
        " equal or less than zero");
}
/*--------------------------------------------
 add time step
 --------------------------------------------*/
void MD::add_boltzmann(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong command");
    boltz=atof(args[1]);
    if(boltz<=0.0)
        error->abort("boltzmann constant cannot"
                     " be equal or less than zero");
}
/*--------------------------------------------
 add time step
 --------------------------------------------*/
void MD::run(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong command");
    int steps=atoi(args[1]);
    
    if(steps<=0)
        error->abort("number of steps cannot "
        "be equal or less than zero");
    
    init();
    run(steps);
    fin();
    
}
