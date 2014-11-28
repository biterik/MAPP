#include "clock.h"
#include "error.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock::Clock(MAPP* mapp):InitPtrs(mapp)
{
    char** args;
    int narg=mapp->parse_line((char*)
    "PE S_xx S_yy S_zz S_yz S_zx S_xy",args);
    
    pe_idx=0;
    stress_idx=1;
    thermo=new ThermoDynamics(mapp,narg,args);
    for(int i=0;i<narg;i++)
        delete [] args[i];
    delete [] args;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock::~Clock()
{
    
}
