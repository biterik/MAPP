
#include "neighbor.h"
#include "ff_eam.h"
#include "atom_types.h"
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 --------------------------------------------*/
ForceField_eam::
ForceField_eam(MAPP* mapp) : ForceField(mapp)
{
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
ForceField_eam::~ForceField_eam()
{

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
 constructor
 --------------------------------------------*/
fl::fl(int nnr,TYPE0 ddr,int nnrho,TYPE0 ddrho
,int nno_types)
{
    nr=nnr;
    dr=ddr;
    nrho=nnrho;
    drho=ddrho;
    no_types=nno_types;
}
/*--------------------------------------------
 constructor
 --------------------------------------------*/
func_fl::func_fl(int nnr,TYPE0 ddr,int nnrho,
TYPE0 ddrho,int nno_types)
:fl(nnr,ddr,nnrho,ddrho,nno_types)
{
    
}




