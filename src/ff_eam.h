
#ifdef FF_Style
    FFStyle(ForceField_eam,EAM)
#else
#ifndef __MAPP__ff_eam__
#define __MAPP__ff_eam__

#include <iostream>
#include "ff.h"
#include "atoms.h"
namespace MAPP_NS
{
    class ForceField_eam : public ForceField
    {
    private:
    protected:
    public:
        ForceField_eam(MAPP *);
        ~ForceField_eam();
        void force_calc(int,TYPE0*);
        TYPE0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
    };
    
    class fl
    {
    private:
    protected:
        int nr;
        TYPE0 dr;
        int nrho;
        TYPE0 drho;
        int no_types;
    public:
        fl(int,TYPE0,int,TYPE0,int);
        void compute_phi(TYPE0,int,int,TYPE0&);
        void compute_phi_dphi(TYPE0,int,int,TYPE0&,TYPE0&);
        void compute_rho(TYPE0,int,int,TYPE0&);
        void compute_rho_drho(TYPE0,int,int,TYPE0&,TYPE0&);
    };
    
    class func_fl:public fl
    {
        func_fl(int,TYPE0,int,TYPE0,int);
    };
    class set_fl:public fl
    {
        set_fl(int,TYPE0,int,TYPE0,int);
    };
    
}
#endif
#endif
