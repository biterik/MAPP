
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
        TYPE0*** F_arr;
        TYPE0*** phi_r_arr;
        TYPE0*** rho_arr;
        
        void set_setfl(int,char**);
        void set_funcfl(int,char**);
        int line_read(FILE*,char*&);
        
        int eam_mode;
        int nr,nrho;
        TYPE0 dr,drho;
        
        void allocate();
        int allocated;
        void clean_up();
        
    public:
        ForceField_eam(MAPP *);
        ~ForceField_eam();
        void force_calc(int,TYPE0*);
        TYPE0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
    };
    

    
}
#endif
#endif
