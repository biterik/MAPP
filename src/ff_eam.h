
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
        void set_fs(int,char**);
        int line_read(FILE*,char*&);
        
        int eam_mode;
        int nr,nrho;
        TYPE0 dr,drho,dr_inv,drho_inv,rho_max;
        
        void allocate();
        int allocated;
        void clean_up();
        
        void set_arrays();
        void interpolate(int,TYPE0,TYPE0**);
        
        
        int** type2rho;
        int** type2phi;
        
        int rho_n,f_n,type_n,x_n,dF_n;
        TYPE0* nrgy_strss;
        TYPE0 cut_sq;
        
        /*--------------------------------------------*/
        TYPE0* drhoi_dr;
        TYPE0* drhoj_dr;
        int max_pairs;
        /*--------------------------------------------*/
        
    public:
        ForceField_eam(MAPP *);
        ~ForceField_eam();
        void force_calc(int,TYPE0*);
        TYPE0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);

        void create_2nd_neigh_lst(){};
        TYPE0 calc_g(int,TYPE0,TYPE0*,TYPE0*){return 0.0;};
        void c_d_calc(){};
    };
    

    
}
#endif
#endif
