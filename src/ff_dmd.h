#ifdef FF_Style
FFStyle(ForceField_DMD,DMD)
#else

#ifndef __MAPP__ff_dmd__
#define __MAPP__ff_dmd__

#include <iostream>
#include "ff.h"
#include "atoms.h"
namespace MAPP_NS
{
    class ForceField_DMD : public ForceField
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
        TYPE0 cut_sq_mod;
        TYPE0 rc,mod_rc,c_0,kbT;
        
        /*--------------------------------------------*/
        TYPE0* drhoi_dr;
        TYPE0* drhoj_dr;
        TYPE0* drhoi_dbeta;
        TYPE0* drhoj_dbeta;
        int max_pairs;
        /*--------------------------------------------*/
        
        
        TYPE0* xi;
        TYPE0* wi_0;
        TYPE0* wi_1;
        TYPE0* wi_2;
        TYPE0 beta_min,beta_max;
        int no_i;
        void set_weight_abs(int);
    public:
        ForceField_DMD(MAPP *);
        ~ForceField_DMD();
        void force_calc(int,TYPE0*);
        TYPE0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
    };
    
    
    
}


#endif 
#endif
