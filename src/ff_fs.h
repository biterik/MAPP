/*--------------------------------------------
 Created by Sina on 01/29/14.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifdef FF_Style
    FFStyle(ForceField_fs,FS)
#else
#ifndef __MAPP__ff_fs__
#define __MAPP__ff_fs__
#include "ff.h"
namespace MAPP_NS {
    class ForceField_fs : public ForceField{
    private:
        int x_n,f_n,type_n,rho_n;
        int arr_size;
        
        TYPE0** mat_t_1;
        TYPE0** mat_t_2;
        TYPE0* mat_A;
        
        TYPE0* cut_phi;
        TYPE0* cut_rho;
        TYPE0* mat_k_1;
        TYPE0* mat_k_2;
        TYPE0* mat_k_3;
        
        int* chk_coef;
        
        TYPE0* nrgy_strss;
        
        /*--------------------------------------------*/
        TYPE0* drhoi_dr;
        TYPE0* drhoj_dr;
        int max_pairs;
        /*--------------------------------------------*/
        
    protected:
    public:
        ForceField_fs(MAPP *);
        ~ForceField_fs();
        void force_calc(int,TYPE0*);
        void force_calc0(int,TYPE0*);
        TYPE0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
        int shift;
    };
}
#endif
#endif
