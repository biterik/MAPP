/*--------------------------------------------
 Created by Sina on 01/29/14.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifdef FF_Style
FFStyle(ForceField_fsm,FSM)
#else
#ifndef __MAPP__ff_fsm__
#define __MAPP__ff_fsm__
#include "ff.h"
namespace MAPP_NS {
    class ForceField_fsm : public ForceField{
    private:
        int x_n,f_n,type_n,rho_n;
        int arr_size;
        
        TYPE0** mat_t_1;
        TYPE0** mat_t_2;
        TYPE0* mat_A;
        
        TYPE0* cut_sq_phi;
        TYPE0* cut_sq_rho;
        TYPE0* mat_k_1;
        TYPE0* mat_k_2;
        TYPE0* mat_k_3;
        
        int* chk_coef;
        
        TYPE0* nrgy_strss;
        
    protected:
    public:
        ForceField_fsm(MAPP *);
        ~ForceField_fsm();
        void force_calc(int,TYPE0*);
        TYPE0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);
        int shift;

        void create_2nd_neigh_lst(){};
        TYPE0 calc_g(int,TYPE0,TYPE0*,TYPE0*){return 0.0;};
        void calc_y(){};
    };
}
#endif
#endif
