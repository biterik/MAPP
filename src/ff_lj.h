/*--------------------------------------------
 Created by Sina on 07/15/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifdef FF_Style
    FFStyle(ForceField_lj,LJ)
#else
#ifndef __MAPP__ff_lj__
#define __MAPP__ff_lj__
#include "ff.h"
#include "atoms.h"
namespace MAPP_NS {
    class ForceField_lj : public ForceField{
    private:
        int x_n,f_n,type_n;
        int arr_size;
        TYPE0* sigma;
        TYPE0* epsilon;
        TYPE0* offset;
        int* chk_coef;
        TYPE0* nrgy_strss;
    protected:
    public:
        ForceField_lj(MAPP *);
        ~ForceField_lj();
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
