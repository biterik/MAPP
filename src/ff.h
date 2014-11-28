/*--------------------------------------------
 Created by Sina on 07/16/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__ff__
#define __MAPP__ff__
#include <iostream>
#include "init.h"
#include "vecmath.h"
#include "atoms.h"
namespace MAPP_NS {
    class ForceField : protected InitPtrs{
    private:
    protected:
    public:
        ForceField(MAPP *);
        virtual ~ForceField();
        virtual void force_calc(int,double*){};
        virtual TYPE0 energy_calc()=0;
        virtual void init(){};
        virtual void fin(){};
        virtual void coef(int,char**){};
        TYPE0* cut_sq;
        TYPE0* cut_sk_sq;
        
        virtual void create_2nd_neigh_lst(){};
        virtual TYPE0 calc_g(int,TYPE0,TYPE0*,TYPE0*)=0;
        virtual void calc_y(){};
    };
}
#endif

