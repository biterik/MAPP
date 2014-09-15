/*--------------------------------------------
 Created by Sina on 2/1/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__min__
#define __MAPP__min__

#include <iostream>
#include "init.h"
#include "atoms.h"
#include "line_search.h"
#include "thermo_dynamics.h"

namespace MAPP_NS {

    class Min : protected InitPtrs
    {
    private:
    protected:
        VecLst* vecs_comm;
        LineSearch* line_search;
        int chng_box;
        int dim;
        int err;
        TYPE0 curr_energy;
        ThermoDynamics* thermo;
        
        int pe_idx;
        int stress_idx;
        char* dof;
        int dof_n;
        int id_n;
    public:
        Min(MAPP *);
        virtual ~Min();
        void errors();
        virtual void run()=0;
        virtual void init()=0;
        virtual void fin()=0;
        void rectify_f(TYPE0*);
        int max_iter;
        TYPE0 energy_tolerance;
        int** H_dof;

    };


}
#endif /* defined(__MAPP__min__) */
