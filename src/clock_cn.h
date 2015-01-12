#ifdef Clock_Style
    ClockStyle(Clock_CN,cn)
#else
#ifndef __MAPP__clock_cn__
#define __MAPP__clock_cn__

#include <stdio.h>
#include "clock.h"
namespace MAPP_NS {
    class Clock_CN :public Clock
    {
    private:
    protected:
        int c_n,c_d_n,dof_tot,dof_lcl;
        int max_iter,no_steps;
        TYPE0 min_gamma,gamma_red,slope;
        TYPE0 m_tol,a_tol,e_tol;
        TYPE0 min_del_t,max_del_t,initial_del_t;
        TYPE0 eq_ratio;

        TYPE0 beta;
        TYPE0 err;
        
        VecLst* vecs_comm;
        
        TYPE0* t;
        
        TYPE0** dy;
        
        TYPE0* y_0;
        TYPE0* y0;
        TYPE0* dy0;
        TYPE0* a;
        TYPE0* g;
        TYPE0* c0;
        TYPE0* g0;
        TYPE0* h;

        TYPE0 solve(TYPE0);
        int interpolate(TYPE0);
        
        
    public:
        Clock_CN(MAPP *,int,char**);
        ~Clock_CN();
        void run();
        void init();
        void fin();
        
    };
}
#endif
#endif 
