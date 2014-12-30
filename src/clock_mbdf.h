#ifdef Clock_Style
    ClockStyle(Clock_MBDF,mbdf)
#else
#ifndef __MAPP__clock_mbdf__
#define __MAPP__clock_mbdf__

#include <stdio.h>
#include "clock.h"
namespace MAPP_NS {
    class Clock_MBDF :public Clock
    {
    private:
    protected:
        TYPE0 min_del_t;
        TYPE0 max_del_t,initial_del_t;
        TYPE0 beta;
        TYPE0 gamma_red,slope;
        
        int c_n,c_d_n;
        int dof_tot,dof_lcl;
        int no_steps;
        
        
        VecLst* vecs_comm;
        TYPE0** y;
        TYPE0* y_0;
        TYPE0* dy;
        TYPE0* a;
        TYPE0* g;
        TYPE0* c0;
        TYPE0* g0;
        TYPE0* h;
        
        
        TYPE0* t;
        TYPE0* mod_alpha;
        TYPE0* mod_d_alpha;
        TYPE0* alph_err;
        
        int interpolate(TYPE0,int);
        //TYPE0 error_clc(TYPE0,int);
        TYPE0 m_tol,a_tol,min_gamma,e_tol;
        TYPE0 eq_ratio;
        int max_iter;
        TYPE0 solve(TYPE0,int);
        int max_order;
        int fac(int);
        TYPE0 step_size(TYPE0,int);
        int err_calc(int,int,TYPE0,TYPE0&);
        TYPE0 err_est(int);
        TYPE0 err;
    public:
        Clock_MBDF(MAPP *,int,char**);
        ~Clock_MBDF();
        void run();
        void init();
        void fin();
        
    };
}

#endif 
#endif