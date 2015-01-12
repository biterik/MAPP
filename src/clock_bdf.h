#ifdef Clock_Style
    ClockStyle(Clock_BDF,bdf)
#else
#ifndef __MAPP__clock_bdf__
#define __MAPP__clock_bdf__

#include <stdio.h>
#include "clock.h"
namespace MAPP_NS {
    class Clock_BDF :public Clock
    {
    private:
    protected:
        int c_n,c_d_n,dof_tot,dof_lcl;
        int max_iter,no_steps,max_order;
        TYPE0 min_gamma,gamma_red,slope;
        TYPE0 m_tol,a_tol,e_tol;
        TYPE0 min_del_t,max_del_t,initial_del_t;
        TYPE0 eq_ratio;
        
        
        TYPE0* y_0;
        
        // stuff for solution
        TYPE0* e_n;
        TYPE0* a;
        TYPE0* g0;
        TYPE0* g;
        TYPE0* h;
        TYPE0* c0;
        TYPE0 beta;
        // stuff for book keeping
        TYPE0* t;
        TYPE0** y;
        TYPE0* dy;
        
        // stuff for coefficient
        TYPE0 alpha_dy_0;
        TYPE0* alpha_y;
        TYPE0 dalpha_dy_0;
        TYPE0* dalpha_y;
        TYPE0* coef;
        
        int lo_ord_avail,hi_ord_avail;
        TYPE0 lo_err,hi_err,err;
        
        
        TYPE0 solve(TYPE0,int);
        int interpolate(TYPE0,int);
        void err_coef(TYPE0*,int);
        void err_est(int);
        void hi_lo_err_est(TYPE0,int);
        TYPE0* lwr_alpha;
        TYPE0 lwr_alpha_dy,lwr_alpha_y;
        
    public:
        Clock_BDF(MAPP *,int,char**);
        ~Clock_BDF();
        void run();
        void init();
        void fin();
    };
}

#endif 
#endif 