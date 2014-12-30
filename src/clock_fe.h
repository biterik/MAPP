#ifdef Clock_Style
    ClockStyle(Clock_FE,fe)
#else
#ifndef __MAPP__clock_fe__
#define __MAPP__clock_fe__

#include <stdio.h>
#include "clock.h"
namespace MAPP_NS {
    class Clock_FE :public Clock
    {
    private:
    protected:
        TYPE0 delta_t;
        VecLst* vecs_comm;
        int dof_tot,dof_lcl;
        int c_n,c_d_n;
        int no_steps;
        
        TYPE0 a_tol,e_tol;
        TYPE0 min_del_t;
        TYPE0 max_del_t;
        TYPE0 eq_ratio;
        TYPE0* y_0;
        TYPE0* dy_0;
        TYPE0* dy_1;
        
    public:
        Clock_FE(MAPP *,int,char**);
        ~Clock_FE();
        void run();
        void init();
        void fin();
        
    };
}
#endif
#endif 
