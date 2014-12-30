#ifdef Clock_Style
    ClockStyle(Clock_Adams,adams)
#else
#ifndef __MAPP__clock_adams__
#define __MAPP__clock_adams__

#include <stdio.h>
#include "clock.h"
namespace MAPP_NS {
    class Clock_Adams :public Clock
    {
    private:
    protected:
        int order;
        TYPE0* beta;
        TYPE0 ave_err;
        
        
        
        TYPE0 delta_t;
        VecLst* vecs_comm;
        
        
        TYPE0** y;
        TYPE0* a;
        TYPE0* g;
        TYPE0* c0;
        TYPE0* g0;
        TYPE0* h;
        int c_n,c_d_n;
        int tot_dim;
        TYPE0 solve(TYPE0);
        void line_search();
        TYPE0 gamma_red,slope;
        int no_steps;
        
        TYPE0* t;
        TYPE0* beta_mod;
        TYPE0 tot_time;
        void interpolate(TYPE0);
        
        
    public:
        Clock_Adams(MAPP *,int,char**);
        ~Clock_Adams();
        void run();
        void init();
        void fin();
        
    };
}

#endif
#endif
