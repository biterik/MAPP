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
        int order;
        TYPE0* alpha;
        TYPE0 beta;
        TYPE0 ave_err;

        
        
        TYPE0 delta_t;
        VecLst* vecs_comm;
        
        
        TYPE0** y;
        TYPE0* a;
        TYPE0* g;
        TYPE0* c0;
        TYPE0* g0;
        TYPE0* h;
        int* y_list;
        int c_n,c_d_n;
        int tot_dim;
        void solve(TYPE0);
        void line_search();
        TYPE0 gamma_red,slope;
        int no_steps;
        int tmp_order;
        
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