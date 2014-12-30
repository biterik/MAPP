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
        TYPE0 min_del_t;
        TYPE0 alpha;
        TYPE0 delta_t;
        VecLst* vecs_comm;
        TYPE0* y_0;
        TYPE0* a;
        TYPE0* g;
        TYPE0* c0;
        TYPE0* g0;
        TYPE0* h;
        int c_n,c_d_n;
        int tot_dim,tot_dof;
        void line_search();
        TYPE0 gamma_red,slope;
        int no_steps;
        
        TYPE0 tot_time;
        TYPE0 RTOL,ATOL,min_gamma;
        int MAX_ITER;
        TYPE0 solve();
        TYPE0 step_size();
        TYPE0 delta_c;
        void prepare();
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
