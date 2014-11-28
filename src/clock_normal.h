#ifdef Clock_Style
    ClockStyle(Clock_NORMAL,normal)
#else
#ifndef __MAPP__clock_normal__
#define __MAPP__clock_normal__

#include <stdio.h>
#include "clock.h"
namespace MAPP_NS {
    class Clock_NORMAL :public Clock
    {
    private:
    protected:
        TYPE0 delta_t;
        VecLst* vecs_comm;
        int c_n,c_d_n;
        int no_steps;

    public:
        Clock_NORMAL(MAPP *,int,char**);
        ~Clock_NORMAL();
        void run();
        void init();
        void fin();
        
    };
}

#endif 
#endif 
