#ifdef Min_Style
    MinStyle(Min_CG,cg)
#else
#ifndef __MAPP__min_cg__
#define __MAPP__min_cg__

#include <iostream>
#include "min.h"
namespace MAPP_NS {
    class Min_CG : public Min {
    private:
    protected:
        int f_prev_n;
        int x_prev_n;
        int x_dim;
        int f_n;
        int h_n;
        
        
        TYPE0** H_prev;
        TYPE0** B_prev;
        TYPE0** f_H_prev;
        TYPE0** f_H;
        TYPE0** h_H;
    public:
        Min_CG(MAPP *,int,char**);
        ~Min_CG();
        void run();
        void init();
        void fin();
        
    };

}
#endif
#endif
