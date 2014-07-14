
#ifdef Min_Style
MinStyle(Min_LBFGS,l-bfgs)

#else
#ifndef __MAPP__min_lbfgs__
#define __MAPP__min_lbfgs__
#include <iostream>
#include "min.h"
namespace MAPP_NS {
    class Min_LBFGS :public Min
    {
    private:
    protected:
        int f_prev_n;
        int x_prev_n;
        
        int m_it;
        int x_dim;
        int f_n;
        int h_n;
        int type_n;
        int* s_list;
        int* y_list;
        
        TYPE0* rho;
        TYPE0* alpha;
        
        TYPE0*** H_y;
        TYPE0*** H_s;
        int* H_s_y_list;
        TYPE0** H_prev;
        TYPE0** B_prev;
        TYPE0** f_H_prev;
        TYPE0** f_H;
        TYPE0** h_H;
    public:
        Min_LBFGS(MAPP *,int,char**);
        ~Min_LBFGS();
        void run();
        void init();
        void fin();

    };
}

#endif
#endif
