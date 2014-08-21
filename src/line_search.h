#ifndef __MAPP__line_search__
#define __MAPP__line_search__

#include <iostream>
#include "atoms.h"
namespace MAPP_NS {
    enum {LS_S,LS_F_DOWNHILL,LS_F_GRAD0,LS_F_ALPHAMIN,MIN_F_MAX_ITER,MIN_F_TOLERANCE};
    
    class LineSearch : protected InitPtrs
    {
    private:
    protected:
        int x_n;
        int x_prev_n;
        int f_n;
        int dim;
        int x_dim;
        
        TYPE0 inner_f_h();
        TYPE0 inner_f_h_s();
        void normalize_h();
        void normalize_h_s();
        TYPE0 energy(TYPE0);
        VecLst* vecs_comm;
        TYPE0 d_max;
        TYPE0 s_max;
        TYPE0** N;
    public:
        TYPE0** H_prev;
        TYPE0** B_prev;
        TYPE0** h_H;
        TYPE0** f_H;
        int chng_box;
        LineSearch(MAPP *,VecLst*);
        virtual ~LineSearch()=0;
        virtual int line_min(TYPE0&,TYPE0&)=0;
        int h_n;
    };
    
    class LineSearch_BackTrack : public LineSearch
    {
    private:
    protected:
        TYPE0 c,rho,alpha_max,alpha_min;
    public:
        LineSearch_BackTrack(MAPP *,VecLst*);
        ~LineSearch_BackTrack();
        int line_min(TYPE0&,TYPE0&);
        int line_min_s(TYPE0&,TYPE0&);
    };

}

#endif /* defined(__MAPP__line_search__) */
