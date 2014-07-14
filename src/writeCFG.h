/*--------------------------------------------
 Created by Sina on 1/29/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__writeCFG__
#define __MAPP__writeCFG__
#include <iostream>
#include "init.h"
namespace MAPP_NS {
    class WriteCFG : protected InitPtrs{
    private:
        char** args;
        int ncomp,tot_dims;
        int* vec_nos;
    protected:
    public:
        WriteCFG(MAPP *,int,char**);
        ~WriteCFG();
        void write_file(int);
    };
}
#endif /* defined(__MAPP__writeCFG__) */
