/*--------------------------------------------
 Created by Sina on 07/30/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__rand_engine__
#define __MAPP__rand_engine__

#include <iostream>
#include "init.h"
namespace MAPP_NS {
    class Random : protected InitPtrs{
    private:
        double reserved;
        int gauss_chk;
        int seed;
    protected:
    public:
        Random(MAPP *,int);
        ~Random();
        double uniform();
        double gaussian();
        double gaussian(double,double);
    };
}
#endif 

