

#ifndef __MAPP__md__
#define __MAPP__md__
#include <iostream>
#include "init.h"
#include "atoms.h"
namespace MAPP_NS {
    
    class MD : protected InitPtrs
    {
    private:
    protected:
    public:
        MD(MAPP *);
        virtual ~MD()=0;
        virtual void init()=0;
        virtual void fin()=0;
        virtual void run(int)=0;

        void add_dt(int,char**);
        void add_boltzmann(int,char**);
        void run(int,char**);
        TYPE0 dt,boltz;
    };
}

#endif