#ifndef __MAPP__clock__
#define __MAPP__clock__
#include "init.h"
#include "atoms.h"
#include <stdio.h>
namespace MAPP_NS {
    
    class Clock : protected InitPtrs
    {
    private:
    protected:
        VecLst* vecs_comm;

        ThermoDynamics* thermo;
        int pe_idx;
        int stress_idx;
        

        
    public:
        Clock(MAPP *);
        virtual ~Clock();
        void errors();
        virtual void run()=0;
        virtual void init()=0;
        virtual void fin()=0;
        
    };
    
    
}

#endif
