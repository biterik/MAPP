#ifdef Command_Style
    CommandStyle(CommandDisplace,displace)
#else
#ifndef __MAPP__command_displace__
#define __MAPP__command_displace__

#include <iostream>
#include "init.h"
#include "atoms.h"
namespace MAPP_NS {
    class CommandDisplace: protected InitPtrs
    {
    private:
        void move(int*,int,TYPE0*);
        int id_n;
    protected:
    public:
        CommandDisplace(MAPP*,int,char**);
        ~CommandDisplace();
        
    };
}

#endif 
#endif
