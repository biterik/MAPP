#ifdef Command_Style
    CommandStyle(CommandFix,fix)
#else
#ifndef __MAPP__command_fix__
#define __MAPP__command_fix__

#include <iostream>
#include "init.h"
#include "atoms.h"
namespace MAPP_NS {
    class CommandFix: protected InitPtrs
    {
    private:
    protected:
        void md(int,char**);
        void dmd(int,char**);
    public:
        CommandFix(MAPP*,int,char**);
        ~CommandFix();
    };
}

#endif 
#endif
