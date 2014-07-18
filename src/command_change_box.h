#ifdef Command_Style
    CommandStyle(CommandChangeBox,change_box)
#else
#ifndef __MAPP__command_change_box__
#define __MAPP__command_change_box__
#include <iostream>
#include "init.h"
namespace MAPP_NS {
    class CommandChangeBox: protected InitPtrs
    {
    private:
    protected:
    public:
        CommandChangeBox(MAPP*,int,char**);
        ~CommandChangeBox();
        
    };
}

#endif
#endif
