/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/

#ifndef __MAPP__error__
#define __MAPP__error__
#include <iostream>
#include "init.h"
namespace MAPP_NS {
    class Error : protected InitPtrs {
    private:
    protected:
        int my_no;
    public:
        Error(MAPP *);
        ~Error();
        void abort(const char*,...);
        void warning(const char*,...);
    };
}
#endif
