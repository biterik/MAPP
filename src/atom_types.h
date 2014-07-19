/*--------------------------------------------
 Created by Sina on 06/27/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__atom_types__
#define __MAPP__atom_types__

#include <iostream>
#include "init.h"
#include "atoms.h"

namespace MAPP_NS {
    class AtomTypes : protected InitPtrs
    {
    private:
    protected:
    public:
        AtomTypes(MAPP *);
        ~AtomTypes();
        
        int no_types;
        char** atom_names;
        TYPE0* mass;
        
        int add_type(TYPE0,char*);
        int find_type(char*);
    };
}
#endif /* defined(__MAPP__atom_types__) */
