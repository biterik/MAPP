/*--------------------------------------------
 Created by Sina on 06/27/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include "atom_types.h"
#include "memory.h"
using namespace std;
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
AtomTypes::AtomTypes(MAPP* mapp):InitPtrs(mapp)
{
    no_types=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
AtomTypes::~AtomTypes()
{

    for(int i=0;i<no_types;i++)
        delete atom_names[i];
    if(no_types)
    {
        delete [] atom_names;
        delete [] mass;
    }
    
}
/*--------------------------------------------
 add a type
 --------------------------------------------*/
int AtomTypes::add_type(TYPE0 m,char* name)
{
    
    
    int type=0;
    for (int i=0;i<no_types;i++)
        if(!strcmp(name,atom_names[i]))
            type=i+1;
    
    if (type)
        return (type-1);
    
    
    
    GROW(atom_names,no_types,no_types+1);
    GROW(mass,no_types,no_types+1);
    int lngth= static_cast<int>(strlen(name))+1;
    CREATE1D(atom_names[no_types],lngth);
    for(int i=0;i<lngth;i++)
        atom_names[no_types][i]=name[i];
    mass[no_types]=m;
    
    no_types++;
    return (no_types-1);
    
    
   
}
/*--------------------------------------------
 find a type
 --------------------------------------------*/
int AtomTypes::find_type(char* name)
{
    int type=0;
    for (int i=0;i<no_types;i++)
        if(!strcmp(name,atom_names[i]))
            type=i+1;
    
    if(type)
        return (type-1);
    else
    {
        error->abort("atom type %s not found",name);
        return -1;
    }
}

