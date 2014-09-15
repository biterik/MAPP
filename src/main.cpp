/*--------------------------------------------
 Created by Sina on 05/12/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "MAPP.h"
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
using namespace MAPP_NS;
using namespace std;
int main(int narg, char** arg)
{

    MPI_Init(&narg,&arg);
    MAPP* mapp = new MAPP(narg,arg,MPI_COMM_WORLD);
    delete mapp;
    MPI_Finalize();    
    
    return EXIT_SUCCESS;
}


