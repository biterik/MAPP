/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/

#ifndef __MAPP__MAPP__
#define __MAPP__MAPP__
#include <iostream>
#include <mpi.h>
#include <stdio.h>

namespace MAPP_NS {
    enum {MD,VG,DMD};
    class MAPP {
    private:
    protected:
    public:
        MAPP(int,char**,MPI_Comm);
        ~MAPP();
        MPI_Comm world;
        int mode;
        int step_no;
        int step_tally;
        class Memory* memory;
        class Error* error;
        class PerAtom* peratom;
        class ForceField* forcefield;
        class Atoms* atoms;
        class Neighbor* neighbor;
        class AtomTypes* atom_types;
        class ThermoDynamics* thermo;
        class Write* write;
        class Clock* clock;
        FILE* output;
        
        FILE* input_file;
        
        class MD* md;
        class Min* min;
        
        int no_commands;
        
        void read_file();
        void command(char*);
        void min_style(int,char**);
        void ff_style(int,char**);
        void md_style(int,char**);
        void clock_style(int,char**);
        void read_style(int,char**);
        void write_style(int,char**);
        void command_style(int,char**);
        void change_mode(int,char**);
        int parse_line(char*,char**&);
        int hash_remover(char*,char*&);
        
        void test();
    };

}


#endif
