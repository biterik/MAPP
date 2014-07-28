/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__neighbor__
#define __MAPP__neighbor__
#include <iostream>
#include "init.h"
#include "atoms.h"
#include "memory.h"

namespace MAPP_NS {
    class Neighbor : protected InitPtrs{
    private:
    protected:
        TYPE0* s_tmp;
        int* bin_denom_list;
        TYPE0* bin_size;
        
        int* tot_bin_grid;
        int tot_bin;
        
        int** bin_neigh_list;
        int* bin_neigh_list_size;
        
        void create_bin_list();
        
        void bin_atoms();
        void bin_atoms_s();
        int x2bin(TYPE0*);
        int s2bin(TYPE0*);
        
        int* first_atom_bin;
        int first_atom_bin_size;
        
        int* atm_bin;
        int atm_bin_size;
        
        int* next_atm;
        int next_atm_size;
       
        int type_n;
        
        void find_bin_no(int,int,int*&,int,int*,int*);
    public:
        int pair_wise;
        int** neighbor_list;
        int* neighbor_list_size;
        int neighbor_list_size_size;
        int no_pairs;
        
        Neighbor(MAPP *);
        ~Neighbor();
        
        void create_list(int,int);
        void init();
        void fin();
        
    };
}

#endif
