/*--------------------------------------------
 Created by Sina on 07/23/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifdef Read_Style
ReadStyle(ReadCFG,cfg)
#else
#ifndef __MAPP__readCFG__
#define __MAPP__readCFG__
#include <iostream>
#include "atoms.h"
#include "read.h"
namespace MAPP_NS {
    class ReadCFG : public Read{
    private:
        FILE* cfgfile;
        
        TYPE0** H0;
        TYPE0** eta;
        TYPE0** eta_sq;
        TYPE0** trns;
        TYPE0** H_x;
        TYPE0** H_x_d;
        
        char* line;
        
        TYPE0 basic_length;
        TYPE0 R;
        int entry_count;
        int ext_cfg;
        int header_cmplt;
        int vel_chk;
        int atom_cmplt;
        int last_type;
        
        void read_header();
        void read_atom();
        void read_atom_dmd();
        void M3sqroot(TYPE0**,TYPE0**);
        void set_box();
        void add_atom_read_x(int,TYPE0*);
        
        int type_n,x_n,x_d_n,id_n,c_n;
        int ch_x_d,ch_type,ch_id,ch_c;
        int dmd_no_types;
        
        int curr_id;
        char* ch_buff;
        VecLst* vec_list;
    protected:
    public:
        ReadCFG(MAPP *,int,char**);
        ~ReadCFG();
    };
}
#endif
#endif

