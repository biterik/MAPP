#ifdef Write_Style
    WriteStyle(WriteCFG,cfg)
#else
#ifndef __MAPP__write_cfg__
#define __MAPP__write_cfg__

#include "write.h"
#include "atoms.h"
#include <iostream>
namespace MAPP_NS {
    class WriteCFG: public Write
    {
    private:
    protected:
        char* file_name;
        int* vec_list;
        int no_vecs;
        int id_xst;
        int id_cmp;
        int type_cmp;
        int sorting;
        int tot_dim;
        int c_n;
        void write_file_dmd(int);
        void write_file_md(int);        
    public:
        WriteCFG(MAPP *,int,char**);
        ~WriteCFG();
        void write_file(int);
    };
    
}
#endif
#endif
