/*--------------------------------------------
 Created by Sina on 2/5/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__xmath__
#define __MAPP__xmath__

#include <iostream>
#include "init.h"

using namespace std;
namespace MAPP_NS {
    class SPARSE : protected InitPtrs
    {
    private:
    protected:
    public:
    
        SPARSE(MAPP*,int,int);
        ~SPARSE();
        int no_elem;
        int h0,h1;
        int* idx0;
        int* idx1;
        double* mtrx;
        void add(int,int,double);
        void add(int,int*,int*,double*);
        void sort(int);
        void vectorize(int);
        
        int no_vec;
        int* idx_vec;
        int* init_pos;
        int* fin_pos;
        
    };
    
    // Parallel sparse decomposed by rows (h0)
    class SPARSE_P : protected InitPtrs
    {
    private:
    protected:
    public:
        SPARSE_P(MAPP*,int,int);
        ~SPARSE_P();
        
        int no_elem;
        int lcl_no_elem;
        int h0,h1;
        int lcl_h0;
        int lo_h0;
        int hi_h0;
        int totp;
        int myno;
        int* lcl_idx0;
        int* lcl_idx1;
        double* lcl_mtrx;
        void add(int,int,double);
        void addd(int,int,double);
        void add(int,int*,int*,double*);
        void sort(int);
        void vectorize(int);
        
        
        int lcl_no_vec;
        int* lcl_idx_vec;
        int* lcl_init_pos;
        int* lcl_fin_pos;
    };
    
    
    class SOLVEAXb : protected InitPtrs
    {
    private:
        SPARSE* A;
        int lsize;
        int tot_lsize;
        int totp;
        int myno;
        int lcl_lo,lcl_hi;
        
        int* lcl_idx0;
        int* lcl_idx1;
        double* lcl_mtrx;
        int lcl_no_elem;
        
        int* comm_snd_size;
        int* comm_rcv_size;
        int** comm_snd_lst;
        double** buff_snd;
        double* c;
        
        double** g;
        double* h;
        double d_sq;
    protected:
    public:
        SOLVEAXb(MAPP*,SPARSE*,double*,int);
        ~SOLVEAXb();
        void xchng(double*);
        void construct_ans();
        int solve(double);
        double* x;
        double* ans;
    };
    
    class COMB : protected InitPtrs
    {
        COMB(MAPP*);
        ~COMB();
        void clean();
        void comb(int*,int);
        void comb_rec(int,int,int,int,int*&,int*& ,int*,int,int);
        int** perm_list;
        int perm_list_size_0;
        int perm_list_size_1;
    };
    
}
#endif /* defined(__MAPP__xmath__) */
