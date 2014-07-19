/*--------------------------------------------
 Created by Sina on 05/22/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#ifndef __MAPP__vecmath__
#define __MAPP__vecmath__
#include <stdio.h> 
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include "error.h"
#include "init.h"
#include "memory.h"
using namespace std;
namespace MAPP_NS {
    
    template <typename VECTYPE>
    class Vec : protected InitPtrs{
    private:
        
    protected:
    public:
        Vec(MAPP* mapp,int size) : InitPtrs(mapp)
        {
            try
            {
                vec=new VECTYPE[size];
            }
            catch(bad_alloc&)
            {
                printf("Error: allocation failure\n");
                exit(EXIT_FAILURE);
            }
            max_size=vec_size=size;
        }
        ~Vec()
        {
            delete [] vec;
        }
        VECTYPE* vec;
        
        int vec_size;
        int max_size;
        
        void resize(int size)
        {
            VECTYPE* new_vec;
            try
            {
                int new_size=MIN(size,vec_size);
                new_vec = new VECTYPE[size];
                memcpy(new_vec,vec,new_size*sizeof(VECTYPE));
            }
            catch (bad_alloc&) {
                printf("Error: allocation failure\n");
                exit(EXIT_FAILURE);
            }
            
            delete [] vec;
            vec = new_vec;
            max_size=vec_size=size;
        }
        
        void resize(int size,int msize)
        {
            if (msize<size)
            {
                printf("Error: max size smaller than vec size\n");
                exit(EXIT_FAILURE);
            }
            VECTYPE* new_vec;
            try
            {
                int new_size=MIN(size,vec_size);
                new_vec = new VECTYPE[msize];
                memcpy(new_vec,vec,new_size*sizeof(VECTYPE));
            }
            catch (bad_alloc&) {
                printf("Error: allocation failure\n");
                exit(EXIT_FAILURE);
            }
            
            delete [] vec;
            vec = new_vec;
            vec_size=size;
            max_size=msize;
        }
        
        
        VECTYPE& operator[](int indx)
        {
            if (indx>=vec_size)
            {
                printf("\nError: larger than size\n");
                exit(EXIT_FAILURE);
            }
                
            return vec[indx];
        }
        
    
    };
    
    
    
    template <typename LISTTYPE>
    class Lst : protected InitPtrs{
    private:
    protected:
    public:
        LISTTYPE* list;
        int* list_cur;
        
        int list_size;
        int no_elem;
        int grow_size;
        Lst(MAPP* mapp) : InitPtrs(mapp)
        {
            list_size=0;
            no_elem=0;
            grow_size=1;
        }
        Lst(MAPP* mapp,int lsize,int noelem)
        {
            CREATE1D(list,lsize);
            CREATE1D(list_cur,noelem);
            list_size=lsize;
            no_elem=noelem;
            grow_size=1;
        }
        ~Lst()
        {
            if(list_size)
                delete [] list;
            if(no_elem)
                delete [] list_cur;
        }
        void newlist(int lsize,int noelem)
        {
            if(list_size)
                delete [] list;
            if(no_elem)
                delete [] list_cur;
            
            CREATE1D(list,lsize);
            CREATE1D(list_cur,noelem);
            
            list_size=lsize;
            no_elem=noelem;
        }
        
        void resize(int lsize,int noelem)
        {

            GROW(list,list_size,lsize);
            GROW(list_cur,no_elem,noelem);
            list_size=lsize;
            no_elem=noelem;
        }
        void resize_elem_list(int noelem)
        {
            if(no_elem && no_elem!=noelem)
                delete [] list_cur;
            if(noelem!=no_elem)
            {
                CREATE1D(list_cur,noelem);
                no_elem=noelem;
            }
        }
        void grow_list(int lsize)
        {
            GROW(list,list_size,lsize);
            list_size=lsize;
        }
        
        LISTTYPE& operator[](int indx)
        {
            if (indx>=list_size-1)
            {
                grow_list(indx+grow_size);
            }
            
            return list[indx];
        }
        LISTTYPE& l_comp(int indx)
        {
            if (indx>=list_size)
            {
                grow_list(indx+grow_size);
            }
            
            return list[indx];
        }
        int& operator()(int indx)
        {
            if (indx>=no_elem)
            {
                printf("\nError: larger than no_elem\n");
                exit(EXIT_FAILURE);
            }
            
            return list_cur[indx];
        }
        int& c_comp(int indx)
        {
            if (indx>=no_elem)
            {
                printf("indx: %d noelem: %d\n",indx,no_elem);
                printf("\nError: larger than no_elem\n");
                exit(EXIT_FAILURE);
            }
            
            return list_cur[indx];
        }
    };
    
    
}
#endif
