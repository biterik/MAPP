#ifndef __MAPP__box_prop__
#define __MAPP__box_prop__
// TYPE0 can be only float, double, or long double
#define TYPE0 double
#define MPI_TYPE0 MPI_DOUBLE
//swap list growth size
#define SWAPGROWTH 10

#include <iostream>
#include <exception>
#include <typeinfo>
#include <cmath>
#include <mpi.h>
#include "init.h"
#include "memory.h"
using namespace std;
namespace MAPP_NS {
/*--------------------------------------------
 atomic vectors
 --------------------------------------------*/
    class AVec
    {
    private:
    protected:
    public:
        int type;
        int ph;
        int dim;
        int byte_size;
        int initiated;
        char* name;
        
        char* vec_0;
        unsigned char* vec_1;
        
        short int* vec_2;
        unsigned short int* vec_3;
        
        int* vec_4;
        unsigned int* vec_5;
        
        long int* vec_6;
        unsigned long int* vec_7;
        long long int* vec_8;
        unsigned long long int* vec_9;
        
        float* vec_10;
        double* vec_11;
        long double* vec_12;
        
        char* vec_0_dump;
        unsigned char* vec_1_dump;
        
        short int* vec_2_dump;
        unsigned short int* vec_3_dump;
        
        int* vec_4_dump;
        unsigned int* vec_5_dump;
        
        long int* vec_6_dump;
        unsigned long int* vec_7_dump;
        long long int* vec_8_dump;
        unsigned long long int* vec_9_dump;
        
        float* vec_10_dump;
        double* vec_11_dump;
        long double* vec_12_dump;
        
        void gather_dump();
        void print_dump(FILE*,int);
        
        Atoms* atms;
        
        ~AVec();
        void init(Atoms*,int,int,int,const char*);
        int unpack(char*&,int);
        int pack(char*&,int,int);
        
        int x_pack(char*&,int,int*,int);
        int x_unpack(char*&,int,int);
        
        void copy(int,int);
        void grow(int);
        void resize(int);
        void* ret_vec(int);
        void* ret_vec();
        void change_dimension(int);
        
        template <typename TYPE>
        TYPE* create(TYPE*&,int);
        
        template <typename TYPE>
        TYPE* grow(TYPE*&,int,int);
        
        void ret(char*&);
        void ret(unsigned char*&);
        
        void ret(short int*&);
        void ret(unsigned short int*&);
        
        void ret(int*&);
        void ret(unsigned int*&);
        
        void ret(long int*&);
        void ret(unsigned long int*&);
        void ret(long long int*&);
        void ret(unsigned long long int*&);
        
        void ret(float*&);
        void ret(double*&);
        void ret(long double*&);
        
        void ret_dump(char*&);
        void ret_dump(unsigned char*&);
        
        void ret_dump(short int*&);
        void ret_dump(unsigned short int*&);
        
        void ret_dump(int*&);
        void ret_dump(unsigned int*&);
        
        void ret_dump(long int*&);
        void ret_dump(unsigned long int*&);
        void ret_dump(long long int*&);
        void ret_dump(unsigned long long int*&);
        
        void ret_dump(float*&);
        void ret_dump(double*&);
        void ret_dump(long double*&);
        
        void del_dump();
    };
/*--------------------------------------------
 collection of all atomic vectors
 --------------------------------------------*/
    template <typename LISTTYPE>
    class SwapLst : protected InitPtrs{
    private:
    protected:
    public:
        
        LISTTYPE** list;
        int* list_size;
        int* pos;
        int grow_size;
        int no_swaps;
        //int* list_cur;
        
        SwapLst(MAPP* mapp):InitPtrs(mapp)
        {
            grow_size=1;
            no_swaps=0;
        }
        SwapLst(MAPP* mapp,int nswaps)
        {
            no_swaps=nswaps;
            CREATE1D(list,no_swaps);
            CREATE1D(list_size,no_swaps);
            CREATE1D(pos,no_swaps);
            grow_size=1;
            
            for(int i=0;i<no_swaps;i++)
                list_size[i]=pos[i]=0;
        }
        ~SwapLst()
        {
            if(no_swaps)
            {
                for(int i=0;i<no_swaps;i++)
                    if(list_size[i])
                        delete [] list[i];
                delete [] list;
                delete [] list_size;
                delete [] pos;
            }
            
        }
        void newlist(int nswaps)
        {
            if(no_swaps)
            {
                for(int i=0;i<no_swaps;i++)
                    if(list_size[i])
                        delete [] list[i];
                delete [] list;
                delete [] list_size;
                delete [] pos;
            }
            no_swaps=nswaps;
            CREATE1D(list,no_swaps);
            CREATE1D(list_size,no_swaps);
            CREATE1D(pos,no_swaps);
            grow_size=1;
            
            for(int i=0;i<no_swaps;i++)
                list_size[i]=pos[i]=0;
            
        }
        void add(int iswap,LISTTYPE val)
        {
            if(pos[iswap]==list_size[iswap])
            {
                GROW(list[iswap],list_size[iswap],list_size[iswap]+grow_size);
                list_size[iswap]+=grow_size;
            }
            list[iswap][pos[iswap]]=val;
            pos[iswap]++;
        }
        void reset()
        {
            for(int iswap=0;iswap<no_swaps;iswap++)
                pos[iswap]=0;
        }
        void reset(int iswap)
        {
            pos[iswap]=0;
        }
        
    };

/*--------------------------------------------
 properties of the box
 --------------------------------------------*/
    class Atoms:protected InitPtrs
    {
    private:
    protected:
        
        // snd and rcv buffers
        char* snd_buff_0;
        int snd_buff_0_capacity;
        char* snd_buff_1;
        int snd_buff_1_capacity;
        char* rcv_buff;
        int rcv_buff_capacity;
        
        char* snd_ph_buff;
        int snd_ph_buff_capacity;
        char* rcv_ph_buff;
        int rcv_ph_buff_capacity;
        
       
        void copy(AVec*,AVec*);
        void move(AVec*,AVec*);
        
        //void grow(int,int);
        void grow(int,int,class VecLst*);
        //void grow(int);
        
        void prime_factorize(int,int*&,int*&);
        
        int* comm_need;
        SwapLst<int>* ph_lst;        
        
    public:
        MPI_Comm comm_world;
        int dimension;
        int tot_natms;
        int natms;
        int natms_ph;
        int atm_vec_ph_size;
        int atm_vec_size;
        
        /* begining of box properties */
        // H matrix
        TYPE0** H;
        // Hinv matrix
        TYPE0** B;
        // lower bond of the local box
        TYPE0* s_lo;
        TYPE0* s_ph_lo;
        // higher bond of the local box
        TYPE0* s_hi;
        TYPE0* s_ph_hi;

        TYPE0 skin;
        TYPE0 tot_cut_ph;
        TYPE0* cut_ph_s;
        
        /* end of box box properties */
        Atoms(MAPP*,MPI_Comm);
        ~Atoms();
        void chng_dim(int);
        void set_ph(TYPE0);
        void set_ph(int);
        void store_0();
        void update_0(int,int,class VecLst*);
        
        // high level griding
        void hard_auto_grid_proc(TYPE0 f);
        void fac(int,int,int,int*&);
        void comb(int*,int);
        void comb_rec(int,int,int,int,int*&,int*&,int*,int,int);
    
        int** fac_list;
        int fac_list_size;
        
        TYPE0 energy_all;
        
        TYPE0** areas;
        TYPE0* vols;
        TYPE0 inter_n_efficency;
        
        int* res_perm;
        int* res_grid;
        
        int** nxt_p;
        int** prv_p;

        // end of high level griding
        
        //communication related parameters
        void auto_grid_proc();
        
        int tot_n;
        int my_n_no;
        int* p_per_n;
        int** n_p_grid;

        int tot_p;
        int* tot_p_grid;
        int my_p_no;
        int* my_loc;
        int** neigh_p;
        
        void ph_setup(int,class VecLst*);
        void ph_xchng(int,int,int*,int,class VecLst*);
        void reset_comm(class VecLst*);
        void reset_comm(int,class VecLst*);
        void xchng_comm(class VecLst*);
        void xchng_comm(int,class VecLst*);
        
        int pack(char*&,int,int&,int,class VecLst*);
        int unpack(char*&,int,int,class VecLst*);
        void x2s(int);
        void s2x(int);
        void x2s_no_correction(int);
        
        void x_unpack(char*&,class VecLst*,int);
        int x_pack(char*&,class VecLst*,int*,int);

        void x_unpack(char*&,int*,int,int);
        int x_pack(char*&,int*,int,int*,int);
        
        void update(class VecLst*);
        void update(int*,int,int);
        void update(int);
        void invert(TYPE0**,TYPE0**,int);
        void invert_lower_triangle(TYPE0**,TYPE0**,int);
        void add_skin(int,char**);
        /*-------------------------------*/
        AVec* vectors;
        int no_vecs;
        int tot_byte_size;
        
        void del(int);
        int find(const char* name);
        int find_exist(const char* name);
        template <typename TYPE>
        int add(int p,int d,const char* name)
        {

            for(int i=0;i<no_vecs;i++)
                if(!strcasecmp(name,vectors[i].name))
                    error->abort("Duplicate vec name");
            
            int type;
            
            if(no_vecs==0)
            {
                if (strcmp(typeid(TYPE).name(),typeid(TYPE0).name())!=0)
                    error->abort("wrong type, zeroth vector shoud be of type %s"
                           ,typeid(TYPE0).name());
                if (p!=1)
                    error->abort("zeroth vector shoud be phantom");
                if(d<dimension)
                    error->abort("dimension of zeroth vector should "
                           "be more or equal to box dimension");
            }
            
            if (strcmp(typeid(TYPE).name(),typeid(char).name())==0)
                type=0;
            else if (strcmp(typeid(TYPE).name(),typeid(unsigned char).name())==0)
                type=1;
            else if (strcmp(typeid(TYPE).name(),typeid(short int).name())==0)
                type=2;
            else if (strcmp(typeid(TYPE).name(),typeid(unsigned short int).name())==0)
                type=3;
            else if (strcmp(typeid(TYPE).name(),typeid(int).name())==0)
                type=4;
            else if (strcmp(typeid(TYPE).name(),typeid(unsigned int).name())==0)
                type=5;
            else if (strcmp(typeid(TYPE).name(),typeid(long int).name())==0)
                type=6;
            else if (strcmp(typeid(TYPE).name(),typeid(unsigned long int).name())==0)
                type=7;
            else if (strcmp(typeid(TYPE).name(),typeid(long long int).name())==0)
                type=8;
            else if (strcmp(typeid(TYPE).name(),typeid(unsigned long long int).name())==0)
                type=9;
            else if (strcmp(typeid(TYPE).name(),typeid(float).name())==0)
                type=10;
            else if (strcmp(typeid(TYPE).name(),typeid(double).name())==0)
                type=11;
            else if (strcmp(typeid(TYPE).name(),typeid(long double).name())==0)
                type=12;
            else
            {
                error->abort("the requested vector type is not provided");
                type=-1;
            }
            AVec* new_vectors;
            new_vectors = new AVec[no_vecs+1];
            for(int i=0;i<no_vecs;i++)
                new_vectors[i].initiated=0;
            
            for(int i=0;i<no_vecs;i++)
                move(&(new_vectors[i]),&(vectors[i]));
            
            new_vectors[no_vecs].init(this,type,p,d,name);
            
            if(no_vecs)
                delete [] vectors;
            vectors=new_vectors;
            no_vecs++;
            tot_byte_size=0;
            for(int i=0;i<no_vecs;i++)
                tot_byte_size+=vectors[i].byte_size;
            
            if(no_vecs==1)
            {
                char* x0name;
                CREATE1D(x0name,100);
                sprintf(x0name,"%s_0",name);
                add<TYPE>(0,dimension,x0name);
                delete [] x0name;
                return (no_vecs-2);
            }
    
            return (no_vecs-1);
        }
        
        void gather_all(VecLst*);
        /*
        template <typename TYPE>
        void ret(TYPE*& x,int ivec)
        {
            if(ivec>=no_vecs)
                error->abort("inconsistent type and vector number");
            
            if(vectors[ivec].type==0)
            {
                if(strcmp(typeid(TYPE).name(),typeid(char).name())==0)
                {
                    x=vectors[ivec].vec_0;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==1)
            {
                if(strcmp(typeid(TYPE).name(),typeid(unsigned char).name())==0)
                {
                    x=vectors[ivec].vec_1;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==2)
            {
                if(strcmp(typeid(TYPE).name(),typeid(short int).name())==0)
                {
                    x=vectors[ivec].vec_2;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==3)
            {
                if(strcmp(typeid(TYPE).name(),typeid(unsigned short int).name())==0)
                {
                    x=vectors[ivec].vec_3;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==4)
            {
                if(strcmp(typeid(TYPE).name(),typeid(int).name())==0)
                {
                    x=vectors[ivec].vec_4;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==5)
            {
                if(strcmp(typeid(TYPE).name(),typeid(unsigned int).name())==0)
                {
                    x=vectors[ivec].vec_5;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==6)
            {
                if(strcmp(typeid(TYPE).name(),typeid(long int).name())==0)
                {
                    x=vectors[ivec].vec_6;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==7)
            {
                if(strcmp(typeid(TYPE).name(),typeid(unsigned long int).name())==0)
                {
                    x=vectors[ivec].vec_7;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==8)
            {
                if(strcmp(typeid(TYPE).name(),typeid(long long int).name())==0)
                {
                    x=vectors[ivec].vec_8;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==9)
            {
                if(strcmp(typeid(TYPE).name(),typeid(unsigned long long int).name())==0)
                {
                    x=vectors[ivec].vec_9;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==10)
            {
                if(strcmp(typeid(TYPE).name(),typeid(float).name())==0)
                {
                    x=vectors[ivec].vec_10;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==11)
            {
                if(strcmp(typeid(TYPE).name(),typeid(double).name())==0)
                {
                    x=vectors[ivec].vec_11;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else if(vectors[ivec].type==12)
            {
                if(strcmp(typeid(TYPE).name(),typeid(long double).name())==0)
                {
                    x=vectors[ivec].vec_12;
                }
                else
                    error->abort("inconsistent type and vector number");
            }
            else
                error->abort("the requested vector type is not provided");
            
        }*/
        
    };
/*--------------------------------------------
 collection of all atomic vectors
 --------------------------------------------*/
    class VecLst: protected InitPtrs
    {
    private:
    protected:
    public:
        int* vec_list;
        int no_vecs;
        int byte_size;
        
        int* ph_vec_list;
        int ph_no_vecs;
        int ph_byte_size;
        
        int* update_every_ph_vec_list;
        int update_every_ph_no_vecs;
        int update_every_ph_byte_size;
        
        VecLst(MAPP*,int,...);
        VecLst(MAPP*,int*,int);
        ~VecLst();
        void add_update(int);
        void del_update(int);
        
    };

}
#endif /* defined(__MAPP__box_prop__) */
