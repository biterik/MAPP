/*--------------------------------------------
 Created by Sina on 2/5/14.
 Copyright (c) 2014 MIT. All rights reserved.
 --------------------------------------------*/
#include <cmath>
#include "xmath.h"
#include "memory.h"
#include "error.h"
#include "atoms.h"

using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
SPARSE::
SPARSE(MAPP* mapp,int n0,int n1): InitPtrs(mapp)
{
    h0=n0;
    h1=n1;
    no_vec=0;
    no_elem=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
SPARSE::~SPARSE()
{
    if(no_elem)
    {
        delete [] idx0;
        delete [] idx1;
        delete [] mtrx;
    }
    
    if(no_vec)
    {
        delete [] fin_pos;
        delete [] idx_vec;
        delete [] init_pos;
    }
    
}
/*--------------------------------------------
 add an element to the sparse
 --------------------------------------------*/
void SPARSE::add(int i0,int i1,double elem)
{
    if(i0>=h0 || i1>=h1 || i0<0 || i1<0)
        error->abort("in sparse the index "
                     "is higher than the number "
                     "columns/rows");
    
    if (elem==0.0)
        return;
    
    int no=0;
    while (no<no_elem && (idx0[no]!=i0
           || idx1[no]!=i1))
        no++;
    
    if(no<no_elem)
    {
        mtrx[no]=elem;
        return;
    }
    
    GROW(idx0,no_elem,no_elem+1);
    GROW(idx1,no_elem,no_elem+1);
    GROW(mtrx,no_elem,no_elem+1);
    
    idx0[no_elem]=i0;
    idx1[no_elem]=i1;
    mtrx[no_elem]=elem;
    no_elem++;
}
/*--------------------------------------------
 add elements
 --------------------------------------------*/
void SPARSE::add(int size,int* i0_lst,int* i1_lst,double* elem_lst)
{
    
    double elem;
    int i0,i1;
    for(int i=0; i<size;i++)
    {
        elem=elem_lst[i];
        i0=i0_lst[i];
        i1=i1_lst[i];
        
        if (elem==0.0)
            continue;

        
        int no=0;
        while (no<no_elem && (idx0[no]!=i0
                || idx1[no]!=i1))
            no++;
        
        if(no<no_elem)
        {
            mtrx[no]=elem;
            continue;
        }
        
        GROW(idx0,no_elem,no_elem+1);
        GROW(idx1,no_elem,no_elem+1);
        GROW(mtrx,no_elem,no_elem+1);
        
        idx0[no_elem]=i0;
        idx1[no_elem]=i1;
        mtrx[no_elem]=elem;
        no_elem++;
    }
}
/*--------------------------------------------
 sort the sparse matrix row/column wise
 --------------------------------------------*/
void SPARSE::sort(int r_c)
{
    int tmp0,tmp1;
    double tmpx;
    if(r_c==0)
    {
        for(int i=0;i<no_elem;i++)
        {
            for(int j=i+1;j<no_elem;j++)
            {
                if(idx0[j]<idx0[i])
                {
                    tmp0=idx0[j];
                    tmp1=idx1[j];
                    tmpx=mtrx[j];
                    
                    idx0[j]=idx0[i];
                    idx1[j]=idx1[i];
                    mtrx[j]=mtrx[i];
                    
                    idx0[i]=tmp0;
                    idx1[i]=tmp1;
                    mtrx[i]=tmpx;
                }
            }
        }
    }
    else if(r_c==1)
    {
        for(int i=0;i<no_elem;i++)
        {
            for(int j=i+1;j<no_elem;j++)
            {
                if(idx1[j]<idx1[i])
                {
                    tmp0=idx0[j];
                    tmp1=idx1[j];
                    tmpx=mtrx[j];
                    
                    idx0[j]=idx0[i];
                    idx1[j]=idx1[i];
                    mtrx[j]=mtrx[i];
                    
                    idx0[i]=tmp0;
                    idx1[i]=tmp1;
                    mtrx[i]=tmpx;
                }
            }
        }
    }
    
    int strt=0;
    int end;
    int pos=strt;
    
    if (r_c==0)
    {
        while(pos<no_elem)
        {
            while(idx0[pos]==idx0[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            end=pos;
            for (int j=strt;j<end;j++)
            {
                for (int k=j+1;k<end;k++)
                {
                    if(idx1[j]>idx1[k])
                    {
                        tmp1=idx1[k];
                        tmpx=mtrx[k];
                        
                        idx1[k]=idx1[j];
                        mtrx[k]=mtrx[j];
                        
                        idx1[j]=tmp1;
                        mtrx[j]=tmpx;
                    }
                }
            }
            strt=pos;
        }
    }
    else if (r_c==1)
    {
        while(pos<no_elem)
        {
            while(idx1[pos]==idx1[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            end=pos;
            for (int j=strt;j<end;j++)
            {
                for (int k=j+1;k<end;k++)
                {
                    if(idx0[j]>idx0[k])
                    {
                        tmp1=idx0[k];
                        tmpx=mtrx[k];
                        
                        idx0[k]=idx0[j];
                        mtrx[k]=mtrx[j];
                        
                        idx0[j]=tmp1;
                        mtrx[j]=tmpx;
                    }
                }
            }
            strt=pos;
        }
    }
}
/*--------------------------------------------
 vectorize the matrix row/column wise
 --------------------------------------------*/
void SPARSE::vectorize(int r_c)
{
    sort(r_c);
    
    if(no_vec)
    {
        delete [] fin_pos;
        delete [] idx_vec;
        delete [] init_pos;
    }
    
    
    no_vec=0;
    int strt=0;
    int pos=strt;
    if (r_c==0)
    {
        while(pos<no_elem)
        {
            while(idx0[pos]==idx0[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            GROW(init_pos,no_vec,no_vec+1);
            GROW(fin_pos,no_vec,no_vec+1);
            GROW(idx_vec,no_vec,no_vec+1);
            
            init_pos[no_vec]=strt;
            fin_pos[no_vec]=pos;
            idx_vec[no_vec]=idx0[strt];
            
            no_vec++;
            strt=pos;
        }
    }
    else if (r_c==1)
    {
        while(pos<no_elem)
        {
            while(idx1[pos]==idx1[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            
            GROW(init_pos,no_vec,no_vec+1);
            GROW(fin_pos,no_vec,no_vec+1);
            GROW(idx_vec,no_vec,no_vec+1);
            
            init_pos[no_vec]=strt;
            fin_pos[no_vec]=pos;
            idx_vec[no_vec]=idx1[strt];
            
            no_vec++;
            strt=pos;
        }
    }
}

/*--------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 --------------------------------------------*/
/*--------------------------------------------
 constructor: uses CG to solve Ax=b
 --------------------------------------------*/
SOLVEAXb::SOLVEAXb(MAPP* mapp,SPARSE* Mat,
double* b,int size): InitPtrs(mapp)
{
    A=Mat;

    if(A->h1!=size)
        error->abort("Dimensions of A and b does not match");

    myno=atoms->my_p_no;
    totp=atoms->tot_p;
    if(size<totp)
    error->abort("Size of vector cannot be "
        "less than number processor");
    
    
    lsize=size/totp;
    
    lcl_lo=myno*lsize;
    if(size%totp!=0)
        if(myno==totp-1)
            lsize+=size%totp;
    lcl_hi=lcl_lo+lsize;
    
    A->vectorize(1);
    
    int no_vec=A->no_vec;
    int* idx_vec=A->idx_vec;
    int* init_pos=A->init_pos;
    int* fin_pos=A->fin_pos;
    
    int pos=0;
    while (idx_vec[pos]<lcl_lo)
    {
        pos++;
        if(pos==no_vec)
            break;
    }
    int idx_pos_lo=pos;
    
    while (idx_vec[pos]<lcl_hi)
    {
        pos++;
        if(pos==no_vec)
            break;
    }
    int idx_pos_hi=pos;
    

    
    int i1,crsi0,crsi1;
    int comp0,comp1;
    double val;
    tot_lsize=lsize;
    int comm_idx;

    CREATE1D(comm_snd_size,totp-1);
    CREATE1D(comm_snd_lst,totp-1);
    CREATE1D(comm_rcv_size,totp-1);
    for(int i=0;i<totp-1;i++)
        comm_snd_size[i]=comm_rcv_size[i]=0;
    
    lcl_no_elem=0;
    for(int i0=idx_pos_lo;i0<idx_pos_hi;i0++)
    {
        for(int j=0;j<no_vec;j++)
        {
            i1=j+idx_pos_lo;
            if(i1>=no_vec)
                i1-=no_vec;
            val=0.0;
            
            crsi0=init_pos[i0];
            crsi1=init_pos[i1];
            while(crsi0<fin_pos[i0])
            {
                while(crsi1<fin_pos[i1] && A->idx0[crsi1]<A->idx0[crsi0])
                    crsi1++;
                if(crsi1<fin_pos[i1] && A->idx0[crsi1]==A->idx0[crsi0])
                    val+=A->mtrx[crsi0]*A->mtrx[crsi1];
                crsi0++;
            }
            if (val!=0.0)
            {
                comp0=idx_vec[i0]-lcl_lo;
                if (idx_vec[i1]>=lcl_lo && idx_vec[i1]<lcl_hi)
                {
                    comp1=idx_vec[i1]-lcl_lo;
                    
                    GROW(lcl_idx0,lcl_no_elem,lcl_no_elem+1);
                    GROW(lcl_idx1,lcl_no_elem,lcl_no_elem+1);
                    GROW(lcl_mtrx,lcl_no_elem,lcl_no_elem+1);
                    
                    lcl_idx0[lcl_no_elem]=comp0;
                    lcl_idx1[lcl_no_elem]=comp1;
                    lcl_mtrx[lcl_no_elem]=val;
                }
                else
                {
                    /*
                    comm_idx=0;
                    while (comm_idx<totp && idx_vec[i1]<comp_lim[comm_idx])
                        comm_idx++;
                    comm_idx--;*/
                    comm_idx=(idx_vec[i1]/(size/totp));
                    if (comm_idx>totp-1)
                        comm_idx=totp-1;
                    comm_idx-=myno+1;
                    if(comm_idx<0)
                        comm_idx+=totp;
                    
                    GROW(comm_snd_lst[comm_idx],comm_snd_size[comm_idx],comm_snd_size[comm_idx]+1);
                    comm_snd_lst[comm_idx][comm_snd_size[comm_idx]]=comp0;
                    comm_snd_size[comm_idx]++;
                    comm_rcv_size[totp-2-comm_idx]++;
                    comp1=tot_lsize;
                    GROW(lcl_idx0,lcl_no_elem,lcl_no_elem+1);
                    GROW(lcl_idx1,lcl_no_elem,lcl_no_elem+1);
                    GROW(lcl_mtrx,lcl_no_elem,lcl_no_elem+1);
                    
                    lcl_idx0[lcl_no_elem]=comp0;
                    lcl_idx1[lcl_no_elem]=comp1;
                    lcl_mtrx[lcl_no_elem]=val;
                    
                    tot_lsize++;
                }
                lcl_no_elem++;
            }
        }
    }
    /*
    int sndp,rcvp;
    MPI_Status status;
    
    for(int i=0;i<totp-1;i++)
    {
        rcvp=myno+i+1;
        sndp=myno-i-1;

        if(rcvp>=totp)
            rcvp-=totp;
        if(sndp<0)
            sndp+=totp;
        
        MPI_Sendrecv(&comm_snd_size[totp-2-i],1,
        MPI_INT,sndp,0,&comm_rcv_size[i],
                     1,MPI_INT,rcvp,0,world,
                     &status);
    }
    */

    CREATE1D(buff_snd,totp-1);
    for(int i=0;i<totp-1;i++)
    {
        CREATE1D(buff_snd[i],comm_snd_size[i]);
    }
 
    CREATE1D(c,lsize);
 
    for(int i=0;i<lsize;i++)
        c[i]=0.0;

    for(int i=idx_pos_lo;i<idx_pos_hi;i++)
        for(int j=init_pos[i];j<fin_pos[i];j++)
            c[idx_vec[i]-lcl_lo]+=b[A->idx0[j]]*(A->mtrx[j]);
    

    d_sq=0.0;
    for(int i=0;i<size;i++)
        d_sq+=b[i]*b[i];

    
    CREATE1D(x,tot_lsize);
    CREATE1D(g,2);
    CREATE1D(g[0],lsize);
    CREATE1D(g[1],lsize);
    CREATE1D(h,tot_lsize);
    
    for(int i=0;i<tot_lsize;i++)
        x[i]=0.0;
    for(int i=lcl_lo;i<lcl_hi;i++)
    {
        //x[i-lcl_lo]=b[i];
        x[i-lcl_lo]=0.0;
    }
    xchng(x);

    CREATE1D(ans,size);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
SOLVEAXb::~SOLVEAXb()
{
    if(totp-1)
    {
        for(int i=0;i<totp-1;i++)
        {
            if(comm_snd_size[i])
            {
                delete [] buff_snd[i];
                delete [] comm_snd_lst[i];
            }
        }
        
        delete [] buff_snd;
        delete [] comm_snd_lst;
        delete [] comm_rcv_size;
        delete [] comm_snd_size;
    }
    if(lcl_no_elem)
    {
        delete [] lcl_idx0;
        delete [] lcl_idx1;
        delete [] lcl_mtrx;
    }
    if(lsize)
    {
        delete [] x;
        delete [] h;
        delete [] c;
        delete [] g[0];
        delete [] g[1];
        delete [] g;
    }
    delete [] ans;
}
/*--------------------------------------------
packing
--------------------------------------------*/
void SOLVEAXb::xchng(double* vec)
{
    MPI_Request request;
    MPI_Status status;
    int rcvp,sndp;
    int curs=lsize;
    for(int i=0;i<totp-1;i++)
    {
        rcvp=myno-i-1;
        sndp=myno+i+1;
        
        if(sndp>=totp)
            sndp-=totp;
        if(rcvp<0)
            rcvp+=totp;
        for(int j=0;j<comm_snd_size[i];j++)
            buff_snd[i][j]=vec[comm_snd_lst[i][j]];
        if (comm_rcv_size[i])
        {
            MPI_Irecv(&vec[curs],comm_rcv_size[i],
                      MPI_DOUBLE,rcvp,0,world,
                      &request);
        }
        if (comm_snd_size[i])
            MPI_Send(&buff_snd[i][0],comm_snd_size[i],
                     MPI_DOUBLE,sndp,0,world);
        if (comm_rcv_size[i])
            MPI_Wait(&request,&status);
        curs+=comm_rcv_size[i];
    }
}
/*--------------------------------------------
 packing
 --------------------------------------------*/
int SOLVEAXb::solve(double tol)
{
    int niter=A->h1+1;
    int comp;
    double g1_sq_l,g01_sq_l,hBh_l,hg_l,xBx_l,cx_l;
    double g0_sq=0.0,g1_sq,g01_sq,hBh,hg,xBx,cx,coef;
    double error=0.0,alpha;
    
    for(int i=0;i<lsize;i++)
        h[i]=0.0;
    for(int iter=0;iter<niter;iter++)
    {
        
        if(iter)
        {
            comp=iter%2;
            for(int i=0;i<lsize;i++)
                g[comp][i]=0.0;
            
            for(int i=0;i<lcl_no_elem;i++)
                g[comp][lcl_idx0[i]]-=2.0*lcl_mtrx[i]*x[lcl_idx1[i]];
            
            
            g1_sq=g1_sq_l=0.0;
            g01_sq=g01_sq_l=0.0;
            for(int i=0;i<lsize;i++)
            {
                g[comp][i]+=2*c[i];
                g1_sq_l+=g[comp][i]*g[comp][i];
                g01_sq_l+=g[1-comp][i]*g[comp][i];
            }
            
            MPI_Allreduce(&g1_sq_l,&g1_sq,1,MPI_DOUBLE,MPI_SUM,world);
            MPI_Allreduce(&g01_sq_l,&g01_sq,1,MPI_DOUBLE,MPI_SUM,world);
            
            
            coef=(g1_sq-g01_sq)/g0_sq;
            for(int i=0;i<lsize;i++)
            {
                h[i]*=coef;
                h[i]+=g[comp][i];
            }
        }
        else
        {
            comp=iter%2;
            for(int i=0;i<lsize;i++)
                h[i]=0.0;
            
            for(int i=0;i<lcl_no_elem;i++)
                h[lcl_idx0[i]]-=2.0*lcl_mtrx[i]*x[lcl_idx1[i]];
            for(int i=0;i<lsize;i++)
                h[i]+=2*c[i];
            for(int i=0;i<lsize;i++)
                g[comp][i]=h[i];
            
            g1_sq=g1_sq_l=0.0;
            for(int i=0;i<lsize;i++)
                g1_sq_l+=h[i]*h[i];
            
            MPI_Allreduce(&g1_sq_l,&g1_sq,1,MPI_DOUBLE,MPI_SUM,world);
        }

        xchng(h);
        
        hBh=hBh_l=0.0;
        for(int i=0;i<lcl_no_elem;i++)
            hBh_l+=h[lcl_idx0[i]]*lcl_mtrx[i]*h[lcl_idx1[i]];
        MPI_Allreduce(&hBh_l,&hBh,1,MPI_DOUBLE,MPI_SUM,world);
        if(hBh<=0)
        {
            printf("The matrix is singular");
            return 1;
        }
        
        hg_l=hg=0.0;
        for(int i=0;i<lsize;i++)
            hg_l+=g[comp][i]*h[i];
        MPI_Allreduce(&hg_l,&hg,1,MPI_DOUBLE,MPI_SUM,world);
        
        alpha=0.5*hg/hBh;
        if(alpha<=0)
        {
            printf("The matrix is singular (alpha became negative)");
            return 1;
        }
        
        for(int i=0;i<lcl_no_elem;i++)
            x[i]+=alpha*h[i];
        
        xchng(x);
        
        xBx=xBx_l=0.0;
        for(int i=0;i<lcl_no_elem;i++)
            xBx_l+=x[lcl_idx0[i]]*lcl_mtrx[i]*x[lcl_idx1[i]];
        cx_l=cx=0.0;
        for(int i=0;i<lsize;i++)
            cx_l+=x[i]*c[i];
        
        MPI_Allreduce(&xBx_l,&xBx,1,MPI_DOUBLE,MPI_SUM,world);
        MPI_Allreduce(&cx_l,&cx,1,MPI_DOUBLE,MPI_SUM,world);
        
        error=abs(xBx-2.0*cx+d_sq);
        if(error<=tol)
        {
            if(atoms->my_p_no==0)
                printf("Converged in %d steps!\n",iter);
            construct_ans();
            return 0;
        }
        g0_sq=g1_sq;
    }
    printf("Did not converge error: %lf!\n",error);
    construct_ans();
    return 1;
}
/*--------------------------------------------
 packing
 --------------------------------------------*/
void SOLVEAXb::construct_ans()
{
    int size=A->h1;
    double* ans_l;
    CREATE1D(ans_l,size);
    for(int i=0;i<size;i++)
        ans_l[i]=ans[i]=0.0;

    for(int i=lcl_lo;i<lcl_hi;i++)
        ans_l[i]=x[i-lcl_lo];
    
    MPI_Allreduce(&ans_l[0],&ans[0],size,MPI_DOUBLE,MPI_SUM,world);
    delete [] ans_l;
}
/*--------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 ---------------------------------------------
 --------------------------------------------*/
/*--------------------------------------------
 constructor
 --------------------------------------------*/
SPARSE_P::
SPARSE_P(MAPP* mapp,int n0,int n1): InitPtrs(mapp)
{
    h0=n0;
    h1=n1;
    no_elem=lcl_no_elem=0;
    myno=atoms->my_p_no;
    totp=atoms->tot_p;
    lcl_h0=static_cast<int>((1.0/totp)*h0+0.5);
    lo_h0=myno*lcl_h0;
    if(myno==totp-1)
        lcl_h0=h0-lcl_h0*(totp-1);
    hi_h0=lcl_h0+lcl_h0;
    
    lcl_no_vec=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
SPARSE_P::~SPARSE_P()
{

    if (lcl_no_elem)
    {
        delete [] lcl_idx0;
        delete [] lcl_idx1;
        delete [] lcl_mtrx;
    }
    if (lcl_no_vec)
    {
        delete [] lcl_idx_vec;
        delete [] lcl_init_pos;
        delete [] lcl_fin_pos;
    }
}
/*--------------------------------------------
 add elements
 --------------------------------------------*/
void SPARSE_P::add(int i0,int i1,double elem)
{
    if(i0>=h0 || i1>=h1 || i0<0 || i1<0)
        error->abort("in sparse the index "
                     "is higher than the number "
                     "columns/rows");
    
    if (elem==0.0)
        return;
    
    if(i0<lo_h0 || i0>=hi_h0)
        return;
    
    int lcl_i0=i0-lo_h0;
    
    int no=0;
    while (no<lcl_no_elem && (lcl_idx0[no]!=lcl_i0
                          || lcl_idx1[no]!=i1))
        no++;
    
    if(no<lcl_no_elem)
    {
        lcl_mtrx[no]=elem;
        return;
    }
    
    GROW(lcl_idx0,lcl_no_elem,no_elem+1);
    GROW(lcl_idx1,lcl_no_elem,no_elem+1);
    GROW(lcl_mtrx,lcl_no_elem,no_elem+1);
    
    lcl_idx0[no_elem]=lcl_i0;
    lcl_idx1[no_elem]=i1;
    lcl_mtrx[no_elem]=elem;
    lcl_no_elem++;
    no_elem++;
    
    MPI_Bcast(&no_elem,1,MPI_INT,0,world);
}
/*--------------------------------------------
 add elements
 --------------------------------------------*/
void SPARSE_P::addd(int i0,int i1,double elem)
{
    if(i0>=h0 || i1>=h1 || i0<0 || i1<0)
        error->abort("in sparse the index "
                     "is higher than the number "
                     "columns/rows");
    
    if (elem==0.0)
        return;
    
    if(i0<lo_h0 || i0>=hi_h0)
        return;
    
    int lcl_i0=i0-lo_h0;
    
    int no=0;
    while (no<lcl_no_elem && (lcl_idx0[no]!=lcl_i0
                              || lcl_idx1[no]!=i1))
        no++;
    
    if(no<lcl_no_elem)
    {
        lcl_mtrx[no]+=elem;
        return;
    }
    
    GROW(lcl_idx0,lcl_no_elem,no_elem+1);
    GROW(lcl_idx1,lcl_no_elem,no_elem+1);
    GROW(lcl_mtrx,lcl_no_elem,no_elem+1);
    
    lcl_idx0[no_elem]=lcl_i0;
    lcl_idx1[no_elem]=i1;
    lcl_mtrx[no_elem]=elem;
    lcl_no_elem++;
    no_elem++;
    
    MPI_Bcast(&no_elem,1,MPI_INT,0,world);
}
/*--------------------------------------------
 add elements
 --------------------------------------------*/
void SPARSE_P::add(int size,int* i0_lst,int* i1_lst,double* elem_lst)
{
    
    double elem;
    int i0,i1,lcl_i0;
    for(int i=0; i<size;i++)
    {
        elem=elem_lst[i];
        i0=i0_lst[i];
        i1=i1_lst[i];
        
        if (elem==0.0)
            continue;
        
        if(i0<lo_h0 || i0>=hi_h0)
            continue;
        
        lcl_i0=i0-lo_h0;
        
        int no=0;
        while (no<lcl_no_elem && (lcl_idx0[no]!=lcl_i0
                                  || lcl_idx1[no]!=i1))
            no++;
        
        if(no<lcl_no_elem)
        {
            lcl_mtrx[no]=elem;
            continue;
        }
        
        GROW(lcl_idx0,lcl_no_elem,no_elem+1);
        GROW(lcl_idx1,lcl_no_elem,no_elem+1);
        GROW(lcl_mtrx,lcl_no_elem,no_elem+1);
        
        lcl_idx0[no_elem]=lcl_i0;
        lcl_idx1[no_elem]=i1;
        lcl_mtrx[no_elem]=elem;
        lcl_no_elem++;
    }

    MPI_Allreduce(&lcl_no_elem,&no_elem,1,MPI_INT,MPI_SUM,world);
}
/*--------------------------------------------
 sort the sparse matrix row/column wise
 --------------------------------------------*/
void SPARSE_P::sort(int r_c)
{
    int tmp0,tmp1;
    double tmpx;
    if(r_c==0)
    {
        for(int i=0;i<lcl_no_elem;i++)
        {
            for(int j=i+1;j<lcl_no_elem;j++)
            {
                if(lcl_idx0[j]<lcl_idx0[i])
                {
                    tmp0=lcl_idx0[j];
                    tmp1=lcl_idx1[j];
                    tmpx=lcl_mtrx[j];
                    
                    lcl_idx0[j]=lcl_idx0[i];
                    lcl_idx1[j]=lcl_idx1[i];
                    lcl_mtrx[j]=lcl_mtrx[i];
                    
                    lcl_idx0[i]=tmp0;
                    lcl_idx1[i]=tmp1;
                    lcl_mtrx[i]=tmpx;
                }
            }
        }
    }
    else if(r_c==1)
    {
        for(int i=0;i<lcl_no_elem;i++)
        {
            for(int j=i+1;j<lcl_no_elem;j++)
            {
                if(lcl_idx1[j]<lcl_idx1[i])
                {
                    tmp0=lcl_idx0[j];
                    tmp1=lcl_idx1[j];
                    tmpx=lcl_mtrx[j];
                    
                    lcl_idx0[j]=lcl_idx0[i];
                    lcl_idx1[j]=lcl_idx1[i];
                    lcl_mtrx[j]=lcl_mtrx[i];
                    
                    lcl_idx0[i]=tmp0;
                    lcl_idx1[i]=tmp1;
                    lcl_mtrx[i]=tmpx;
                }
            }
        }
    }
    
    int strt=0;
    int end;
    int pos=strt;
    
    if (r_c==0)
    {
        while(pos<lcl_no_elem)
        {
            while(lcl_idx0[pos]==lcl_idx0[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            end=pos;
            for (int j=strt;j<end;j++)
            {
                for (int k=j+1;k<end;k++)
                {
                    if(lcl_idx1[j]>lcl_idx1[k])
                    {
                        tmp1=lcl_idx1[k];
                        tmpx=lcl_mtrx[k];
                        
                        lcl_idx1[k]=lcl_idx1[j];
                        lcl_mtrx[k]=lcl_mtrx[j];
                        
                        lcl_idx1[j]=tmp1;
                        lcl_mtrx[j]=tmpx;
                    }
                }
            }
            strt=pos;
        }
    }
    else if (r_c==1)
    {
        while(pos<no_elem)
        {
            while(lcl_idx1[pos]==lcl_idx1[strt])
            {
                pos++;
                if(pos==no_elem)
                    break;
            }
            end=pos;
            for (int j=strt;j<end;j++)
            {
                for (int k=j+1;k<end;k++)
                {
                    if(lcl_idx0[j]>lcl_idx0[k])
                    {
                        tmp1=lcl_idx0[k];
                        tmpx=lcl_mtrx[k];
                        
                        lcl_idx0[k]=lcl_idx0[j];
                        lcl_mtrx[k]=lcl_mtrx[j];
                        
                        lcl_idx0[j]=tmp1;
                        lcl_mtrx[j]=tmpx;
                    }
                }
            }
            strt=pos;
        }
    }
}
/*--------------------------------------------
 vectorize the matrix row/column wise
 --------------------------------------------*/
void SPARSE_P::vectorize(int r_c)
{
    sort(r_c);
    
    if(lcl_no_vec)
    {
        delete [] lcl_fin_pos;
        delete [] lcl_idx_vec;
        delete [] lcl_init_pos;
    }
    
    
    lcl_no_vec=0;
    int strt=0;
    int pos=strt;
    if (r_c==0)
    {
        while(pos<lcl_no_elem)
        {
            while(lcl_idx0[pos]==lcl_idx0[strt])
            {
                pos++;
                if(pos==lcl_no_elem)
                    break;
            }
            GROW(lcl_init_pos,lcl_no_vec,lcl_no_vec+1);
            GROW(lcl_fin_pos,lcl_no_vec,lcl_no_vec+1);
            GROW(lcl_idx_vec,lcl_no_vec,lcl_no_vec+1);
            
            lcl_init_pos[lcl_no_vec]=strt;
            lcl_fin_pos[lcl_no_vec]=pos;
            lcl_idx_vec[lcl_no_vec]=lcl_idx0[strt];
            
            lcl_no_vec++;
            strt=pos;
        }
    }
    else if (r_c==1)
    {
        while(pos<lcl_no_elem)
        {
            while(lcl_idx1[pos]==lcl_idx1[strt])
            {
                pos++;
                if(pos==lcl_no_elem)
                    break;
            }
            
            GROW(lcl_init_pos,lcl_no_vec,lcl_no_vec+1);
            GROW(lcl_fin_pos,lcl_no_vec,lcl_no_vec+1);
            GROW(lcl_idx_vec,lcl_no_vec,lcl_no_vec+1);
            
            lcl_init_pos[lcl_no_vec]=strt;
            lcl_fin_pos[lcl_no_vec]=pos;
            lcl_idx_vec[lcl_no_vec]=lcl_idx1[strt];
            
            lcl_no_vec++;
            strt=pos;
        }
    }
    
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
COMB::COMB(MAPP* mapp):InitPtrs(mapp)
{
    perm_list_size_0=0;
    perm_list_size_1=0;
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
COMB::~COMB()
{
    clean();
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
void COMB::clean()
{
    if(perm_list_size_0)
    {
        for(int i=0;i<perm_list_size_1;i++)
            delete [] perm_list[i];
    }
    if(perm_list_size_1)
        delete [] perm_list;
    perm_list_size_0=0;
    perm_list_size_1=0;
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
void COMB::comb(int* no_list, int list_size)
{
    int tot=0;
    for(int i=0;i<list_size;i++)
        tot+=no_list[i];

    clean();
    
    perm_list_size_0=tot;
    perm_list_size_1=0;
    
    int* tmp0;
    CREATE1D(tmp0,tot);
    int* tmp1;
    CREATE1D(tmp1,tot);
    
    comb_rec(no_list[0],tot,0,no_list[0]
    ,tmp0,tmp1,no_list,0,tot);
    delete [] tmp0;
    delete [] tmp1;
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
void COMB::comb_rec(int no,int tot,int pos
,int max_no,int*& lvl,int*& tmp,int* y
,int y_pos,int y_tot)
{
    if (no!=0)
    {
        for(int i=pos;i<tot;i++)
        {
            lvl[no-1]=i;
            comb_rec(no-1,tot,i+1,max_no,lvl,tmp,y,y_pos,y_tot);
            
        }
    }
    else
    {
        if(y_pos==0)
            for (int i=0;i<tot;i++)
                tmp[i]=-1;
        int p=0;
        int pp=0;
        int ppp=max_no-1;
        while(ppp>=0)
        {
            if(tmp[p]==-1)
            {
                if(pp==lvl[ppp])
                {
                    tmp[p]=y_pos;
                    ppp--;
                }
                pp++;
            }
            p++;
        }
        
        if(y[y_pos]==tot)
        {

            GROW(perm_list,perm_list_size_1,perm_list_size_1+1);
            CREATE1D(perm_list[perm_list_size_1], perm_list_size_0);
            perm_list_size_1++;
            
            for(int i=0;i<y_tot;i++)
                if(tmp[i]==y_pos)
                    tmp[i]=-1;
        }
        else
        {
            int noo=y[y_pos+1];
            int tott=tot-y[y_pos];
            int* h=lvl+y[y_pos];
            comb_rec(noo,tott,0,noo,h,tmp,y,y_pos+1,y_tot);
            for(int i=0;i<y_tot;i++)
                if(tmp[i]==y_pos)
                    tmp[i]=-1;
        }
    }
}

