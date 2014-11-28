#include "neighbor.h"
#include "atoms.h"
#include <cstdlib>
using namespace std;
using namespace MAPP_NS;

/*--------------------------------------------
 
 --------------------------------------------*/
Atoms::Atoms(MAPP* mapp,MPI_Comm communicator)
:InitPtrs(mapp)
{
    comm_world=communicator;
    dimension=3;
    natms=0;
    natms_ph=0;
    atm_vec_ph_size=0;
    atm_vec_size=0;
    tot_cut_ph=0.0;
    skin=0.5;
    
    CREATE2D(H,dimension,dimension);
    CREATE2D(B,dimension,dimension);
    
    CREATE1D(s_lo,dimension);
    CREATE1D(s_hi,dimension);
    CREATE1D(s_ph_lo,dimension);
    CREATE1D(s_ph_hi,dimension);
    CREATE1D(cut_ph_s,dimension);
    CREATE1D(comm_need,dimension);
    
    CREATE1D(tot_p_grid,dimension);
    CREATE1D(my_loc,dimension);
    CREATE2D(neigh_p,dimension,2);
    
    
    //communication related parameters
    MPI_Comm_rank(comm_world,&my_p_no);
    MPI_Comm_size(comm_world,&tot_p);
    ph_lst= new SwapLst<int>(mapp);
    
    snd_buff_0_capacity=0;
    snd_buff_1_capacity=0;
    rcv_buff_capacity=0;
    snd_ph_buff_capacity=0;
    rcv_ph_buff_capacity=0;
    
    no_vecs=0;
    tot_byte_size=0;
    
    
    int name_length;
    char* node_name;
    CREATE1D(node_name,MPI_MAX_PROCESSOR_NAME);
    MPI_Get_processor_name(node_name,&name_length);
    node_name[name_length] = '\0';
    name_length++;
    

    int* name_lenghts;
    CREATE1D(name_lenghts,tot_p);
    for(int i=0;i<tot_p;i++)
        name_lenghts[i]=0;
    
    
    name_lenghts[my_p_no]=name_length;
    int* all_name_lenghts;
    CREATE1D(all_name_lenghts,tot_p);
    MPI_Allreduce(name_lenghts,all_name_lenghts,tot_p, MPI_INT,MPI_SUM,comm_world);
    delete [] name_lenghts;
    
    
    
    
    char** all_names;
    CREATE1D(all_names,tot_p);
    for(int i=0;i<tot_p;i++)
        CREATE1D(all_names[i],all_name_lenghts[i]);
    
    
    
    for(int i=0;i<name_length;i++)
        all_names[my_p_no][i]=node_name[i];
    delete [] node_name;
    
    
    for(int i=0;i<tot_p;i++)
        MPI_Bcast(all_names[i],all_name_lenghts[i],MPI_CHAR,i,comm_world);
    
    
    int* node_no;
    CREATE1D(node_no,tot_p);
    
    for(int i=0;i<tot_p;i++)
        node_no[i]=-1;  
    
    tot_n=0;
    for(int i=0;i<tot_p;i++)
    {
        if(node_no[i]==-1)
        {
            node_no[i]=tot_n;
            tot_n++;
        }
        for(int j=i+1;j<tot_p;j++)
            if(strcmp(all_names[i],all_names[j])==0)
                node_no[j]=node_no[i];
    }
    
    
    for(int i=0;i<tot_p;i++)
        if(all_name_lenghts[i])
            delete [] all_names[i];
    delete [] all_names;
    
    
    CREATE1D(p_per_n,tot_n);
    
    for(int i=0;i<tot_n;i++)
        p_per_n[i]=0;
    
    for(int i=0;i<tot_p;i++)
        p_per_n[node_no[i]]++;
    
    my_n_no=node_no[my_p_no];

    
    CREATE1D(n_p_grid,tot_n);
    for(int i=0;i<tot_n;i++)
        CREATE1D(n_p_grid[i],p_per_n[i]);
    
    int pos;
    for(int inode=0;inode<tot_n;inode++)
    {
        pos=0;
        for(int iproc=0;iproc<tot_p;iproc++)
        {
            if(node_no[iproc]==inode)
            {
                n_p_grid[inode][pos]=iproc;
                pos++;
            }
        }
    }
    delete [] node_no;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
Atoms::~Atoms()
{

    if(dimension)
    {
        for(int i=0;i<dimension;i++)
        {
            delete [] H[i];
            delete [] B[i];
            delete [] neigh_p[i];
        }
        delete [] H;
        delete [] B;
        delete [] s_lo;
        delete [] s_hi;
        delete [] s_ph_lo;
        delete [] s_ph_hi;
        delete [] cut_ph_s;
        delete [] comm_need;
        
        delete [] tot_p_grid;
        delete [] my_loc;
        delete [] neigh_p;
    }
    
    delete ph_lst;
    
    
    if(snd_buff_0_capacity)
        delete [] snd_buff_0;
    if(snd_buff_1_capacity)
        delete [] snd_buff_1;
    if(rcv_buff_capacity)
        delete [] rcv_buff;
    if(snd_ph_buff_capacity)
        delete [] snd_ph_buff;
    if(rcv_ph_buff_capacity)
        delete [] rcv_ph_buff;
    
    if(no_vecs)
        delete [] vectors;
    
    for(int i=0;i<tot_n;i++)
        delete n_p_grid[i];
    delete [] n_p_grid;
    delete [] p_per_n;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::chng_dim(int dim)
{
    if(dim<1)
        error->abort("box dimension should be more than 0");
    
    if(dimension)
    {
        for(int i=0;i<dimension;i++)
        {
            delete [] H[i];
            delete [] B[i];
            delete [] neigh_p[i];
        }
        delete [] H;
        delete [] B;
        delete [] s_lo;
        delete [] s_hi;
        delete [] s_ph_lo;
        delete [] s_ph_hi;
        delete [] cut_ph_s;
        delete [] comm_need;
        
        delete [] tot_p_grid;
        delete [] my_loc;
        delete [] neigh_p;
    }
    
    
    dimension=dim;
    
    CREATE2D(H,dimension,dimension);
    CREATE2D(B,dimension,dimension);
    
    CREATE1D(s_lo,dimension);
    CREATE1D(s_hi,dimension);
    CREATE1D(s_ph_lo,dimension);
    CREATE1D(s_ph_hi,dimension);
    CREATE1D(cut_ph_s,dimension);
    CREATE1D(comm_need,dimension);
    
    CREATE1D(tot_p_grid,dimension);
    CREATE1D(my_loc,dimension);
    CREATE2D(neigh_p,dimension,2);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::set_ph(TYPE0 cut)
{
    tot_cut_ph=cut+skin;
    TYPE0 tmp;
    for (int i=0;i<dimension;i++)
    {
        tmp=0.0;
        for(int j=i;j<dimension;j++)
            tmp+=B[j][i]*B[j][i];
        
        cut_ph_s[i]=tot_cut_ph*sqrt(tmp);
    }
       
    for (int i=0;i<dimension;i++)
    {
        s_ph_hi[i]=s_hi[i]-cut_ph_s[i];
        s_ph_lo[i]=s_lo[i]+cut_ph_s[i];
    }
    for (int i=0;i<dimension;i++)
        comm_need[i]=static_cast<int>
        (cut_ph_s[i]/(s_hi[i]-s_lo[i]))+1;
    
    int nswap=0;
    for(int i=0;i<dimension;i++)
        nswap+=comm_need[i];
    nswap*=2;
    
    ph_lst->newlist(nswap);
    ph_lst->grow_size=SWAPGROWTH;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::set_ph(int box_change)
{
    if(box_change==0)
        return;
    
    TYPE0 tmp;
    for (int i=0;i<dimension;i++)
    {
        tmp=0.0;
        for(int j=i;j<dimension;j++)
            tmp+=B[j][i]*B[j][i];
        
        cut_ph_s[i]=tot_cut_ph*sqrt(tmp);
    }
    
    for (int i=0;i<dimension;i++)
    {
        s_ph_hi[i]=s_hi[i]-cut_ph_s[i];
        s_ph_lo[i]=s_lo[i]+cut_ph_s[i];
    }
    for (int i=0;i<dimension;i++)
        comm_need[i]=static_cast<int>
        (cut_ph_s[i]/(s_hi[i]-s_lo[i]))+1;
    
    int nswap=0;
    for(int i=0;i<dimension;i++)
        nswap+=comm_need[i];
    nswap*=2;
    
    ph_lst->newlist(nswap);
    ph_lst->grow_size=SWAPGROWTH;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::ph_setup(int box_chng,class VecLst* list)
{
    set_ph(box_chng);
    
    natms_ph=0;
    ph_lst->reset();
    int lo_r_dir,hi_r_dir,lo_r_dim,hi_r_dim,tmp_r;
    
    TYPE0* s;
    int s_dim=vectors[0].dim;
    int iswap=0;
    
    lo_r_dim=0;
    hi_r_dim=natms+natms_ph;
    for(int idim=0;idim<dimension;idim++)
    {
        lo_r_dir=lo_r_dim;
        hi_r_dir=hi_r_dim;
                
        for(int i=0;i<comm_need[idim];i++)
        {
            tmp_r=natms+natms_ph;
            //s=(TYPE0*)vectors[0].ret_vec();
            vectors[0].ret(s);
            for(int iatm=lo_r_dir;iatm<hi_r_dir;iatm++)
                if(s[iatm*s_dim+idim]<=s_ph_lo[idim])
                    ph_lst->add(iswap,iatm);
            
            ph_xchng(0,idim,ph_lst->list[iswap],ph_lst->pos[iswap],list);
            hi_r_dir=natms+natms_ph;
            lo_r_dir=tmp_r;
            iswap++;
        }

        lo_r_dir=lo_r_dim;
        hi_r_dir=hi_r_dim;
        
        for(int i=0;i<comm_need[idim];i++)
        {
            tmp_r=natms+natms_ph;
            //s=(TYPE0*)vectors[0].ret_vec();
            vectors[0].ret(s);
            for(int iatm=lo_r_dir;iatm<hi_r_dir;iatm++)
                if(s_ph_hi[idim]<=s[iatm*s_dim+idim])
                    ph_lst->add(iswap,iatm);
            ph_xchng(1,idim,ph_lst->list[iswap],ph_lst->pos[iswap],list);
            hi_r_dir=natms+natms_ph;
            lo_r_dir=tmp_r;
            iswap++;
        }
        
        lo_r_dim=0;
        hi_r_dim=natms+natms_ph;
    }

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::ph_xchng(int dir,int dim,
int* atm_list,int atm_list_size,class VecLst* list)
{
    
    if(tot_p_grid[dim]!=1)
    {


        int snd_proc=neigh_p[dim][dir];
        int rcv_proc=neigh_p[dim][1-dir];
        int snd_buff_size=list->ph_byte_size
        *atm_list_size;
        int rcv_buff_size;
        if(snd_buff_size>snd_ph_buff_capacity)
        {
            if(snd_ph_buff_capacity)
                delete [] snd_ph_buff;
            CREATE1D(snd_ph_buff,snd_buff_size);
            snd_ph_buff_capacity=snd_buff_size;
        }
        
        MPI_Request request[2];
        MPI_Status status[2];
        MPI_Sendrecv(&snd_buff_size,1,
            MPI_INT,snd_proc,0,&rcv_buff_size,
            1,MPI_INT,rcv_proc,0,comm_world,
            &status[0]);

        
        if(rcv_buff_size>rcv_ph_buff_capacity)
        {
            if(rcv_ph_buff_capacity)
                delete [] rcv_ph_buff;
            CREATE1D(rcv_ph_buff,rcv_buff_size);
            rcv_ph_buff_capacity=rcv_buff_size;
        }
        x_pack(snd_ph_buff,list,atm_list,atm_list_size);
 
        
        if (rcv_buff_size)
        {
            MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                      MPI_BYTE,rcv_proc,0,comm_world,
                      &request[0]);
        }
        if (snd_buff_size)
        {
            MPI_Isend(snd_ph_buff,snd_buff_size,
                      MPI_BYTE,snd_proc,0,comm_world,
                      &request[1]);
        }

        int no_new_atms=rcv_buff_size/list->ph_byte_size;
        
        if (rcv_buff_size&&snd_buff_size==0)
        {
            MPI_Wait(&request[0],&status[0]);
            x_unpack(rcv_ph_buff,list, no_new_atms);
        }
        else if (rcv_buff_size==0&&snd_buff_size)
        {
            MPI_Wait(&request[1],&status[1]);
        }
        else if (rcv_buff_size&&snd_buff_size)
        {
            MPI_Waitall(2,request,status);
            x_unpack(rcv_ph_buff,list, no_new_atms);
        }
        
        //correction due to periodic boundary conditions
        
        if((dir==0&&my_loc[dim]==tot_p_grid[dim]-1) ||
           (dir==1&&my_loc[dim]==0))
        {
            int last_atm_no=natms+natms_ph;
            int strt_atm_no=last_atm_no-no_new_atms;
            //TYPE0* s=(TYPE0*)vectors[0].ret_vec();
            TYPE0* s;
            vectors[0].ret(s);
            int s_dim=vectors[0].dim;
            if(dir==0)
                for (int i=strt_atm_no;i<last_atm_no;i++)
                    s[i*s_dim+dim]++;
            else if(dir==1)
                for (int i=strt_atm_no;i<last_atm_no;i++)
                    s[i*s_dim+dim]--;
        }
    }
    else
    {
        int new_size=natms_ph+natms
        +atm_list_size;
        int old_size=atm_vec_ph_size;
        
        if(new_size>old_size)
        {
            grow(1,new_size-old_size,list);
        }
        
        int last_atm_no=natms+natms_ph;
        
        for(int ivec=0;ivec<list->ph_no_vecs;ivec++)
            for(int i=0;i<atm_list_size;i++)
                vectors[list->ph_vec_list[ivec]]
                .copy(last_atm_no+i,atm_list[i]);
    
        
        //correction due to periodic boundary conditions
        //TYPE0* s=(TYPE0*)vectors[0].ret_vec();
        TYPE0* s;
        vectors[0].ret(s);
        int s_dim=vectors[0].dim;
        
        if(dir==0)
            for (int i=natms+natms_ph;i<natms+natms_ph+atm_list_size;i++)
                s[i*s_dim+dim]++;
        else if(dir==1)
            for (int i=natms+natms_ph;i<natms+natms_ph+atm_list_size;i++)
                s[i*s_dim+dim]--;
        
        natms_ph+=atm_list_size;
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::reset_comm(class VecLst* list)
{
    for(int idim=0;idim<dimension;idim++)
        reset_comm(idim,list);

}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::reset_comm(int idim,class VecLst* list)
{
    if(tot_p_grid[idim]==1)
        return;
    int icomp;

    TYPE0* s;
    vectors[0].ret(s);
    int s_dim=vectors[0].dim;
    
    int snd_buff_0_size=0;
    int snd_buff_1_size=0;
    
    int iatm=0;
    
    while (iatm<natms)
    {
        icomp=s_dim*iatm;
        if(s[icomp+idim]<s_lo[idim])
        {
            if(my_loc[idim]!=tot_p_grid[idim]-1)
                snd_buff_0_size=pack(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            else if(my_loc[idim]==tot_p_grid[idim]-1&&0.5<=s[icomp+idim])
                snd_buff_0_size=pack(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            else if(my_loc[idim]==tot_p_grid[idim]-1&&s[icomp+idim]<0.5)
                snd_buff_1_size=pack(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
        }
        else if(s[icomp+idim]>=s_hi[idim])
        {
            if(my_loc[idim]!=0)
                snd_buff_1_size=pack(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
            else if (my_loc[idim]==0&&s[icomp+idim]<0.5)
            {
                snd_buff_1_size=pack(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
            }
            else if (my_loc[idim]==0&&0.5<=s[icomp+idim])
            {
                snd_buff_0_size=pack(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            }
        }
        else iatm++;
    }
    
    int rcv_buff_size;
    MPI_Request request;
    MPI_Status status;
    
    int no_atms;
    int d_size=sizeof(TYPE0);
    
    int max_snd_buff_0_size;
    MPI_Allreduce(&snd_buff_0_size,&max_snd_buff_0_size,1,
                  MPI_INT,MPI_MAX,comm_world);
    while(max_snd_buff_0_size)
    {
        int rcv_p=neigh_p[idim][1];
        int snd_p=neigh_p[idim][0];
        MPI_Sendrecv(&snd_buff_0_size,1,
                     MPI_INT,snd_p,0,&rcv_buff_size,
                     1,MPI_INT,rcv_p,0,comm_world,
                     &status);
        
        if(rcv_buff_capacity<rcv_buff_size)
        {
            if(rcv_buff_capacity)
                delete [] rcv_buff;
            CREATE1D(rcv_buff,rcv_buff_size);
            rcv_buff_capacity=rcv_buff_size;
        }
        if(rcv_buff_size)
            MPI_Irecv(rcv_buff,rcv_buff_size,
                      MPI_BYTE,rcv_p,0,comm_world,
                      &request);
        if(snd_buff_0_size)
            MPI_Send(snd_buff_0,snd_buff_0_size,
                     MPI_BYTE,snd_p,0,comm_world);
        
        snd_buff_0_size=0;
        if (rcv_buff_size)
        {
            MPI_Wait(&request,&status);
            
            no_atms=rcv_buff_size/list->byte_size;
            
            for(int i=0;i<no_atms;i++)
            {
                TYPE0 s_i;
                memcpy(&s_i,&rcv_buff[i*list->byte_size+idim*d_size],d_size);
                if(s_i<s_lo[idim])
                {
                    if(snd_buff_0_capacity<snd_buff_0_size+list->byte_size)
                    {
                        GROW(snd_buff_0,snd_buff_0_capacity,snd_buff_0_size+list->byte_size);
                        snd_buff_0_capacity=snd_buff_0_size+list->byte_size;
                    }
                    memcpy(&(snd_buff_0[snd_buff_0_size]),&(rcv_buff[i*list->byte_size]),list->byte_size);
                    snd_buff_0_size+=list->byte_size;
                }
                else
                    unpack(rcv_buff,i*list->byte_size,1,list);
            }
        }
        MPI_Allreduce(&snd_buff_0_size,&max_snd_buff_0_size,1,
                      MPI_INT,MPI_MAX,comm_world);
    }
    
    int max_snd_buff_1_size;
    MPI_Allreduce(&snd_buff_1_size,&max_snd_buff_1_size,1,
                  MPI_INT,MPI_MAX,comm_world);
    while(max_snd_buff_1_size)
    {
        int rcv_p=neigh_p[idim][0];
        int snd_p=neigh_p[idim][1];
        MPI_Sendrecv(&snd_buff_1_size,1,
                     MPI_INT,snd_p,0,&rcv_buff_size,
                     1,MPI_INT,rcv_p,0,comm_world,
                     &status);
        
        if(rcv_buff_capacity<rcv_buff_size)
        {
            if(rcv_buff_capacity)
                delete [] rcv_buff;
            CREATE1D(rcv_buff,rcv_buff_size);
            rcv_buff_capacity=rcv_buff_size;
        }
        if(rcv_buff_size)
            MPI_Irecv(rcv_buff,rcv_buff_size,
                      MPI_BYTE,rcv_p,0,comm_world,
                      &request);
        if(snd_buff_1_size)
            MPI_Send(snd_buff_1,snd_buff_1_size,
                     MPI_BYTE,snd_p,0,comm_world);
        
        snd_buff_1_size=0;
        if (rcv_buff_size)
        {
            MPI_Wait(&request,&status);
            
            no_atms=rcv_buff_size/list->byte_size;
            for(int i=0;i<no_atms;i++)
            {
                TYPE0 s_i;
                memcpy(&s_i,&rcv_buff[i*list->byte_size+idim*d_size],d_size);
                if(s_hi[idim]<=s_i)
                {
                    if(snd_buff_1_capacity<snd_buff_1_size+list->byte_size)
                    {
                        GROW(snd_buff_1,snd_buff_1_capacity,snd_buff_1_size+list->byte_size);
                        snd_buff_1_capacity=snd_buff_1_size+list->byte_size;
                    }
                    memcpy(&snd_buff_1[snd_buff_1_size],&rcv_buff[i*list->byte_size],list->byte_size);
                    snd_buff_1_size+=list->byte_size;
                }
                else
                    unpack(rcv_buff,i*list->byte_size,1,list);
            }
        }
        MPI_Allreduce(&snd_buff_1_size,&max_snd_buff_1_size,1,
                      MPI_INT,MPI_MAX,comm_world);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::xchng_comm(class VecLst* list)
{
    for(int idim=0;idim<dimension;idim++)
        xchng_comm(idim,list);


}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::xchng_comm(int idim,class VecLst* list)
{
    if(tot_p_grid[idim]==1)
        return;
    int icomp;
    //TYPE0* s=(TYPE0*)vectors[0].ret_vec();
    TYPE0* s;
    vectors[0].ret(s);
    int s_dim=vectors[0].dim;
    
    int snd_buff_0_size=0;
    int snd_buff_1_size=0;
    
    
    int iatm=0;
    while (iatm<natms)
    {
        icomp=s_dim*iatm;
        if(s[icomp+idim]<s_lo[idim])
        {
            
            if(my_loc[idim]!=tot_p_grid[idim]-1)
                snd_buff_0_size=pack(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            else if(my_loc[idim]==tot_p_grid[idim]-1&&0.5<=s[icomp+idim])
                snd_buff_0_size=pack(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            else if(my_loc[idim]==tot_p_grid[idim]-1&&s[icomp+idim]<0.5)
                snd_buff_1_size=pack(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
        }
        else if(s[icomp+idim]>=s_hi[idim])
        {
            if(my_loc[idim]!=0)
                snd_buff_1_size=pack(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
            else if (my_loc[idim]==0&&s[icomp+idim]<0.5)
            {
                snd_buff_1_size=pack(snd_buff_1,snd_buff_1_size,snd_buff_1_capacity,iatm,list);
            }
            else if (my_loc[idim]==0&&0.5<=s[icomp+idim])
            {
                snd_buff_0_size=pack(snd_buff_0,snd_buff_0_size,snd_buff_0_capacity,iatm,list);
            }
        }
        else iatm++;
    }
    int rcv_buff_size,rcv_p,snd_p;
    MPI_Request request;
    MPI_Status status;
    int no_atms;    

    rcv_p=neigh_p[idim][1];
    snd_p=neigh_p[idim][0];
    MPI_Sendrecv(&snd_buff_0_size,1,
                 MPI_INT,snd_p,0,&rcv_buff_size,
                 1,MPI_INT,rcv_p,0,comm_world,
                 &status);
    if(rcv_buff_capacity<rcv_buff_size)
    {
        if(rcv_buff_capacity)
            delete [] rcv_buff;
        CREATE1D(rcv_buff,rcv_buff_size);
        rcv_buff_capacity=rcv_buff_size;
    }
    
    if(rcv_buff_size)
        MPI_Irecv(rcv_buff,rcv_buff_size,
                  MPI_BYTE,rcv_p,0,comm_world,
                  &request);
    
    if(snd_buff_0_size)
        MPI_Send(snd_buff_0,snd_buff_0_size,
                 MPI_BYTE,snd_p,0,comm_world);
    snd_buff_0_size=0;
    if (rcv_buff_size)
    {
        
        MPI_Wait(&request,&status);
        no_atms=rcv_buff_size/list->byte_size;
        
        unpack(rcv_buff,0,no_atms,list);
        
        //s=(TYPE0*)vectors[0].ret_vec();
        vectors[0].ret(s);
        icomp=(natms-no_atms)*s_dim+idim;
        for(int i=0;i<no_atms;i++)
        {
            if(s[icomp]<s_lo[idim])
                error->abort("atom moved more than one processor");
            icomp+=s_dim;
        }
        
    }
    
    
    rcv_p=neigh_p[idim][0];
    snd_p=neigh_p[idim][1];
    MPI_Sendrecv(&snd_buff_1_size,1,
                 MPI_INT,snd_p,0,&rcv_buff_size,
                 1,MPI_INT,rcv_p,0,comm_world,
                 &status);
    if(rcv_buff_capacity<rcv_buff_size)
    {
        if(rcv_buff_capacity)
            delete [] rcv_buff;
        CREATE1D(rcv_buff,rcv_buff_size);
        rcv_buff_capacity=rcv_buff_size;
    }
    if(rcv_buff_size)
        MPI_Irecv(rcv_buff,rcv_buff_size,
                  MPI_BYTE,rcv_p,0,comm_world,
                  &request);
    if(snd_buff_1_size)
        MPI_Send(snd_buff_1,snd_buff_1_size,
                 MPI_BYTE,snd_p,0,comm_world);
    snd_buff_1_size=0;
    if (rcv_buff_size)
    {
        MPI_Wait(&request,&status);
        
        no_atms=rcv_buff_size/list->byte_size;
        
        unpack(rcv_buff,0,no_atms,list);
        
        //s=(TYPE0*)vectors[0].ret_vec();
        vectors[0].ret(s);
        icomp=(natms-no_atms)*s_dim+idim;
        for(int i=0;i<no_atms;i++)
        {
            if(s_hi[idim]<=s[icomp])
                error->abort("atom moved more than one processor");
            icomp+=s_dim;
        }

    }
 
}
/*--------------------------------------------
 update number of vectors
 --------------------------------------------*/
void Atoms::update(class VecLst* list)
{
    int iswap=0;
    natms_ph=0;
    
    if(list->ph_no_vecs==1)
    {
        for(int idim=0;idim<dimension;idim++)
        {
            if(tot_p_grid[idim]!=1)
            {
                for(int idir=0;idir<2;idir++)
                {
                    for(int icomm=0;icomm<comm_need[idim];icomm++)
                    {
                        int snd_proc=neigh_p[idim][idir];
                        int rcv_proc=neigh_p[idim][1-idir];
                        int snd_buff_size=list->ph_byte_size
                        *ph_lst->pos[iswap];
                        int rcv_buff_size;
                        if(snd_buff_size>snd_ph_buff_capacity)
                        {
                            if(snd_ph_buff_capacity)
                                delete [] snd_ph_buff;
                            CREATE1D(snd_ph_buff,snd_buff_size);
                            snd_ph_buff_capacity=snd_buff_size;
                        }
                        
                        x_pack(snd_ph_buff,list,ph_lst->list[iswap],ph_lst->pos[iswap]);
                        
                        MPI_Request request[2];
                        MPI_Status status[2];
                        
                        MPI_Sendrecv(&snd_buff_size,1,MPI_INT,snd_proc,0,&rcv_buff_size,
                        1,MPI_INT,rcv_proc,0,comm_world,&status[0]);
                        
                        if (rcv_buff_size)
                        {
                            MPI_Irecv(vectors[list->ph_vec_list[0]].ret_vec(natms+natms_ph)
                            ,rcv_buff_size,MPI_BYTE,rcv_proc,0,comm_world,&request[0]);
                            
                        }
                        if (snd_buff_size)
                        {
                            MPI_Isend(snd_ph_buff,snd_buff_size,
                            MPI_BYTE,snd_proc,0,comm_world,&request[1]);
                        }
                        
                        int no_new_atms=rcv_buff_size/list->ph_byte_size;
                        
                        if (rcv_buff_size&&snd_buff_size==0)
                            MPI_Wait(&request[0],&status[0]);
                        else if (rcv_buff_size==0&&snd_buff_size)
                            MPI_Wait(&request[1],&status[1]);
                        else if (rcv_buff_size&&snd_buff_size)
                            MPI_Waitall(2,request,status);
                        
                        natms_ph+=no_new_atms;
                        if(list->ph_vec_list[0]==0)
                        {
                            if((idir==0&&my_loc[idim]==tot_p_grid[idim]-1) ||
                               (idir==1&&my_loc[idim]==0))
                            {
                                int last_atm_no=natms+natms_ph;
                                int strt_atm_no=last_atm_no-no_new_atms;
                                //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
                                TYPE0* x;
                                vectors[0].ret(x);
                                int x_dim=vectors[0].dim;
                                
                                if(idir==0)
                                    for (int i=strt_atm_no;i<last_atm_no;i++)
                                        for(int j=0;j<=idim;j++)
                                            x[i*x_dim+j]+=H[idim][j];
                                else if(idir==1)
                                    for (int i=strt_atm_no;i<last_atm_no;i++)
                                        for(int j=0;j<=idim;j++)
                                            x[i*x_dim+j]-=H[idim][j];
                            }
                        }
                        iswap++;
                    }
                }
            }
            else
            {
                for(int idir=0;idir<2;idir++)
                {
                    for(int icomm=0;icomm<comm_need[idim];icomm++)
                    {
                        int last_atm_no=natms+natms_ph;
                        
                        for(int ivec=0;ivec<list->ph_no_vecs;ivec++)
                            for(int i=0;i<ph_lst->pos[iswap];i++)
                                vectors[list->ph_vec_list[ivec]]
                                .copy(last_atm_no+i,ph_lst->list[iswap][i]);
                        
                        if(list->ph_vec_list[0]==0)
                        {
                            //correction due to periodic boundary conditions
                            //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
                            TYPE0* x;
                            vectors[0].ret(x);
                            int x_dim=vectors[0].dim;
                            
                            if(idir==0)
                                for (int i=natms+natms_ph;i<natms+natms_ph+ph_lst->pos[iswap];i++)
                                    for(int j=0;j<=idim;j++)
                                        x[i*x_dim+j]+=H[idim][j];

                            else if(idir==1)
                                for (int i=natms+natms_ph;i<natms+natms_ph+ph_lst->pos[iswap];i++)
                                    for(int j=0;j<=idim;j++)
                                        x[i*x_dim+j]-=H[idim][j];
                            
                        }
                        natms_ph+=ph_lst->pos[iswap];
                        iswap++;
                    }
                }
            }
            
        }
    }
    else
    {
        for(int idim=0;idim<dimension;idim++)
        {
            if(tot_p_grid[idim]!=1)
            {
                for(int idir=0;idir<2;idir++)
                {
                    for(int icomm=0;icomm<comm_need[idir];icomm++)
                    {
                        int snd_proc=neigh_p[idim][idir];
                        int rcv_proc=neigh_p[idim][1-idir];
                        int snd_buff_size=list->ph_byte_size
                        *ph_lst->pos[iswap];
                        int rcv_buff_size;
                        if(snd_buff_size>snd_ph_buff_capacity)
                        {
                            if(snd_ph_buff_capacity)
                                delete [] snd_ph_buff;
                            CREATE1D(snd_ph_buff,snd_buff_size);
                            snd_ph_buff_capacity=snd_buff_size;
                        }
                        
                        x_pack(snd_ph_buff,list,ph_lst->list[iswap],ph_lst->pos[iswap]);
                        
                        MPI_Request request[2];
                        MPI_Status status[2];
                        
                        MPI_Sendrecv(&snd_buff_size,1,
                        MPI_INT,snd_proc,0,&rcv_buff_size,
                        1,MPI_INT,rcv_proc,0,comm_world,&status[0]);
                        
                        if(rcv_buff_size>rcv_ph_buff_capacity)
                        {
                            if(rcv_ph_buff_capacity)
                                delete [] rcv_ph_buff;
                            CREATE1D(rcv_ph_buff,rcv_buff_size);
                            rcv_ph_buff_capacity=rcv_buff_size;
                        }
                        
                        if (rcv_buff_size)
                        {
                            MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                                      MPI_BYTE,rcv_proc,0,comm_world,
                                      &request[0]);
                            
                        }
                        if (snd_buff_size)
                        {
                            MPI_Isend(snd_ph_buff,snd_buff_size,
                                      MPI_BYTE,snd_proc,0,comm_world,
                                      &request[1]);
                        }
                        
                        int no_new_atms=rcv_buff_size/list->ph_byte_size;
                        
                        if (rcv_buff_size&&snd_buff_size==0)
                        {
                            MPI_Wait(&request[0],&status[0]);
                            x_unpack(rcv_ph_buff,list,no_new_atms);
                            
                        }
                        else if (rcv_buff_size==0&&snd_buff_size)
                        {
                            MPI_Wait(&request[1],&status[1]);
                        }
                        else if (rcv_buff_size&&snd_buff_size)
                        {
                            MPI_Waitall(2,request,status);
                            x_unpack(rcv_ph_buff,list,no_new_atms);
                        }
                        
                        if(list->ph_vec_list[0]==0)
                        {
                            if((idir==0&&my_loc[idim]==tot_p_grid[idim]-1) ||
                               (idir==1&&my_loc[idim]==0))
                            {
                                int last_atm_no=natms+natms_ph;
                                int strt_atm_no=last_atm_no-no_new_atms;
                                //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
                                TYPE0* x;
                                vectors[0].ret(x);
                                int x_dim=vectors[0].dim;
                                
                                if(idir==0)
                                    for (int i=strt_atm_no;i<last_atm_no;i++)
                                        for(int j=0;j<=idim;j++)
                                            x[i*x_dim+j]+=H[idim][j];
                                else if(idir==1)
                                    for (int i=strt_atm_no;i<last_atm_no;i++)
                                        for(int j=0;j<=idim;j++)
                                            x[i*x_dim+j]-=H[idim][j];
                            }
                        }
                        iswap++;
                    }
                }
            }
            else
            {
                for(int idir=0;idir<2;idir++)
                {
                    for(int icomm=0;icomm<comm_need[idim];icomm++)
                    {
                        int last_atm_no=natms+natms_ph;
                        
                        for(int ivec=0;ivec<list->ph_no_vecs;ivec++)
                            for(int i=0;i<ph_lst->pos[iswap];i++)
                                vectors[list->ph_vec_list[ivec]]
                                .copy(last_atm_no+i,ph_lst->list[iswap][i]);
                        
                        if(list->ph_vec_list[0]==0)
                        {
                            //correction due to periodic boundary conditions
                            //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
                            TYPE0* x;
                            vectors[0].ret(x);
                            int x_dim=vectors[0].dim;
                            
                            if(idir==0)
                                for (int i=natms+natms_ph;i<natms+natms_ph+ph_lst->pos[iswap];i++)
                                    for(int j=0;j<=idim;j++)
                                        x[i*x_dim+j]+=H[idim][j];
                            else if(idir==1)
                                for (int i=natms+natms_ph;i<natms+natms_ph+ph_lst->pos[iswap];i++)
                                    for(int j=0;j<=idim;j++)
                                        x[i*x_dim+j]-=H[idim][j];
                            
                        }
                        natms_ph+=ph_lst->pos[iswap];
                        iswap++;
                    }
                }
            }
        }
    }
    
}
/*--------------------------------------------
 update number of vectors
 --------------------------------------------*/
void Atoms::update(int ivec)
{
    int iivec=ivec;
    update(&iivec,1,vectors[iivec].byte_size);
}
/*--------------------------------------------
 update number of vectors
 --------------------------------------------*/
void Atoms::update(int* vec_list,int no_vecs,int vec_byte_size)
{
    int iswap=0;
    natms_ph=0;
    
    if(no_vecs==1)
    {
        for(int idim=0;idim<dimension;idim++)
        {
            if(tot_p_grid[idim]!=1)
            {
                for(int idir=0;idir<2;idir++)
                {
                    for(int icomm=0;icomm<comm_need[idim];icomm++)
                    {
                        int snd_proc=neigh_p[idim][idir];
                        int rcv_proc=neigh_p[idim][1-idir];
                        int snd_buff_size=vec_byte_size
                        *ph_lst->pos[iswap];
                        int rcv_buff_size;
                        if(snd_buff_size>snd_ph_buff_capacity)
                        {
                            if(snd_ph_buff_capacity)
                                delete [] snd_ph_buff;
                            
                            CREATE1D(snd_ph_buff,snd_buff_size);
                            snd_ph_buff_capacity=snd_buff_size;
                        }
                        
                        x_pack(snd_ph_buff,vec_list,no_vecs,ph_lst->list[iswap],ph_lst->pos[iswap]);
                        
                        MPI_Request request[2];
                        MPI_Status status[2];
                        
                        MPI_Sendrecv(&snd_buff_size,1,MPI_INT,snd_proc,0,&rcv_buff_size,
                                     1,MPI_INT,rcv_proc,0,comm_world,&status[0]);
                        
                        if (rcv_buff_size)
                        {
                            MPI_Irecv(vectors[vec_list[0]].ret_vec(natms+natms_ph)
                                      ,rcv_buff_size,MPI_BYTE,rcv_proc,0,comm_world,&request[0]);
                            
                        }
                        if (snd_buff_size)
                        {
                            MPI_Isend(snd_ph_buff,snd_buff_size,
                                      MPI_BYTE,snd_proc,0,comm_world,&request[1]);
                        }
                        
                        int no_new_atms=rcv_buff_size/vec_byte_size;
                        
                        if (rcv_buff_size&&snd_buff_size==0)
                            MPI_Wait(&request[0],&status[0]);
                        else if (rcv_buff_size==0&&snd_buff_size)
                            MPI_Wait(&request[1],&status[1]);
                        else if (rcv_buff_size&&snd_buff_size)
                            MPI_Waitall(2,request,status);
                        
                        natms_ph+=no_new_atms;
                        if(vec_list[0]==0)
                        {
                            if((idir==0&&my_loc[idim]==tot_p_grid[idim]-1) ||
                               (idir==1&&my_loc[idim]==0))
                            {
                                int last_atm_no=natms+natms_ph;
                                int strt_atm_no=last_atm_no-no_new_atms;
                                //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
                                TYPE0* x;
                                vectors[0].ret(x);
                                int x_dim=vectors[0].dim;
                                
                                if(idir==0)
                                    for (int i=strt_atm_no;i<last_atm_no;i++)
                                        for(int j=0;j<=idim;j++)
                                            x[i*x_dim+j]+=H[idim][j];
                                else if(idir==1)
                                    for (int i=strt_atm_no;i<last_atm_no;i++)
                                        for(int j=0;j<=idim;j++)
                                            x[i*x_dim+j]-=H[idim][j];
                            }
                        }
                        iswap++;
                    }
                }
            }
            else
            {
                for(int idir=0;idir<2;idir++)
                {
                    for(int icomm=0;icomm<comm_need[idim];icomm++)
                    {
                        int last_atm_no=natms+natms_ph;
                        
                        for(int ivec=0;ivec<no_vecs;ivec++)
                            for(int i=0;i<ph_lst->pos[iswap];i++)
                                vectors[vec_list[ivec]]
                                .copy(last_atm_no+i,ph_lst->list[iswap][i]);
                        
                        if(vec_list[0]==0)
                        {
                            //correction due to periodic boundary conditions
                            //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
                            TYPE0* x;
                            vectors[0].ret(x);
                            int x_dim=vectors[0].dim;
                            
                            if(idir==0)
                                for (int i=natms+natms_ph;i<natms+natms_ph+ph_lst->pos[iswap];i++)
                                    for(int j=0;j<=idim;j++)
                                        x[i*x_dim+j]+=H[idim][j];
                            
                            else if(idir==1)
                                for (int i=natms+natms_ph;i<natms+natms_ph+ph_lst->pos[iswap];i++)
                                    for(int j=0;j<=idim;j++)
                                        x[i*x_dim+j]-=H[idim][j];
                            
                        }
                        natms_ph+=ph_lst->pos[iswap];
                        iswap++;
                    }
                }
            }
            
        }
    }
    else
    {
        for(int idim=0;idim<dimension;idim++)
        {
            if(tot_p_grid[idim]!=1)
            {
                for(int idir=0;idir<2;idir++)
                {
                    for(int icomm=0;icomm<comm_need[idir];icomm++)
                    {
                        int snd_proc=neigh_p[idim][idir];
                        int rcv_proc=neigh_p[idim][1-idir];
                        int snd_buff_size=vec_byte_size
                        *ph_lst->pos[iswap];
                        int rcv_buff_size;
                        if(snd_buff_size>snd_ph_buff_capacity)
                        {
                            if(snd_ph_buff_capacity)
                                delete [] snd_ph_buff;
                            CREATE1D(snd_ph_buff,snd_buff_size);
                            snd_ph_buff_capacity=snd_buff_size;
                        }
                        
                        x_pack(snd_ph_buff,vec_list,no_vecs,ph_lst->list[iswap],ph_lst->pos[iswap]);
                        
                        MPI_Request request[2];
                        MPI_Status status[2];
                        
                        MPI_Sendrecv(&snd_buff_size,1,
                                     MPI_INT,snd_proc,0,&rcv_buff_size,
                                     1,MPI_INT,rcv_proc,0,comm_world,&status[0]);
                        
                        if(rcv_buff_size>rcv_ph_buff_capacity)
                        {
                            if(rcv_ph_buff_capacity)
                                delete [] rcv_ph_buff;
                            CREATE1D(rcv_ph_buff,rcv_buff_size);
                            rcv_ph_buff_capacity=rcv_buff_size;
                        }
                        
                        if (rcv_buff_size)
                        {
                            MPI_Irecv(rcv_ph_buff,rcv_buff_size,
                                      MPI_BYTE,rcv_proc,0,comm_world,
                                      &request[0]);
                            
                        }
                        if (snd_buff_size)
                        {
                            MPI_Isend(snd_ph_buff,snd_buff_size,
                                      MPI_BYTE,snd_proc,0,comm_world,
                                      &request[1]);
                        }
                        
                        int no_new_atms=rcv_buff_size/vec_byte_size;
                        
                        if (rcv_buff_size&&snd_buff_size==0)
                        {
                            MPI_Wait(&request[0],&status[0]);
                            x_unpack(rcv_ph_buff,vec_list,no_vecs,no_new_atms);
                            
                        }
                        else if (rcv_buff_size==0&&snd_buff_size)
                        {
                            MPI_Wait(&request[1],&status[1]);
                        }
                        else if (rcv_buff_size&&snd_buff_size)
                        {
                            MPI_Waitall(2,request,status);
                            x_unpack(rcv_ph_buff,vec_list,no_vecs,no_new_atms);
                        }
                        
                        if(vec_list[0]==0)
                        {
                            if((idir==0&&my_loc[idim]==tot_p_grid[idim]-1) ||
                               (idir==1&&my_loc[idim]==0))
                            {
                                int last_atm_no=natms+natms_ph;
                                int strt_atm_no=last_atm_no-no_new_atms;
                                //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
                                TYPE0* x;
                                vectors[0].ret(x);
                                int x_dim=vectors[0].dim;
                                
                                if(idir==0)
                                    for (int i=strt_atm_no;i<last_atm_no;i++)
                                        for(int j=0;j<=idim;j++)
                                            x[i*x_dim+j]+=H[idim][j];
                                else if(idir==1)
                                    for (int i=strt_atm_no;i<last_atm_no;i++)
                                        for(int j=0;j<=idim;j++)
                                            x[i*x_dim+j]-=H[idim][j];
                            }
                        }
                        iswap++;
                    }
                }
            }
            else
            {
                for(int idir=0;idir<2;idir++)
                {
                    for(int icomm=0;icomm<comm_need[idim];icomm++)
                    {
                        int last_atm_no=natms+natms_ph;
                        
                        for(int ivec=0;ivec<no_vecs;ivec++)
                            for(int i=0;i<ph_lst->pos[iswap];i++)
                                vectors[vec_list[ivec]]
                                .copy(last_atm_no+i,ph_lst->list[iswap][i]);
                        
                        if(vec_list[0]==0)
                        {
                            //correction due to periodic boundary conditions
                            //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
                            TYPE0* x;
                            vectors[0].ret(x);
                            int x_dim=vectors[0].dim;
                            
                            if(idir==0)
                                for (int i=natms+natms_ph;i<natms+natms_ph+ph_lst->pos[iswap];i++)
                                    for(int j=0;j<=idim;j++)
                                        x[i*x_dim+j]+=H[idim][j];
                            else if(idir==1)
                                for (int i=natms+natms_ph;i<natms+natms_ph+ph_lst->pos[iswap];i++)
                                    for(int j=0;j<=idim;j++)
                                        x[i*x_dim+j]-=H[idim][j];
                            
                        }
                        natms_ph+=ph_lst->pos[iswap];
                        iswap++;
                    }
                }
            }
        }
    }
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::del(int idx)
{
    if(idx<0 || no_vecs-1<idx)
    {
        printf("Error: index should be between 0 and no_vecs-1");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    tot_byte_size-=vectors[idx].byte_size;
    
    AVec* new_vectors;
    new_vectors = new AVec[no_vecs-1];
    for(int i=0;i<no_vecs-1;i++)
        new_vectors[i].initiated=0;
    
    int k=0;
    for(int i=0;i<no_vecs;i++)
    {
        if(i!=idx)
        {
            move(&(new_vectors[k]),&(vectors[i]));
            k++;
        }
    }
    delete [] vectors;
    vectors=new_vectors;
    
    no_vecs--;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Atoms::find(const char* name)
{
    for (int i=0;i<no_vecs;i++)
        if(strcmp(vectors[i].name,name)==0)
            return i;
    
    error->abort("no such name: %s!!!",name);
    
    return -1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Atoms::find_exist(const char* name)
{
    for (int i=0;i<no_vecs;i++)
        if(strcmp(vectors[i].name,name)==0)
            return i;
    return -1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::copy(AVec* avec_dst,AVec* avec_src)
{
    if(avec_dst->initiated)
    {
        if(avec_dst->type==0)
            delete [] avec_dst->vec_0;
        else if(avec_dst->type==1)
            delete [] avec_dst->vec_1;
        else if(avec_dst->type==2)
            delete [] avec_dst->vec_2;
        else if(avec_dst->type==3)
            delete [] avec_dst->vec_3;
        else if(avec_dst->type==4)
            delete [] avec_dst->vec_4;
        else if(avec_dst->type==5)
            delete [] avec_dst->vec_5;
        else if(avec_dst->type==6)
            delete [] avec_dst->vec_6;
        else if(avec_dst->type==7)
            delete [] avec_dst->vec_7;
        else if(avec_dst->type==8)
            delete [] avec_dst->vec_8;
        else if(avec_dst->type==9)
            delete [] avec_dst->vec_9;
        else if(avec_dst->type==10)
            delete [] avec_dst->vec_10;
        else if(avec_dst->type==11)
            delete [] avec_dst->vec_11;
        else if(avec_dst->type==12)
            delete [] avec_dst->vec_12;
        
        delete [] avec_dst->name;
    }
    
    avec_dst->type=avec_src->type;
    avec_dst->ph=avec_src->ph;
    avec_dst->dim=avec_src->dim;
    avec_dst->byte_size=avec_src->byte_size;
    int size;
    int max_size;
    
    if (avec_dst->ph)
    {
        size=avec_dst->byte_size*(natms);
        max_size=avec_dst->dim*atm_vec_size;
    }
    else
    {
        size=avec_dst->byte_size*(natms_ph+natms);
        max_size=avec_dst->dim*atm_vec_ph_size;
    }
    
    if(avec_dst->type==0)
    {
        CREATE1D(avec_dst->vec_0,max_size);
        memcpy(avec_dst->vec_0,avec_src->vec_0,size);
    }
    else if(avec_dst->type==1)
    {
        CREATE1D(avec_dst->vec_1,max_size);
        memcpy(avec_dst->vec_1,avec_src->vec_1,size);
    }
    else if(avec_dst->type==2)
    {
        CREATE1D(avec_dst->vec_2,max_size);
        memcpy(avec_dst->vec_2,avec_src->vec_2,size);
    }
    else if(avec_dst->type==3)
    {
        CREATE1D(avec_dst->vec_3,max_size);
        memcpy(avec_dst->vec_3,avec_src->vec_3,size);
    }
    else if(avec_dst->type==4)
    {
        CREATE1D(avec_dst->vec_4,max_size);
        memcpy(avec_dst->vec_4,avec_src->vec_4,size);
    }
    else if(avec_dst->type==5)
    {
        CREATE1D(avec_dst->vec_5,max_size);
        memcpy(avec_dst->vec_5,avec_src->vec_5,size);
    }
    else if(avec_dst->type==6)
    {
        CREATE1D(avec_dst->vec_6,max_size);
        memcpy(avec_dst->vec_6,avec_src->vec_6,size);
    }
    else if(avec_dst->type==7)
    {
        CREATE1D(avec_dst->vec_7,max_size);
        memcpy(avec_dst->vec_7,avec_src->vec_7,size);
    }
    else if(avec_dst->type==8)
    {
        CREATE1D(avec_dst->vec_8,max_size);
        memcpy(avec_dst->vec_8,avec_src->vec_8,size);
    }
    else if(avec_dst->type==9)
    {
        CREATE1D(avec_dst->vec_9,max_size);
        memcpy(avec_dst->vec_9,avec_src->vec_9,size);
    }
    else if(avec_dst->type==10)
    {
        CREATE1D(avec_dst->vec_10,max_size);
        memcpy(avec_dst->vec_10,avec_src->vec_10,size);
    }
    else if(avec_dst->type==11)
    {
        CREATE1D(avec_dst->vec_11,max_size);
        memcpy(avec_dst->vec_11,avec_src->vec_11,size);
    }
    else if(avec_dst->type==12)
    {
        CREATE1D(avec_dst->vec_12,max_size);
        memcpy(avec_dst->vec_12,avec_src->vec_12,size);
    }
    
    int lngth=static_cast<int>(strlen(avec_src->name))+1;
    CREATE1D(avec_dst->name,lngth);
    memcpy(avec_dst->name,avec_src->name,lngth*sizeof(char));
    avec_dst->atms=avec_src->atms;
    avec_dst->initiated=1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::move(AVec* avec_dst,AVec* avec_src)
{
    if(avec_dst->initiated)
    {
        if(avec_dst->type==0)
            delete [] avec_dst->vec_0;
        else if(avec_dst->type==1)
            delete [] avec_dst->vec_1;
        else if(avec_dst->type==2)
            delete [] avec_dst->vec_2;
        else if(avec_dst->type==3)
            delete [] avec_dst->vec_3;
        else if(avec_dst->type==4)
            delete [] avec_dst->vec_4;
        else if(avec_dst->type==5)
            delete [] avec_dst->vec_5;
        else if(avec_dst->type==6)
            delete [] avec_dst->vec_6;
        else if(avec_dst->type==7)
            delete [] avec_dst->vec_7;
        else if(avec_dst->type==8)
            delete [] avec_dst->vec_8;
        else if(avec_dst->type==9)
            delete [] avec_dst->vec_9;
        else if(avec_dst->type==10)
            delete [] avec_dst->vec_10;
        else if(avec_dst->type==11)
            delete [] avec_dst->vec_11;
        else if(avec_dst->type==12)
            delete [] avec_dst->vec_12;
        
        delete [] avec_dst->name;
    }
    
    avec_dst->type=avec_src->type;
    avec_dst->ph=avec_src->ph;
    avec_dst->dim=avec_src->dim;
    avec_dst->byte_size=avec_src->byte_size;
    avec_dst->vec_0=avec_src->vec_0;
    avec_dst->vec_1=avec_src->vec_1;
    avec_dst->vec_2=avec_src->vec_2;
    avec_dst->vec_3=avec_src->vec_3;
    avec_dst->vec_4=avec_src->vec_4;
    avec_dst->vec_5=avec_src->vec_5;
    avec_dst->vec_6=avec_src->vec_6;
    avec_dst->vec_7=avec_src->vec_7;
    avec_dst->vec_8=avec_src->vec_8;
    avec_dst->vec_9=avec_src->vec_9;
    avec_dst->vec_10=avec_src->vec_10;
    avec_dst->vec_11=avec_src->vec_11;
    avec_dst->vec_12=avec_src->vec_12;
    avec_dst->name=avec_src->name;
    avec_dst->atms=avec_src->atms;
    avec_dst->initiated=1;
    avec_src->initiated=0;

}
/*--------------------------------------------
 this fucntion packs individual atoms from the
 vectors in VecLst; it also replaces the said
 atom with the last one
 --------------------------------------------*/
int Atoms::pack(char*& buff,int buff_pos,int& buff_capacity
,int iatm,class VecLst* list)
{
    
    if(buff_pos+list->byte_size>buff_capacity)
    {
        GROW(buff,buff_capacity,buff_pos+list->byte_size);
        buff_capacity=buff_pos+list->byte_size;
    }
    
    for(int i=0;i<list->no_vecs;i++)
    {
        buff_pos=vectors[list->vec_list[i]].pack(buff,buff_pos,iatm);
        vectors[list->vec_list[i]].copy(iatm,natms-1);
    }
    natms--;
    
    return buff_pos;
}
/*--------------------------------------------
 note: I have changed the order of list of
 vecs and atoms; in order to make it more
 consistant with pack which does it atom by
 atom; it might have been ok up to this point
 I have unpacked the atoms one by one.
 
 also the growing seems a bit out of place; it
 seems that it also grows the unneccessary
 vectors; the vectors that are not in the
 VecLst need not to be grown they only need to
 be resized;
 
 I fixed the growth; however, keep an eye on
 that.
 
 These whole changes might cause a shit storm.
 --------------------------------------------*/
int Atoms::unpack(char*& buff,int buff_pos
,int atm_list_size,class VecLst* list)
{
    int new_size=natms
    +atm_list_size;
    int old_size=atm_vec_size;
    
    if(new_size>old_size)
        grow(0,new_size-old_size,list);
    
    for(int iatm=0;iatm<atm_list_size;iatm++)
    {
        
        for(int i=0;i<list->no_vecs;i++)
        {

            buff_pos=vectors[list->vec_list[i]].unpack(buff,buff_pos);
        }
        natms++;
    }

    /*
    if(new_size>old_size)
        grow(new_size-old_size);
    for(int i=0;i<list->no_vecs;i++)
    {
        for(int iatm=0;iatm<atm_list_size;iatm++)
        {
            buff_pos=vectors[list->vec_list[i]].unpack(buff,buff_pos);
            natms++;
        }
        natms-=atm_list_size;
    }
    natms+=atm_list_size;
    */
    return buff_pos;
}
/*--------------------------------------------
 this function packs a list of atoms from the
 atomic vectors inside VecLst (from phantom
 vectors).
 --------------------------------------------*/
int Atoms::x_pack(char*& buff,class VecLst* list
,int* atm_list,int atm_list_size)
{
    int buff_pos=0;
    
    for(int i=0;i<list->ph_no_vecs;i++)
        buff_pos=vectors[list->ph_vec_list[i]]
        .x_pack(buff,buff_pos,atm_list,atm_list_size);
    
    return buff_pos;
}
/*--------------------------------------------
 this function unpacks the phantom atoms
 gathered by x_pack()
 --------------------------------------------*/
void Atoms::x_unpack(char*& buff
,class VecLst* list,int atm_list_size)
{
    int new_size=natms_ph+natms
    +atm_list_size;
    int old_size=atm_vec_ph_size;
    
    if(new_size>old_size)
    {
        grow(1,new_size-old_size,list);
        atm_vec_ph_size=new_size;
    }
    
    int buff_pose=0;
    
    for(int i=0;i<list->ph_no_vecs;i++)
        buff_pose=vectors[list->ph_vec_list[i]]
        .x_unpack(buff,buff_pose,atm_list_size);
    
    natms_ph+=atm_list_size;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int Atoms::x_pack(char*& buff,int* vec_list,int no_vecs
,int* atm_list,int atm_list_size)
{
    int buff_pos=0;
    for(int i=0;i<no_vecs;i++)
        buff_pos=vectors[vec_list[i]]
        .x_pack(buff,buff_pos,atm_list,atm_list_size);
    
    return buff_pos;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::x_unpack(char*& buff,int* vec_list
,int no_vecs,int atm_list_size)
{
    int buff_pose=0;
    
    for(int i=0;i<no_vecs;i++)
        buff_pose=vectors[vec_list[i]]
        .x_unpack(buff,buff_pose,atm_list_size);
    
    natms_ph+=atm_list_size;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::store_0()
{
    TYPE0* x;
    vectors[0].ret(x);
    TYPE0* x_0;
    vectors[1].ret(x_0);
    
    int x_dim=vectors[0].dim;
    int x_0_dim=vectors[1].dim;
    int x_comp=0;
    int x_0_comp=0;
    for(int i=0;i<natms;i++)
    {
        for(int idim=0;idim<dimension;idim++)
            x_0[x_0_comp+idim]=x[x_comp+idim];
        x_comp+=x_dim;
        x_0_comp+=x_0_dim;
    }
}
/*--------------------------------------------
 please note that if the box's size changes,
 this function assumes the you have already 
 inversed the H matrix and calculated the B
 matrix; if you have truble inersing your H
 matrix, you can use the invert(double**,
 double**,int). however, please be advised 
 this procedure is numerical and you need to
 adjust TOLERANCE value.
 --------------------------------------------*/
void Atoms::update_0(int box_change,int neigh_chk
,class VecLst* list)
{
    //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
    TYPE0* x;
    vectors[0].ret(x);
    //TYPE0* x_0=(TYPE0*)vectors[1].ret_vec();
    TYPE0* x_0;
    vectors[1].ret(x_0);
    
    int x_dim=vectors[0].dim;
    int x_0_dim=vectors[1].dim;


    TYPE0 dx;
    TYPE0 sq_dx;
    
    int iatm=0;
    int check=0;
    TYPE0 sq_skin=skin*skin;
    
    int x_comp=0;
    int x_0_comp=0;
    while (iatm<natms&& check==0)
    {
        
        sq_dx=0.0;
        for(int idim=0;idim<dimension;idim++)
        {
            dx=x[x_comp+idim]-x_0[x_0_comp+idim];
            sq_dx+=dx*dx;
        }
        if(sq_dx>sq_skin)
            check=1;
        
        iatm++;
        x_comp+=x_dim;
        x_0_comp+=x_0_dim;
    }
    
    int check_all;
    MPI_Allreduce(&check,&check_all,1,
    MPI_INT,MPI_MAX,comm_world);

    
    if(check_all)
    {
        x2s(natms);
        xchng_comm(list);
        ph_setup(box_change,list);
        if(neigh_chk)
        {
            neighbor->create_list(box_change,1);
            /*no need to convert s back to x,
             the neighbor list takes care of it*/
        }
        else
        {
            s2x(natms+natms_ph);
        }

        x_comp=0;
        x_0_comp=0;
 
        vectors[0].ret(x);
        vectors[1].ret(x_0);
        for(int i=0;i<natms;i++)
        {
            for(int idim=0;idim<dimension;idim++)
                x_0[x_0_comp+idim]=x[x_comp+idim];
            x_comp+=x_dim;
            x_0_comp+=x_0_dim;
        }
    }
    else
    {
        

        update(list->update_every_ph_vec_list
               ,list->update_every_ph_no_vecs
               ,list->update_every_ph_byte_size);
        

    }
}
/*--------------------------------------------
 transform x 2 s
 --------------------------------------------*/
void Atoms::x2s(int no)
{
    //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
    TYPE0* x;
    vectors[0].ret(x);
    int x_dim=vectors[0].dim;
    int icomp=0;
    for(int i=0;i<no;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            x[icomp+j]=x[icomp+j]*B[j][j];
            for(int k=j+1;k<dimension;k++)
                x[icomp+j]+=B[k][j]*x[icomp+k];
            
        }
        
        for(int j=0;j<dimension;j++)
        {
            while(x[icomp+j]<0.0)
                x[icomp+j]++;
            while(x[icomp+j]>=1.0)
                x[icomp+j]--;
            
        }
               
        icomp+=x_dim;
    }
}
/*--------------------------------------------
 transform x 2 s
 --------------------------------------------*/
void Atoms::x2s_no_correction(int no)
{
    //TYPE0* x=(TYPE0*)vectors[0].ret_vec();
    TYPE0* x;
    vectors[0].ret(x);
    int x_dim=vectors[0].dim;
    int icomp=0;
    for(int i=0;i<no;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            x[icomp+j]=x[icomp+j]*B[j][j];
            for(int k=j+1;k<dimension;k++)
                x[icomp+j]+=B[k][j]*x[icomp+k];
            
        }
        
        icomp+=x_dim;
    }
}
/*--------------------------------------------
 transform s 2 x
 --------------------------------------------*/
void Atoms::s2x(int no)
{
    //TYPE0* s=(TYPE0*)vectors[0].ret_vec();
    TYPE0* s;
    vectors[0].ret(s);
    int s_dim=vectors[0].dim;
    int icomp=0;
    for(int i=0;i<no;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            s[icomp+j]=s[icomp+j]*H[j][j];
            for(int k=j+1;k<dimension;k++)
                s[icomp+j]+=H[k][j]*s[icomp+k];
        }
        icomp+=s_dim;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
/*
void Atoms::grow(int no_atms)
{
    for(int i=0;i<no_vecs;i++)
        vectors[i].grow(no_atms);
    atm_vec_size+=no_atms;
    atm_vec_ph_size+=no_atms;
}*/
/*--------------------------------------------
 
 --------------------------------------------*/
/*void Atoms::grow(int ph,int no_atms)
{
    for(int i=0;i<no_vecs;i++)
        if(vectors[i].ph==ph)
            vectors[i].grow(no_atms);
    
    if(ph)
        atm_vec_ph_size+=no_atms;
    else
        atm_vec_size+=no_atms;
}*/
/*--------------------------------------------
 
 --------------------------------------------*/
void Atoms::grow(int ph,int no_atms
,class VecLst* list)
{
    int* vec_list;
    int vec_list_size;
    
    if(ph)
    {
        vec_list=list->ph_vec_list;
        vec_list_size=list->ph_no_vecs;
    }
    else
    {
        vec_list=list->vec_list;
        vec_list_size=list->no_vecs;
    }
    
    int ivec=0;
    for(int i=0;i<no_vecs;i++)
    {
        if(i==vec_list[ivec])
        {
            vectors[vec_list[ivec]].grow(no_atms);
            ivec++;
        }
        else
        {
            if(ph)
            {
                if(vectors[i].ph)
                    vectors[i].resize(no_atms);
            }
            else
                vectors[i].resize(no_atms);
        }
    }
    
    if(ph)
        atm_vec_ph_size+=no_atms;
    else
    {
        atm_vec_size+=no_atms;
        atm_vec_ph_size+=no_atms;
    }
    

}

/*--------------------------------------------
 autogrid the domain
 --------------------------------------------*/
void Atoms::hard_auto_grid_proc(TYPE0 f)
{
    
    if(f<=0 || f>1)
        error->abort("inter-node efficiency "
                     "should be between 0.0 and 1.0");
    
    int eq_p_per_n=1;
    for(int i=0;i<tot_n;i++)
        if(p_per_n[i]!=p_per_n[0])
            eq_p_per_n=0;
    
    if(eq_p_per_n || f==1)
    {
        auto_grid_proc();
        return;
    }
    
    inter_n_efficency=f;
    
    energy_all=-1.0;
    
    int* list;
    CREATE1D(list,dimension);
    
    fac_list_size=0;
    fac(tot_p,dimension,0,list);
    delete [] list;
    
    int* denom;
    CREATE1D(denom,dimension);
    int* my_tmp_loc;
    CREATE1D(my_tmp_loc,dimension);
    CREATE1D(vols,fac_list_size);
    
    CREATE1D(areas,fac_list_size);
    CREATE1D(nxt_p,fac_list_size);
    CREATE1D(prv_p,fac_list_size);
    for(int i=0;i<fac_list_size;i++)
    {
        CREATE1D(areas[i],dimension);
        CREATE1D(nxt_p[i],dimension);
        CREATE1D(prv_p[i],dimension);
    }
    
    for(int ifac=0;ifac<fac_list_size;ifac++)
    {
        vols[ifac]=1.0;
        for(int i=0;i<dimension;i++)
        {
            vols[ifac]*=H[i][i]/static_cast<TYPE0>(fac_list[ifac][i]);
            areas[ifac][i]=1.0;
            for(int j=0;j<dimension;j++)
                if(i!=j)
                    areas[ifac][i]*=H[j][j]/static_cast<TYPE0>(fac_list[ifac][j]);
        }
        
        for(int i=0;i<dimension;i++)
        {
            int no=1;
            for(int j=0;j<i;j++)
                no*=fac_list[ifac][j];
            denom[i]=no;
        }
        
        int ttmp=my_p_no;
        for(int i=dimension-1;i>-1;i--)
        {
            my_tmp_loc[i]=ttmp/denom[i];
            ttmp-=my_tmp_loc[i]*denom[i];
        }
        
        int tmp0;
        int tmp1;
        
        for(int i=0;i<dimension;i++)
        {
            
            tmp0=my_tmp_loc[i]-1;
            if(tmp0==-1)
                tmp0=fac_list[ifac][i]-1;
            
            tmp1=my_tmp_loc[i]+1;
            if(tmp1==fac_list[ifac][i])
                tmp1=0;
            
            int nop=0;
            int nom=0;
            
            for(int j=0;j<dimension;j++)
            {
                if(i==j)
                {
                    nom+=tmp0*denom[j];
                    nop+=tmp1*denom[j];
                }
                else
                {
                    nom+=my_tmp_loc[j]*denom[j];
                    nop+=my_tmp_loc[j]*denom[j];
                }
            }
            prv_p[ifac][i]=nom;
            nxt_p[ifac][i]=nop;
        }
    }
    
    delete [] denom;
    delete [] my_tmp_loc;
    
    
    CREATE1D(res_perm,tot_p);
    CREATE1D(res_grid,dimension);
    
    
    
    comb(p_per_n,tot_n);
    
    
    
    for(int i=0;i<fac_list_size;i++)
    {
        delete [] areas[i];
        delete [] fac_list[i];
        delete [] nxt_p[i];
        delete [] prv_p[i];
    }
    delete [] vols;
    delete [] areas;
    delete [] fac_list;
    delete [] nxt_p;
    delete [] prv_p;
    fac_list_size=0;
    
    
    
    int my_p_in_my_node=-1;
    for(int i=0;i<p_per_n[my_n_no];i++)
        if(n_p_grid[my_n_no][i]==my_p_no)
            my_p_in_my_node=i;
    
    int my_place_in_perm=0;
    int ino=0;
    while (res_perm[my_place_in_perm]!=my_n_no
           ||
           my_p_in_my_node!=ino)
    {
        if(res_perm[my_place_in_perm]==my_n_no)
            ino++;
        my_place_in_perm++;
    }
    
    for(int i=0;i<dimension;i++)
        tot_p_grid[i]=res_grid[i];
    
    CREATE1D(denom,dimension);
    for(int i=0;i<dimension;i++)
    {
        int no=1;
        for(int j=0;j<i;j++)
            no*=tot_p_grid[j];
        denom[i]=no;
    }
    
    int ttmp=my_place_in_perm;
    for(int i=dimension-1;i>-1;i--)
    {
        my_loc[i]=ttmp/denom[i];
        ttmp-=my_loc[i]*denom[i];
    }
    
    int tmp0,tmp1;
    
    for(int i=0;i<dimension;i++)
    {
        
        tmp0=my_loc[i]-1;
        if(tmp0==-1)
            tmp0=tot_p_grid[i]-1;
        
        tmp1=my_loc[i]+1;
        if(tmp1==tot_p_grid[i])
            tmp1=0;
        
        int nop=0;
        int nom=0;
        
        for(int j=0;j<dimension;j++)
        {
            if(i==j)
            {
                nom+=tmp0*denom[j];
                nop+=tmp1*denom[j];
            }
            else
            {
                nom+=my_loc[j]*denom[j];
                nop+=my_loc[j]*denom[j];
            }
        }
        
        
        ino=0;
        for(int k=0;k<res_perm[nom];k++)
            if(res_perm[k]==res_perm[nom])
                ino++;
        neigh_p[i][0]=n_p_grid[res_perm[nom]][ino];
        
        ino=0;
        for(int k=0;k<res_perm[nop];k++)
            if(res_perm[k]==res_perm[nop])
                ino++;
        neigh_p[i][1]=n_p_grid[res_perm[nop]][ino];
        
    }
    
    delete [] denom;
    delete [] res_perm;
    delete [] res_grid;
    
    for(int i=0;i<dimension;i++)
    {
        s_lo[i]=
        static_cast<TYPE0>(my_loc[i])
        /static_cast<TYPE0>(tot_p_grid[i]);
        
        s_hi[i]=
        (static_cast<TYPE0>(my_loc[i])+1.0)
        /static_cast<TYPE0>(tot_p_grid[i]);
    }
    
    if(my_p_no==0)
        fprintf(output,"Autogrid performed:");
    
    for(int i=0;i<dimension-1;i++)
        if(my_p_no==0)
            fprintf(output," %d by",tot_p_grid[i]);
    if(my_p_no==0)
        fprintf(output," %d",tot_p_grid[dimension-1]);
    if(my_p_no==0)
        fprintf(output,"\n\n");
}
/*--------------------------------------------
 autogrid the domain
 --------------------------------------------*/
void Atoms::auto_grid_proc()
{
    
    int eq_p_per_n=1;
    for(int i=0;i<tot_n;i++)
        if(p_per_n[i]!=p_per_n[0])
            eq_p_per_n=0;

    TYPE0 vol;
    TYPE0* area;
    CREATE1D(area,dimension);
    int prin_dimension=0;
    
    if(eq_p_per_n)
    {
        for(int i=0;i<dimension;i++)
        {
            area[i]=1.0;
            for(int j=0;j<dimension;j++)
                if(i!=j)
                    area[i]*=H[j][j];
        }
        
        TYPE0 min_area=-1.0;
        for(int i=0;i<dimension;i++)
            if(min_area<0.0 || area[i]<min_area)
            {
                min_area=area[i];
                prin_dimension=i;
            }
        
    }
    
    int* list;
    CREATE1D(list,dimension);
    
    fac_list_size=0;
    fac(tot_p,dimension,0,list);
    
    TYPE0 ratio=-1.0;
    TYPE0 tmp_ratio;

    if(eq_p_per_n && tot_n>1)
    {
        for(int ifac=0;ifac<fac_list_size;ifac++)
        {
            if(fac_list[ifac][prin_dimension]==tot_n)
            {
                vol=1.0;
                for(int i=0;i<dimension;i++)
                {
                    vol*=H[i][i]/static_cast<TYPE0>(fac_list[ifac][i]);
                    area[i]=1.0;
                    for(int j=0;j<dimension;j++)
                        if(i!=j)
                            area[i]*=H[j][j]/static_cast<TYPE0>(fac_list[ifac][j]);
                }
                
                tmp_ratio=0.0;
                for(int i=0;i<dimension;i++)
                    tmp_ratio+=2.0*area[i];
                
                tmp_ratio=tmp_ratio/vol;
                if(tmp_ratio<ratio||ratio<0.0)
                {
                    ratio=tmp_ratio;
                    for(int i=0;i<dimension;i++)
                        tot_p_grid[i]=fac_list[ifac][i];
                }
            }
            
        }
        
        
        for(int i=0;i<fac_list_size;i++)
            delete [] fac_list[i];
        delete [] fac_list;
        
        
        if(dimension>1)
        {
            int* tmp_tot_p_grid;
            CREATE1D(tmp_tot_p_grid,dimension-1);
            
            int* dim_indx;
            CREATE1D(dim_indx,dimension-1);
            
            int* my_tmp_loc;
            CREATE1D(my_tmp_loc,dimension-1);
            
            int* my_tmp_denom;
            CREATE1D(my_tmp_denom,dimension-1);
            
            
            int pos=0;
            
            for(int i=0;i<dimension;i++)
            {
                if(i!=prin_dimension)
                {
                    tmp_tot_p_grid[pos]=tot_p_grid[i];
                    dim_indx[pos]=i;
                    pos++;
                }
            }
            
            int my_p_in_my_node=-1;
            for(int i=0;i<p_per_n[my_n_no];i++)
                if(n_p_grid[my_n_no][i]==my_p_no)
                    my_p_in_my_node=i;

            for(int i=0;i<dimension-1;i++)
            {
                
                int no=1;
                for(int j=0;j<i;j++)
                    no*=tmp_tot_p_grid[j];
                my_tmp_denom[i]=no;
            }
            
            
            int ttmp=my_p_in_my_node;
            for(int i=dimension-2;i>-1;i--)
            {
                my_tmp_loc[i]=ttmp/my_tmp_denom[i];
                ttmp-=my_tmp_loc[i]*my_tmp_denom[i];
            }
            
            
            for(int i=0;i<dimension-1;i++)
            {
                my_loc[dim_indx[i]]=my_tmp_loc[i];
            }
            my_loc[prin_dimension]=my_n_no;
            
            
            
            int nop;
            int nom;
            int tmp0,tmp1;
            for(int i=0;i<dimension-1;i++)
            {
                
                tmp0=my_tmp_loc[i]-1;
                if(tmp0==-1)
                    tmp0=tmp_tot_p_grid[i]-1;
                
                tmp1=my_tmp_loc[i]+1;
                if(tmp1==tmp_tot_p_grid[i])
                    tmp1=0;
                
                nop=0;
                nom=0;
                
                for(int j=0;j<dimension-1;j++)
                {
                    if(i==j)
                    {
                        nom+=tmp0*my_tmp_denom[j];
                        nop+=tmp1*my_tmp_denom[j];
                    }
                    else
                    {
                        nom+=my_tmp_loc[j]*my_tmp_denom[j];
                        nop+=my_tmp_loc[j]*my_tmp_denom[j];
                    }
                }
                
                neigh_p[dim_indx[i]][1]=n_p_grid[my_n_no][nop];
                neigh_p[dim_indx[i]][0]=n_p_grid[my_n_no][nom];
            }
            
            tmp0=my_n_no-1;
            if(tmp0==-1)
                tmp0=tot_n-1;
            neigh_p[prin_dimension][0]=n_p_grid[tmp0][my_p_in_my_node];
            
            tmp1=my_n_no+1;
            if(tmp1==tot_n)
                tmp1=0;
            neigh_p[prin_dimension][1]=n_p_grid[tmp1][my_p_in_my_node];
            
            
            delete [] dim_indx;
            delete [] tmp_tot_p_grid;
            delete [] my_tmp_loc;
            delete [] my_tmp_denom;
            delete [] area;
            
            
        }
        else
        {
            tot_p_grid[0]=tot_p;
            
            int my_p_in_my_node=-1;
            for(int i=0;i<p_per_n[my_n_no];i++)
                if(n_p_grid[my_n_no][i]==my_p_no)
                    my_p_in_my_node=i;
            
            
            my_loc[0]=p_per_n[my_n_no]*my_n_no;
            my_loc[0]+=my_p_in_my_node;
            
            int tmp0p,tmp1p,tmp0n,tmp1n;
            tmp0p=my_p_in_my_node;
            tmp1p=my_p_in_my_node;
            tmp0n=my_n_no;
            tmp1n=my_n_no;
            
            tmp0p--;
            if(tmp0p==-1)
            {
                tmp0p=p_per_n[my_n_no]-1;
                tmp0n--;
                if(tmp0n==-1)
                    tmp0n=tot_n-1;
            }
            
            tmp1p++;
            if(tmp1p==p_per_n[my_n_no])
            {
                tmp1p=0;
                tmp1n++;
                if(tmp1n==tot_n)
                    tmp1n=0;
            }
            
            neigh_p[0][0]=n_p_grid[tmp0n][tmp0p];
            neigh_p[0][1]=n_p_grid[tmp1n][tmp1p];
            
        }

    }
    else
    {
        for(int ifac=0;ifac<fac_list_size;ifac++)
        {
            vol=1.0;
            for(int i=0;i<dimension;i++)
            {
                vol*=H[i][i]/static_cast<TYPE0>(fac_list[ifac][i]);
                area[i]=1.0;
                for(int j=0;j<dimension;j++)
                    if(i!=j)
                        area[i]*=H[j][j]/static_cast<TYPE0>(fac_list[ifac][j]);
            }
                
            tmp_ratio=0.0;
            for(int i=0;i<dimension;i++)
                tmp_ratio+=2.0*area[i];
            
            tmp_ratio=tmp_ratio/vol;
            if(tmp_ratio<ratio||ratio<0.0)
            {
                ratio=tmp_ratio;
                for(int i=0;i<dimension;i++)
                    tot_p_grid[i]=fac_list[ifac][i];
            }
            
        }
        
        delete [] area;
        
        for(int i=0;i<fac_list_size;i++)
            delete [] fac_list[i];
        delete [] fac_list;
        
        
        MPI_Comm cartesian;
        for(int i=0;i<dimension;i++)
            list[i]=1;
        MPI_Cart_create(comm_world,dimension,tot_p_grid,list,1,&cartesian);
        MPI_Cart_get(cartesian,dimension,tot_p_grid,list,my_loc);
        for(int i=0;i<dimension;i++)
            MPI_Cart_shift(cartesian,i,1,&neigh_p[i][0],&neigh_p[i][1]);
        
        MPI_Comm_free(&cartesian);
        
    }
    
    delete [] list;
    
    for(int i=0;i<dimension;i++)
    {
        s_lo[i]=
        static_cast<TYPE0>(my_loc[i])
        /static_cast<TYPE0>(tot_p_grid[i]);
        
        s_hi[i]=
        (static_cast<TYPE0>(my_loc[i])+1.0)
        /static_cast<TYPE0>(tot_p_grid[i]);
    }
    if(my_p_no==0)
        fprintf(output,"Autogrid performed:");
    
    for(int i=0;i<dimension-1;i++)
        if(my_p_no==0)
            fprintf(output," %d by",tot_p_grid[i]);
    if(my_p_no==0)
        fprintf(output," %d",tot_p_grid[dimension-1]);
    if(my_p_no==0)
        fprintf(output,"\n\n");
    
    
}
/*--------------------------------------------
 autogrid the domain
 --------------------------------------------*/
void Atoms::fac(int no,int dim,int pos,int*& tmp_list)
{
    if(dim>1)
    {
        for(int i=1;i<=no;i++)
        {
            if(no%i==0)
            {
                tmp_list[pos]=i;
                fac(no/i,dim-1,pos+1,tmp_list);
            }
        }
    }
    else
    {
        tmp_list[pos]=no;
        GROW(fac_list,fac_list_size,fac_list_size+1);
        CREATE1D(fac_list[fac_list_size],pos+1);
        for(int i=0;i<=pos;i++)
            fac_list[fac_list_size][i]=tmp_list[i];
        
        fac_list_size++;
    }
}
/*--------------------------------------------
 combination
 --------------------------------------------*/
void Atoms::comb(int* no_list, int list_size)
{
    int tot=0;
    for(int i=0;i<list_size;i++)
        tot+=no_list[i];
    
    int* tmp0;
    CREATE1D(tmp0,tot);
    int* tmp1;
    CREATE1D(tmp1,tot);
    comb_rec(no_list[0],tot,0,no_list[0],tmp0,tmp1,no_list,0,tot);
    delete [] tmp0;
    delete [] tmp1;
}
/*--------------------------------------------
 recursive combination
 --------------------------------------------*/
void Atoms::comb_rec(int no,int tot,int pos
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
        
        if(y[y_pos]!=tot)
        {
            int noo=y[y_pos+1];
            int tott=tot-y[y_pos];
            
            int* h=lvl+y[y_pos];
            comb_rec(noo,tott,0,noo,h,tmp,y,y_pos+1,y_tot);
            for(int i=0;i<y_tot;i++)
                if(tmp[i]==y_pos)
                    tmp[i]=-1;
        }
        else
        {
            
            
            for(int ifac=0;ifac<fac_list_size;ifac++)
            {
                TYPE0 energy=0.0;
                
                energy=0.0;
                for(int i=0;i<dimension;i++)
                {
                    if(tmp[my_p_no]-tmp[prv_p[ifac][i]])
                        energy+=static_cast<TYPE0>(areas[ifac][i])/inter_n_efficency;
                    else
                        energy+=static_cast<TYPE0>(areas[ifac][i]);
                    
                    if(tmp[my_p_no]-tmp[nxt_p[ifac][i]])
                        energy+=static_cast<TYPE0>(areas[ifac][i])/inter_n_efficency;
                    else
                        energy+=static_cast<TYPE0>(areas[ifac][i]);
                    
                }
                energy=energy/vols[ifac];
                
                TYPE0 energy_tot=0.0;
                MPI_Allreduce(&energy,&energy_tot,1,MPI_TYPE0,MPI_SUM,comm_world);
                
                if(energy_tot<energy_all||energy_all<0.0)
                {
                    energy_all=energy_tot;
                    for(int i=0;i<dimension;i++)
                        res_grid[i]=fac_list[ifac][i];
                    
                    for(int i=0;i<tot_p;i++)
                        res_perm[i]=tmp[i];
                }
            }
            
            for(int i=0;i<y_tot;i++)
                if(tmp[i]==y_pos)
                    tmp[i]=-1;
        }
        
    }
}
/*--------------------------------------------
 find the prime factors and powers of a number
 --------------------------------------------*/
void Atoms::prime_factorize(int no,int*& primes,int*& pow)
{
    int chk;
    int j,k,l;
    int size=0;
    
    for(int i=2;i<=sqrt(no);i++)
    {
        if(no%i==0)
        {
            chk=1;
            
            j=2;
            while (chk && j<=sqrt(i))
            {
                if(i%j==0)
                    chk=0;
                j++;
            }
            
            if(chk)
            {
                k=0;
                l=no;
                while(l%i==0)
                {
                    l=l/i;
                    k++;
                }
                
                GROW(primes,size,size+1);
                GROW(pow,size,size+1);
                primes[size]=i;
                pow[size]=k;
                size++;
            }
        }
    }
    
    if(size==0)
    {
        GROW(primes,size,size+1);
        GROW(pow,size,size+1);
        primes[size]=no;
        pow[size]=1;
        size++;
    }
}
/*--------------------------------------------
 inversion funtion to calculate B whenever H
 is changed
 --------------------------------------------*/
void Atoms::invert(TYPE0** A,TYPE0** Ainv,int dim)
{
    if(dim==0)
        return;
    
    TYPE0** ATA;
    CREATE2D(ATA,dim,dim);
    TYPE0* c;
    CREATE1D(c,dim);
    TYPE0* x;
    CREATE1D(x,dim);
    TYPE0* g;
    CREATE1D(g,dim);
    TYPE0* g0;
    CREATE1D(g0,dim);
    TYPE0* h;
    CREATE1D(h,dim);
    TYPE0 a0,a1,alpha;
    TYPE0 g0g0,gg,gg0,ratio;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            Ainv[i][j]=ATA[i][j]=0.0;
    
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            for(int k=0;k<dim;k++)
                ATA[i][j]+=A[k][i]*A[k][j];
    
    for(int itry=0;itry<dim;itry++)
    {
        for(int i=0;i<dim;i++)
        {
            c[i]=A[itry][i];
            x[i]=c[i];
        }
        
        
        g0g0=0.0;
        for(int i=0;i<dim;i++)
        {
            h[i]=2.0*c[i];
            for(int j=0;j<dim;j++)
                h[i]-=2.0*ATA[i][j]*x[j];
            g[i]=h[i];
            g0g0+=h[i]*h[i];
        }
        
        int jtry=0;
        double error=1.0;
        while(jtry<dim+1 && error!=0.0)
        {
            
            if(g0g0==0.0)
            {
                error=0.0;
                continue;
            }
            
            
            a0=0.0;
            a1=0.0;
            for(int i=0;i<dim;i++)
            {
                a0+=h[i]*g[i];
                for(int j=0;j<dim;j++)
                    a1+=h[i]*ATA[i][j]*h[j];
            }
            if(a1==0.0)
            {
                error=0.0;
                continue;
            }
            alpha=0.5*a0/a1;
            
            for(int i=0;i<dim;i++)
                x[i]+=alpha*h[i];
            
            //cout << "chk 3" << endl;
            
            gg=0.0;
            gg0=0.0;
            for(int i=0;i<dim;i++)
            {
                g[i]=2.0*c[i];
                for(int j=0;j<dim;j++)
                    g[i]-=2.0*ATA[i][j]*x[j];
                gg+=g[i]*g[i];
                gg0+=g0[i]*g[i];
            }
            
            //cout << "chk 4" << endl;
            ratio=(gg-gg0)/g0g0;
            g0g0=gg;
            
            
            for(int i=0;i<dim;i++)
            {
                h[i]=ratio*h[i]+g[i];
                g0[i]=g[i];
            }

            
            error=0.0;
            for(int i=0;i<dim;i++)
            {
                for(int j=0;j<dim;j++)
                    error+=x[i]*ATA[i][j]*x[i];
                error-=2*c[i]*x[i];
            }
            error++;
            
            jtry++;
        }
        
        for(int i=0;i<dim;i++)
            Ainv[i][itry]=x[i];
    }
    
    for(int i=0;i<dim;i++)
        delete [] ATA[i];
    delete [] ATA;
    delete [] c;
    delete [] x;
    delete [] g;
    delete [] g0;
    delete [] h;
}
/*--------------------------------------------
 inversion funtion to calculate B whenever H
 is changed
 --------------------------------------------*/
void Atoms::invert_lower_triangle(TYPE0** A,TYPE0** Ainv,int dim)
{
    invert(A,Ainv,dim);
    for(int i=0;i<dim;i++)
        for(int j=i+1;j<dim;j++)
            Ainv[i][j]=0.0;
    
}
/*--------------------------------------------
 add skin command
 --------------------------------------------*/
void Atoms::add_skin(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong command");
    TYPE0 s=atof(args[1]);
    if(s<=0.0)
        error->abort("skin cannot be equal or less than zero");
    skin=s;
}

/*----------------------------------------------------------------------------------------*/

/*--------------------------------------------
 
 --------------------------------------------*/
VecLst::VecLst(MAPP* mapp,int no_args,...)
:InitPtrs(mapp)
{
    va_list args;
    va_start (args,no_args);
    
    no_vecs=no_args;
    CREATE1D(vec_list,no_vecs);
    
    for(int i=0;i<no_args;i++)
    {
        vec_list[i]=va_arg(args,int);
        if(vec_list[i]<0 || vec_list[i] >(atoms->no_vecs-1))
        {
            error->abort("wrong index for the vectors: %i",vec_list[i]);
        }
        
    }
    
    va_end (args);
    
    int tmp;
    for(int i=0;i<no_vecs;i++)
    {
        for(int j=i+1;j<no_vecs;j++)
        {
            if(vec_list[i]==vec_list[j])
                error->abort("Error: in vector_list two "
                "vectors with the same index cannot exist\n");
            
            if(vec_list[i]>vec_list[j])
            {
                tmp=vec_list[i];
                vec_list[i]=vec_list[j];
                vec_list[j]=tmp;
            }
        }
    }
    
    byte_size=0;
    for(int i=0;i<no_vecs;i++)
        byte_size+=atoms->vectors[vec_list[i]].byte_size;
    
    ph_no_vecs=0;
    for(int i=0;i<no_vecs;i++)
        if(atoms->vectors[vec_list[i]].ph)
            ph_no_vecs++;
    
    CREATE1D(ph_vec_list,ph_no_vecs);
    
    int pos=0;
    ph_byte_size=0;
    for(int i=0;i<no_vecs;i++)
        if(atoms->vectors[vec_list[i]].ph)
        {
            ph_vec_list[pos++]=vec_list[i];
            ph_byte_size+=atoms->vectors[vec_list[i]].byte_size;
        }
    
    update_every_ph_no_vecs=0;
    update_every_ph_byte_size=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
VecLst::VecLst(MAPP* mapp,int* list,int no)
:InitPtrs(mapp)
{
    
    no_vecs=no;
    CREATE1D(vec_list,no_vecs);
    for(int i=0;i<no_vecs;i++)
        vec_list[i]=list[i];
    
    int tmp;
    for(int i=0;i<no_vecs;i++)
    {
        for(int j=i+1;j<no_vecs;j++)
        {
            if(vec_list[i]==vec_list[j])
                error->abort("Error: in vector_list two "
                             "vectors with the same index cannot exist\n");
            
            if(vec_list[i]>vec_list[j])
            {
                tmp=vec_list[i];
                vec_list[i]=vec_list[j];
                vec_list[j]=tmp;
            }
        }
    }
    
    byte_size=0;
    for(int i=0;i<no_vecs;i++)
        byte_size+=atoms->vectors[vec_list[i]].byte_size;
    
    ph_no_vecs=0;
    for(int i=0;i<no_vecs;i++)
        if(atoms->vectors[vec_list[i]].ph)
            ph_no_vecs++;
    
    CREATE1D(ph_vec_list,ph_no_vecs);
    
    int pos=0;
    ph_byte_size=0;
    for(int i=0;i<no_vecs;i++)
        if(atoms->vectors[vec_list[i]].ph)
        {
            ph_vec_list[pos++]=vec_list[i];
            ph_byte_size+=atoms->vectors[vec_list[i]].byte_size;
        }
    
    update_every_ph_no_vecs=0;
    update_every_ph_byte_size=0;
}
/*--------------------------------------------
 
 --------------------------------------------*/
VecLst::~VecLst()
{
    if(no_vecs)
        delete [] vec_list;
    if(ph_no_vecs)
        delete [] ph_vec_list;
    if(update_every_ph_no_vecs)
        delete [] update_every_ph_vec_list;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::add_update(int vec_no)
{
    int chk=0;
    for(int i=0;i<ph_no_vecs;i++)
        if(ph_vec_list[i]==vec_no)
            chk=1;
    
    if(chk==0)
        error->abort("update vector should "
        "be of phantom kind");
    
    
    for(int i=0;i<update_every_ph_no_vecs;i++)
        if(update_every_ph_vec_list[i]==vec_no)
            error->abort("duplicate update vector");
    
    
    GROW(update_every_ph_vec_list
    ,update_every_ph_no_vecs
    ,update_every_ph_no_vecs+1);
    update_every_ph_vec_list[update_every_ph_no_vecs]=vec_no;
    update_every_ph_no_vecs++;
    
    
    int tmp;
    for(int i=0;i<update_every_ph_no_vecs;i++)
    {
        for(int j=i+1;j<update_every_ph_no_vecs;j++)
        {
            if(update_every_ph_vec_list[i]==update_every_ph_vec_list[j])
                error->abort("Error: in vector_list two "
                "vectors with the same index cannot exist\n");
            
            if(update_every_ph_vec_list[i]>update_every_ph_vec_list[j])
            {
                tmp=update_every_ph_vec_list[i];
                update_every_ph_vec_list[i]=update_every_ph_vec_list[j];
                update_every_ph_vec_list[j]=tmp;
            }
        }
    }
    
    update_every_ph_byte_size+=atoms->vectors[vec_no].byte_size;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void VecLst::del_update(int vec_no)
{
    int chk=0;
    for(int i=0;i<ph_no_vecs;i++)
        if(ph_vec_list[i]==vec_no)
            chk=1;
    
    if(chk==0)
        error->abort("update vector not found");

    int* vec_tmp;
    CREATE1D(vec_tmp,update_every_ph_no_vecs-1);
    
    int k=0;
    for(int i=0;i<update_every_ph_no_vecs;i++)
    {
        if(ph_vec_list[i]!=vec_no)
        {
            vec_tmp[k]=ph_vec_list[i];
            k++;
        }
    }
    
    delete [] ph_vec_list;
    
    ph_vec_list=vec_tmp;
    
    update_every_ph_byte_size-=atoms->vectors[vec_no].byte_size;
    update_every_ph_no_vecs--;
}
/*----------------------------------------------------------------------------------------*/

/*--------------------------------------------
 
 --------------------------------------------*/
AVec::~AVec()
{
    if(initiated)
    {
        if(type==0)
            delete [] vec_0;
        else if(type==1)
            delete [] vec_1;
        else if(type==2)
            delete [] vec_2;
        else if(type==3)
            delete [] vec_3;
        else if(type==4)
            delete [] vec_4;
        else if(type==5)
            delete [] vec_5;
        else if(type==6)
            delete [] vec_6;
        else if(type==7)
            delete [] vec_7;
        else if(type==8)
            delete [] vec_8;
        else if(type==9)
            delete [] vec_9;
        else if(type==10)
            delete [] vec_10;
        else if(type==11)
            delete [] vec_11;
        else if(type==12)
            delete [] vec_12;
        delete [] name;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <typename TYPE>
TYPE* AVec::create(TYPE*& array,int d0)
{
    array = NULL;
    try
    {
        array = new TYPE [d0];
    }
    catch(bad_alloc&)
    {
        printf("Error: memory allocation failure\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    return array;
}
/*--------------------------------------------
 
 --------------------------------------------*/
template <typename TYPE>
TYPE* AVec::grow(TYPE*& array,int oldsize,int newsize)
{
    if (oldsize==0)
    {
        return create(array,newsize);
    }
    else if (oldsize==newsize)
    {
        return array;
    }
    else
    {
        
        TYPE* newarray=array;
        try
        {
            int size1=newsize;
            int size2=oldsize;
            int size=MIN(size1,size2);
            newarray = new TYPE[newsize];
            memcpy(newarray,array,size*sizeof(TYPE));
            delete [] array;
            array=newarray;
        }
        catch (bad_alloc&)
        {
            printf("Error: memory allocation failure");
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        return array;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::init(Atoms* at,int tp,int p,int d
                ,const char* vec_name)
{
    atms=at;
    type=tp;
    ph=p;
    dim=d;
    
    if(type<0 || type>12)
    {
        printf("Error: Wrong type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    int max_size;
    
    if(ph)
        max_size=dim*atms->atm_vec_ph_size;
    else
        max_size=dim*atms->atm_vec_size;
    
    if(type==0)
    {
        create(vec_0,max_size);
        byte_size=sizeof(char)*dim;
    }
    else if(type==1)
    {
        create(vec_1,max_size);
        byte_size=sizeof(unsigned char)*dim;
    }
    else if(type==2)
    {
        create(vec_2,max_size);
        byte_size=sizeof(short int)*dim;
    }
    else if(type==3)
    {
        create(vec_3,max_size);
        byte_size=sizeof(unsigned short int)*dim;
    }
    else if(type==4)
    {
        create(vec_4,max_size);
        byte_size=sizeof(int)*dim;
    }
    else if(type==5)
    {
        create(vec_5,max_size);
        byte_size=sizeof(unsigned int)*dim;
    }
    else if(type==6)
    {
        create(vec_6,max_size);
        byte_size=sizeof(long int)*dim;
    }
    else if(type==7)
    {
        create(vec_7,max_size);
        byte_size=sizeof(unsigned long int)*dim;
    }
    else if(type==8)
    {
        create(vec_8,max_size);
        byte_size=sizeof(long long int)*dim;
    }
    else if(type==9)
    {
        create(vec_9,max_size);
        byte_size=sizeof(unsigned long long int)*dim;
    }
    else if(type==10)
    {
        create(vec_10,max_size);
        byte_size=sizeof(float)*dim;
    }
    else if(type==11)
    {
        create(vec_11,max_size);
        byte_size=sizeof(double)*dim;
    }
    else if(type==12)
    {
        create(vec_12,max_size);
        byte_size=sizeof(long double)*dim;
    }
    
    int lngth=static_cast<int>(strlen(vec_name))+1;
    create(name,lngth);
    
    memcpy(name,vec_name,lngth*sizeof(char));
    
    
    initiated=1;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int AVec::unpack(char*& buff,int buff_pos)
{
    int no=atms->natms;
    
    if(type==0)
        memcpy(&vec_0[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==1)
        memcpy(&vec_1[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==2)
        memcpy(&vec_2[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==3)
        memcpy(&vec_3[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==4)
        memcpy(&vec_4[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==5)
        memcpy(&vec_5[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==6)
        memcpy(&vec_6[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==7)
        memcpy(&vec_7[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==8)
        memcpy(&vec_8[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==9)
        memcpy(&vec_9[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==10)
        memcpy(&vec_10[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==11)
        memcpy(&vec_11[no*dim],&(buff[buff_pos]),byte_size);
    else if(type==12)
        memcpy(&vec_12[no*dim],&(buff[buff_pos]),byte_size);
    
    buff_pos+=byte_size;
    return buff_pos;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int AVec::pack(char*& buff,int buff_pos,int i)
{
    if(type==0)
        memcpy(&(buff[buff_pos]),&vec_0[i*dim],byte_size);
    else if(type==1)
        memcpy(&(buff[buff_pos]),&vec_1[i*dim],byte_size);
    else if(type==2)
        memcpy(&(buff[buff_pos]),&vec_2[i*dim],byte_size);
    else if(type==3)
        memcpy(&(buff[buff_pos]),&vec_3[i*dim],byte_size);
    else if(type==4)
        memcpy(&(buff[buff_pos]),&vec_4[i*dim],byte_size);
    else if(type==5)
        memcpy(&(buff[buff_pos]),&vec_5[i*dim],byte_size);
    else if(type==6)
        memcpy(&(buff[buff_pos]),&vec_6[i*dim],byte_size);
    else if(type==7)
        memcpy(&(buff[buff_pos]),&vec_7[i*dim],byte_size);
    else if(type==8)
        memcpy(&(buff[buff_pos]),&vec_8[i*dim],byte_size);
    else if(type==9)
        memcpy(&(buff[buff_pos]),&vec_9[i*dim],byte_size);
    else if(type==10)
        memcpy(&(buff[buff_pos]),&vec_10[i*dim],byte_size);
    else if(type==11)
        memcpy(&(buff[buff_pos]),&vec_11[i*dim],byte_size);
    else if(type==12)
        memcpy(&(buff[buff_pos]),&vec_12[i*dim],byte_size);
    
    buff_pos+=byte_size;
    return buff_pos;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int AVec::x_pack(char*& buff,int buff_pos
                 ,int* atm_list,int atm_list_size)
{
    
    if(type==0)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_0[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==1)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_1[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==2)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_2[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==3)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_3[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==4)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_4[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==5)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_5[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==6)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_6[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==7)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_7[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==8)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_8[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==9)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_9[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==10)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_10[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==11)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_11[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    else if(type==12)
    {
        for(int i=0;i<atm_list_size;i++)
        {
            memcpy(&(buff[buff_pos]),&vec_12[atm_list[i]*dim],byte_size);
            buff_pos+=byte_size;
        }
    }
    

    return buff_pos;
}
/*--------------------------------------------
 
 --------------------------------------------*/
int AVec::x_unpack(char*& buff,int buff_pos
                   ,int atm_list_size)
{
    int no=atms->natms+atms->natms_ph;
    int tot_byte=byte_size*atm_list_size;
    
    if(type==0)
        memcpy(&vec_0[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==1)
        memcpy(&vec_1[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==2)
        memcpy(&vec_2[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==3)
        memcpy(&vec_3[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==4)
        memcpy(&vec_4[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==5)
        memcpy(&vec_5[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==6)
        memcpy(&vec_6[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==7)
        memcpy(&vec_7[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==8)
        memcpy(&vec_8[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==9)
        memcpy(&vec_9[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==10)
        memcpy(&vec_10[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==11)
        memcpy(&vec_11[no*dim],&(buff[buff_pos]),tot_byte);
    else if(type==12)
        memcpy(&vec_12[no*dim],&(buff[buff_pos]),tot_byte);
    
    buff_pos+=tot_byte;
    return buff_pos;
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::copy(int dst,int src)
{
    
    src*=dim;
    dst*=dim;
    
    if(type==0)
        memcpy(&vec_0[dst],&vec_0[src],byte_size);
    else if(type==1)
        memcpy(&vec_1[dst],&vec_1[src],byte_size);
    else if(type==2)
        memcpy(&vec_2[dst],&vec_2[src],byte_size);
    else if(type==3)
        memcpy(&vec_3[dst],&vec_3[src],byte_size);
    else if(type==4)
        memcpy(&vec_4[dst],&vec_4[src],byte_size);
    else if(type==5)
        memcpy(&vec_5[dst],&vec_5[src],byte_size);
    else if(type==6)
        memcpy(&vec_6[dst],&vec_6[src],byte_size);
    else if(type==1)
        memcpy(&vec_7[dst],&vec_7[src],byte_size);
    else if(type==8)
        memcpy(&vec_8[dst],&vec_8[src],byte_size);
    else if(type==9)
        memcpy(&vec_9[dst],&vec_9[src],byte_size);
    else if(type==10)
        memcpy(&vec_10[dst],&vec_10[src],byte_size);
    else if(type==11)
        memcpy(&vec_11[dst],&vec_11[src],byte_size);
    else if(type==12)
        memcpy(&vec_12[dst],&vec_12[src],byte_size);
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::grow(int no_atms)
{
    int old_size;
    int new_size;
    if(ph)
    {
        old_size=dim*(atms->atm_vec_ph_size);
        new_size=dim*(atms->atm_vec_ph_size+no_atms);
    }
    else
    {
        old_size=dim*(atms->atm_vec_size);
        new_size=dim*(atms->atm_vec_size+no_atms);
    }
    
    if(type==0)
        grow(vec_0,old_size,new_size);
    else if(type==1)
        grow(vec_1,old_size,new_size);
    else if(type==2)
        grow(vec_2,old_size,new_size);
    else if(type==3)
        grow(vec_3,old_size,new_size);
    else if(type==4)
        grow(vec_4,old_size,new_size);
    else if(type==5)
        grow(vec_5,old_size,new_size);
    else if(type==6)
        grow(vec_6,old_size,new_size);
    else if(type==7)
        grow(vec_7,old_size,new_size);
    else if(type==8)
        grow(vec_8,old_size,new_size);
    else if(type==9)
        grow(vec_9,old_size,new_size);
    else if(type==10)
        grow(vec_10,old_size,new_size);
    else if(type==11)
        grow(vec_11,old_size,new_size);
    else if(type==12)
        grow(vec_12,old_size,new_size);
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::resize(int no_atms)
{
    int new_size;
    if(ph)
        new_size=dim*(atms->atm_vec_ph_size+no_atms);
    else
        new_size=dim*(atms->atm_vec_size+no_atms);
    
    if(type==0)
    {
        delete [] vec_0;
        create(vec_0,new_size);
    }
    else if(type==1)
    {
        delete [] vec_1;
        create(vec_1,new_size);
    }
    else if(type==2)
    {
        delete [] vec_2;
        create(vec_2,new_size);
    }
    else if(type==3)
    {
        delete [] vec_3;
        create(vec_3,new_size);
    }
    else if(type==4)
    {
        delete [] vec_4;
        create(vec_4,new_size);
    }
    else if(type==5)
    {
        delete [] vec_5;
        create(vec_5,new_size);
    }
    else if(type==6)
    {
        delete [] vec_6;
        create(vec_6,new_size);
    }
    else if(type==7)
    {
        delete [] vec_7;
        create(vec_7,new_size);
    }
    else if(type==8)
    {
        delete [] vec_8;
        create(vec_8,new_size);
    }
    else if(type==9)
    {
        delete [] vec_9;
        create(vec_9,new_size);
    }
    else if(type==10)
    {
        delete [] vec_10;
        create(vec_10,new_size);
    }
    else if(type==11)
    {
        delete [] vec_11;
        create(vec_11,new_size);
    }
    else if(type==12)
    {
        delete [] vec_12;
        create(vec_12,new_size);
    }
}

/*--------------------------------------------
 
 --------------------------------------------*/
void* AVec::ret_vec(int no)
{
    no*=dim;
    if(type==0)
        return &vec_0[no];
    else if(type==1)
        return &vec_1[no];
    else if(type==2)
        return &vec_2[no];
    else if(type==3)
        return &vec_3[no];
    else if(type==4)
        return &vec_4[no];
    else if(type==5)
        return &vec_5[no];
    else if(type==6)
        return &vec_6[no];
    else if(type==7)
        return &vec_7[no];
    else if(type==8)
        return &vec_8[no];
    else if(type==9)
        return &vec_9[no];
    else if(type==10)
        return &vec_10[no];
    else if(type==11)
        return &vec_11[no];
    else
        return &vec_12[no];
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void* AVec::ret_vec()
{
    if(type==0)
        return vec_0;
    else if(type==1)
        return vec_1;
    else if(type==2)
        return vec_2;
    else if(type==3)
        return vec_3;
    else if(type==4)
        return vec_4;
    else if(type==5)
        return vec_5;
    else if(type==6)
        return vec_6;
    else if(type==7)
        return vec_7;
    else if(type==8)
        return vec_8;
    else if(type==9)
        return vec_9;
    else if(type==10)
        return vec_10;
    else if(type==11)
        return vec_11;
    else
        return vec_12;
    
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::gather_dump()
{
    MPI_Comm comm_world=atms->comm_world;
    int my_p_no=atms->my_p_no;
    int tot_p=atms->tot_p;
    int tot_natms=atms->tot_natms;
    int natms=atms->natms;
    int* rcv_size=NULL;
    
    if(my_p_no!=0)
        MPI_Send(&natms,1,MPI_INT,0,my_p_no,comm_world);
    
    if(my_p_no==0)
    {
        create(rcv_size,tot_p);
        MPI_Status status;
        for(int iproc=1;iproc<tot_p;iproc++)
            MPI_Recv(&rcv_size[iproc],1,MPI_INT,iproc,iproc,comm_world,&status);
    }
    
    

    if(type==0)
    {
        if(my_p_no==0)
        {
            create(vec_0_dump,dim*tot_natms);
            memcpy(vec_0_dump,vec_0,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_0,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_0_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==1)
    {
        if(my_p_no==0)
        {
            create(vec_1_dump,dim*tot_natms);
            memcpy(vec_1_dump,vec_1,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_1,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_1_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==2)
    {
        if(my_p_no==0)
        {
            create(vec_2_dump,dim*tot_natms);
            memcpy(vec_2_dump,vec_2,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_2,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_2_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==3)
    {
        if(my_p_no==0)
        {
            create(vec_3_dump,dim*tot_natms);
            memcpy(vec_3_dump,vec_3,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_3,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_3_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==4)
    {
        if(my_p_no==0)
        {
            create(vec_4_dump,dim*tot_natms);
            memcpy(vec_4_dump,vec_4,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_4,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_4_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==5)
    {
        if(my_p_no==0)
        {
            create(vec_5_dump,dim*tot_natms);
            memcpy(vec_5_dump,vec_5,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_5,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_5_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==6)
    {
        if(my_p_no==0)
        {
            create(vec_6_dump,dim*tot_natms);
            memcpy(vec_6_dump,vec_6,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_6,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_6_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==7)
    {
        if(my_p_no==0)
        {
            create(vec_7_dump,dim*tot_natms);
            memcpy(vec_7_dump,vec_7,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_7,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_7_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==8)
    {
        if(my_p_no==0)
        {
            create(vec_8_dump,dim*tot_natms);
            memcpy(vec_8_dump,vec_8,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_8,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_8_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==9)
    {
        if(my_p_no==0)
        {
            create(vec_9_dump,dim*tot_natms);
            memcpy(vec_9_dump,vec_9,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_9,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_9_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==10)
    {
        if(my_p_no==0)
        {
            create(vec_10_dump,dim*tot_natms);
            memcpy(vec_10_dump,vec_10,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_10,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_10_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==11)
    {
        if(my_p_no==0)
        {
            create(vec_11_dump,dim*tot_natms);
            memcpy(vec_11_dump,vec_11,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_11,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_11_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    else if(type==12)
    {
        if(my_p_no==0)
        {
            create(vec_12_dump,dim*tot_natms);
            memcpy(vec_12_dump,vec_12,natms*byte_size);
        }
        
        if(my_p_no!=0)
            MPI_Send(vec_12,natms*byte_size,MPI_BYTE,0,my_p_no,comm_world);
        
        if(my_p_no==0)
        {
            MPI_Status status;
            int tot_atoms=natms;
            for(int iproc=1;iproc<tot_p;iproc++)
            {
                MPI_Recv(&vec_12_dump[tot_atoms*dim],(byte_size*rcv_size[iproc]),MPI_BYTE,iproc,iproc,comm_world,&status);
                tot_atoms+=rcv_size[iproc];
            }
        }
    }
    
    if(my_p_no==0)
    {
        delete [] rcv_size;
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::print_dump(FILE* fp,int iatm)
{
    if(atms->my_p_no!=0)
        return;
    
    int icmp=iatm*dim;
    if(type==0)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%c ",vec_0_dump[icmp+i]);
    }
    else if(type==1)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%c ",vec_1_dump[icmp+i]);
    }
    else if(type==2)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%d ",vec_2_dump[icmp+i]);
    }
    else if(type==3)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%d ",vec_3_dump[icmp+i]);
    }
    else if(type==4)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%d ",vec_4_dump[icmp+i]);
    }
    else if(type==5)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%d ",vec_5_dump[icmp+i]);
    }
    else if(type==6)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%ld ",vec_6_dump[icmp+i]);
    }
    else if(type==7)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%ld ",vec_7_dump[icmp+i]);
    }
    else if(type==8)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%lld ",vec_8_dump[icmp+i]);
    }
    else if(type==9)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%lld ",vec_9_dump[icmp+i]);
    }
    else if(type==10)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%f ",vec_10_dump[icmp+i]);
    }
    else if(type==11)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%lf ",vec_11_dump[icmp+i]);
    }
    else if(type==12)
    {
        for(int i=0;i<dim;i++)
            fprintf(fp,"%Lf ",vec_12_dump[icmp+i]);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::ret(char*& x)
{
    if(type==0)
        x=vec_0;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(unsigned char*& x)
{
    if(type==1)
        x=vec_1;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(short int*& x)
{
    if(type==2)
        x=vec_2;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(unsigned short int*& x)
{
    if(type==3)
        x=vec_3;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
}
void AVec::ret(int*& x)
{
    if(type==4)
        x=vec_4;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(unsigned int*& x)
{
    if(type==5)
        x=vec_5;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(long int*& x)
{
    if(type==6)
        x=vec_6;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(unsigned long int*& x)
{
    if(type==7)
        x=vec_7;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(long long int*& x)
{
    if(type==8)
        x=vec_8;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(unsigned long long int*& x)
{
    if(type==9)
        x=vec_9;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(float*& x)
{
    if(type==10)
        x=vec_10;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(double*& x)
{
    if(type==11)
        x=vec_11;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret(long double*& x)
{
    if(type==12)
        x=vec_12;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::ret_dump(char*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==0)
        x=vec_0_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(unsigned char*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==1)
        x=vec_1_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(short int*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==2)
        x=vec_2_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(unsigned short int*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==3)
        x=vec_3_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
}
void AVec::ret_dump(int*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==4)
        x=vec_4_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(unsigned int*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==5)
        x=vec_5_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(long int*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==6)
        x=vec_6_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(unsigned long int*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==7)
        x=vec_7_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(long long int*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==8)
        x=vec_8_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(unsigned long long int*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==9)
        x=vec_9_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(float*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==10)
        x=vec_10_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(double*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==11)
        x=vec_11_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
void AVec::ret_dump(long double*& x)
{
    if(atms->my_p_no!=0)
        return;
    if(type==12)
        x=vec_12_dump;
    else
    {
        printf("Error: incorrect type");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
}
/*--------------------------------------------
 
 --------------------------------------------*/
void AVec::del_dump()
{
    if(atms->my_p_no!=0)
        return;
    if(type==0)
        delete [] vec_0_dump;
    else if(type==1)
        delete [] vec_1_dump;
    else if(type==2)
        delete [] vec_2_dump;
    else if(type==3)
        delete [] vec_3_dump;
    else if(type==4)
        delete [] vec_4_dump;
    else if(type==5)
        delete [] vec_5_dump;
    else if(type==6)
        delete [] vec_6_dump;
    else if(type==7)
        delete [] vec_7_dump;
    else if(type==8)
        delete [] vec_8_dump;
    else if(type==9)
        delete [] vec_9_dump;
    else if(type==10)
        delete [] vec_10_dump;
    else if(type==11)
        delete [] vec_11_dump;
    else if(type==12)
        delete [] vec_12_dump;
    
}
/*--------------------------------------------

 --------------------------------------------*/
void AVec::change_dimension(int d)
{
    
    
    if(d==dim)
        return;
    
    int d_min=MIN(d,dim);
    
    int tot;
    if(ph==0)
        tot=atms->natms;
    else
        tot=atms->natms+atms->natms_ph;
    
    if(type==0)
    {
        char* tmp_vec_0;
        create(tmp_vec_0,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_0[i*d+j]=vec_0[i*dim+j];
        
        delete [] vec_0;
        vec_0=tmp_vec_0;
        byte_size=sizeof(char)*d;
    }
    else if(type==1)
    {
        unsigned char* tmp_vec_1;
        create(tmp_vec_1,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_1[i*d+j]=vec_1[i*dim+j];
        
        delete [] vec_1;
        vec_1=tmp_vec_1;
        byte_size=sizeof(unsigned char)*d;
    }
    else if(type==2)
    {
        short int* tmp_vec_2;
        create(tmp_vec_2,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_2[i*d+j]=vec_2[i*dim+j];
        
        delete [] vec_2;
        vec_2=tmp_vec_2;
        byte_size=sizeof(short int)*d;
    }
    else if(type==3)
    {
        unsigned short int* tmp_vec_3;
        create(tmp_vec_3,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_3[i*d+j]=vec_3[i*dim+j];
        
        delete [] vec_3;
        vec_3=tmp_vec_3;
        byte_size=sizeof(unsigned short int)*d;
    }
    else if(type==4)
    {
        int* tmp_vec_4;
        create(tmp_vec_4,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_4[i*d+j]=vec_4[i*dim+j];
        
        delete [] vec_4;
        vec_4=tmp_vec_4;
        byte_size=sizeof(int)*d;
    }
    else if(type==5)
    {
        unsigned int* tmp_vec_5;
        create(tmp_vec_5,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_5[i*d+j]=vec_5[i*dim+j];
        
        delete [] vec_5;
        vec_5=tmp_vec_5;
        byte_size=sizeof(unsigned int)*d;
    }
    else if(type==6)
    {
        long int* tmp_vec_6;
        create(tmp_vec_6,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_6[i*d+j]=vec_6[i*dim+j];
        
        delete [] vec_6;
        vec_6=tmp_vec_6;
        byte_size=sizeof(long int)*d;
    }
    else if(type==7)
    {
        unsigned long int* tmp_vec_7;
        create(tmp_vec_7,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_7[i*d+j]=vec_7[i*dim+j];
        
        delete [] vec_7;
        vec_7=tmp_vec_7;
        byte_size=sizeof(unsigned long int)*d;
    }
    else if(type==8)
    {
        long long int* tmp_vec_8;
        create(tmp_vec_8,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_8[i*d+j]=vec_8[i*dim+j];
        
        delete [] vec_8;
        vec_8=tmp_vec_8;
        byte_size=sizeof(long long int)*d;
    }
    else if(type==9)
    {
        unsigned long long int* tmp_vec_9;
        create(tmp_vec_9,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_9[i*d+j]=vec_9[i*dim+j];
        
        delete [] vec_9;
        vec_9=tmp_vec_9;
        byte_size=sizeof(unsigned long long int)*d;
    }
    else if(type==10)
    {
        float* tmp_vec_10;
        create(tmp_vec_10,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_10[i*d+j]=vec_10[i*dim+j];
        
        delete [] vec_10;
        vec_10=tmp_vec_10;
        byte_size=sizeof(float)*d;
    }
    else if(type==11)
    {
        double* tmp_vec_11;
        create(tmp_vec_11,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_11[i*d+j]=vec_11[i*dim+j];
        
        delete [] vec_11;
        vec_11=tmp_vec_11;
        byte_size=sizeof(double)*d;
    }
    else if(type==12)
    {
        long double* tmp_vec_12;
        create(tmp_vec_12,d*tot);
        for(int i=0;i<tot;i++)
            for(int j=0;j<d_min;j++)
                tmp_vec_12[i*d+j]=vec_12[i*dim+j];
        
        delete [] vec_12;
        vec_12=tmp_vec_12;
        byte_size=sizeof(long double)*d;
    }
    dim=d;
    
}

