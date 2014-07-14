/*--------------------------------------------
 Created by Sina on 06/05/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "neighbor.h"
#include "FF.h"
using namespace std;
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Neighbor::Neighbor(MAPP* mapp):InitPtrs(mapp)
{
    neighbor_list_size_size=0;
    tot_bin=0;
    first_atom_bin_size=0;
    next_atm_size=0;
    pair_wise=1;
    atm_bin_size=0;
    /*
    if(atoms->dimension)
    {
        int d=atoms->dimension;
        CREATE1D(tot_bin_grid,d);
        CREATE1D(bin_size,d);
        CREATE1D(bin_denom_list,d);
        CREATE1D(s_tmp,d);
    }*/
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Neighbor::~Neighbor()
{
    /*
    if(atoms->dimension)
    {
        delete [] tot_bin_grid;
        delete [] bin_size;
        delete [] bin_denom_list;
        delete [] s_tmp;
    }
     */
    
    if(tot_bin)
    {
        for(int i=0;i<tot_bin;i++)
            if(bin_neigh_list_size[i])
                delete [] bin_neigh_list[i];
        
        delete [] bin_neigh_list;
        delete [] bin_neigh_list_size;
    }
    
    if(neighbor_list_size_size)
    {
        for(int i=0;i<neighbor_list_size_size;i++)
            if(neighbor_list_size[i])
                delete [] neighbor_list[i];
        delete [] neighbor_list_size;
    }
    
    if(first_atom_bin_size)
        delete [] first_atom_bin;
    if(next_atm_size)
        delete [] next_atm;
    if(atm_bin_size)
        delete [] atm_bin;
    
}
/*--------------------------------------------
 initiation before MD
 --------------------------------------------*/
void Neighbor::init()
{
       
    if(atoms->dimension)
    {
        int d=atoms->dimension;
        CREATE1D(tot_bin_grid,d);
        CREATE1D(bin_size,d);
        CREATE1D(bin_denom_list,d);
        CREATE1D(s_tmp,d);
    }

    type_n=atoms->find("type");
    create_bin_list();
}
/*--------------------------------------------
 finalize MD
 --------------------------------------------*/
void Neighbor::fin()
{
    if(atoms->dimension)
    {
        delete [] tot_bin_grid;
        delete [] bin_size;
        delete [] bin_denom_list;
        delete [] s_tmp;
    }

    if(tot_bin)
    {
        for(int i=0;i<tot_bin;i++)
            if(bin_neigh_list_size[i])
                delete [] bin_neigh_list[i];
        
        delete [] bin_neigh_list;
        delete [] bin_neigh_list_size;
    }
    
    if(neighbor_list_size_size)
    {
        for(int i=0;i<neighbor_list_size_size;i++)
            if(neighbor_list_size[i])
                delete [] neighbor_list[i];
        delete [] neighbor_list_size;
    }
    if(first_atom_bin_size)
        delete [] first_atom_bin;
    if(next_atm_size)
        delete [] next_atm;
    if(atm_bin_size)
        delete [] atm_bin;
    
    neighbor_list_size_size=0;
    tot_bin=0;
    first_atom_bin_size=0;
    next_atm_size=0;
    atm_bin_size=0;
    
}
/*--------------------------------------------
 create the neighbopr list
 --------------------------------------------*/
void Neighbor::create_list(int box_change,int s_or_x)
{
    if(box_change)
        create_bin_list();

    if(s_or_x)
        bin_atoms_s();
    else
        bin_atoms();
 
    
    if(neighbor_list_size_size)
    {
        for(int i=0;i<neighbor_list_size_size;i++)
            if(neighbor_list_size[i])
                delete [] neighbor_list[i];

    
        delete [] neighbor_list_size;
        delete [] neighbor_list;

    }
    
    neighbor_list_size_size=atoms->natms;
   
    neighbor_list=CREATE1D(neighbor_list,neighbor_list_size_size);
    CREATE1D(neighbor_list_size,neighbor_list_size_size);
    for(int i=0;i<neighbor_list_size_size;i++)
        neighbor_list_size[i]=0;
    
    
    //int* type=(int*)atoms->vectors[type_no].ret_vec();
    //TYPE0* x=(TYPE0*)atoms->vectors[0].ret_vec();
    
    TYPE0* x;
    atoms->vectors[0].ret(x);
    int* type;
    atoms->vectors[type_n].ret(type);
    int x_dim=atoms->vectors[0].dim;

    int dim=atoms->dimension;
    TYPE0 cut_sq,rsq;
    int icomp,jcomp;
    int ibin,jbin;
    int iatm,jatm;
    TYPE0* cut_sk_sq=forcefield->cut_sk_sq;
    
    int* tmp_neigh_list;
    int tmp_neigh_list_size=1024;
    int tmp_neigh_list_grow=50;
    CREATE1D(tmp_neigh_list,tmp_neigh_list_size);
    /*
    for(iatm=0;iatm<atoms->natms;iatm++)
    {
        icomp=x_dim*iatm;
        for(jatm=iatm+1;jatm<atoms->natms+atoms->natms_ph;jatm++)
        {
            jcomp=x_dim*jatm;
            rsq=0.0;
            for(int idim=0;idim<dim;idim++)
                rsq+=(x[icomp+idim]-x[jcomp+idim])
                *(x[icomp+idim]-x[jcomp+idim]);
            cut_sq=cut_sk_sq[COMP(type[iatm],type[jatm])];
            
            if(rsq<cut_sq)
            {
                if(neighbor_list_size[iatm]+1>tmp_neigh_list_size)
                {
                    GROW(tmp_neigh_list, tmp_neigh_list_size,tmp_neigh_list_size+tmp_neigh_list_grow);
                    tmp_neigh_list_size+=tmp_neigh_list_grow;
                }
                tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                neighbor_list_size[iatm]++;
            }
            
        }
        if(neighbor_list_size[iatm])
        {
            CREATE1D(neighbor_list[iatm],neighbor_list_size[iatm]);
            memcpy(neighbor_list[iatm],tmp_neigh_list,neighbor_list_size[iatm]*sizeof(int));
        }
    }
    if(atm_bin_size)
        delete [] atm_bin;
    atm_bin_size=0;
    delete [] tmp_neigh_list;
    return;
    
    */
  
    if(pair_wise)
    {
        for(iatm=0;iatm<atoms->natms;iatm++)
        {
            icomp=x_dim*iatm;
            ibin=atm_bin[iatm];
            
            for(int j=0;j<bin_neigh_list_size[ibin];j++)
            {
                
                jbin=bin_neigh_list[ibin][j];
                jatm=first_atom_bin[jbin];
                while(jatm!=-1)
                {
                    if(jatm>iatm)
                    {
                        cut_sq=cut_sk_sq[COMP(type[iatm],type[jatm])];
                        jcomp=x_dim*jatm;
                        
                        rsq=0.0;
                        for(int idim=0;idim<dim;idim++)
                            rsq+=(x[icomp+idim]-x[jcomp+idim])
                            *(x[icomp+idim]-x[jcomp+idim]);
                        
                        if(rsq<cut_sq)
                        {
                            if(neighbor_list_size[iatm]+1>tmp_neigh_list_size)
                            {
                                GROW(tmp_neigh_list, tmp_neigh_list_size,tmp_neigh_list_size+tmp_neigh_list_grow);
                                tmp_neigh_list_size+=tmp_neigh_list_grow;
                            }
                            tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                            neighbor_list_size[iatm]++;
                        }
                    }
                    jatm=next_atm[jatm];
                }
            }
            if(neighbor_list_size[iatm])
            {
                CREATE1D(neighbor_list[iatm],neighbor_list_size[iatm]);
                memcpy(neighbor_list[iatm],tmp_neigh_list,neighbor_list_size[iatm]*sizeof(int));
            }
        }
    }
    else
    {
        for(iatm=0;iatm<atoms->natms;iatm++)
        {
            icomp=x_dim*iatm;
            ibin=atm_bin[iatm];
            
            for(int j=0;j<bin_neigh_list_size[ibin];j++)
            {
                jbin=bin_neigh_list[ibin][j];
                jatm=first_atom_bin[jbin];
                while(jatm!=-1)
                {
                    if(jatm!=iatm)
                    {
                        cut_sq=cut_sk_sq[COMP(type[iatm],type[jatm])];
                        jcomp=x_dim*jatm;
                        
                        rsq=0.0;
                        for(int idim=0;idim<dim;idim++)
                            rsq+=(x[icomp+idim]-x[jcomp+idim])
                            *(x[icomp+idim]-x[jcomp+idim]);
                        
                        if(rsq<cut_sq)
                        {
                            if(neighbor_list_size[iatm]+1>tmp_neigh_list_size)
                            {
                                GROW(tmp_neigh_list, tmp_neigh_list_size,tmp_neigh_list_size+tmp_neigh_list_grow);
                                tmp_neigh_list_size+=tmp_neigh_list_grow;
                            }
                            tmp_neigh_list[neighbor_list_size[iatm]]=jatm;
                            neighbor_list_size[iatm]++;
                        }
                    }
                    
                    jatm=next_atm[jatm];
                }
                
            }
            if(neighbor_list_size[iatm])
            {
                CREATE1D(neighbor_list[iatm],neighbor_list_size[iatm]);
                memcpy(neighbor_list[iatm],tmp_neigh_list,neighbor_list_size[iatm]*sizeof(int));
            }
        }
    }
    delete [] tmp_neigh_list;
    if(atm_bin_size)
        delete [] atm_bin;
    atm_bin_size=0;
    
    
}
/*--------------------------------------------
 create the bin neighbor list
 --------------------------------------------*/
void Neighbor::create_bin_list()
{
    TYPE0* cut_ph_s=atoms->cut_ph_s;
    TYPE0* s_lo=atoms->s_lo;
    TYPE0* s_hi=atoms->s_hi;
    int dim=atoms->dimension;
    if(tot_bin)
    {
        for(int i=0;i<tot_bin;i++)
            if(bin_neigh_list_size[i])
                delete [] bin_neigh_list[i];
        
        delete [] bin_neigh_list;
        delete [] bin_neigh_list_size;
    }
    tot_bin=1;
    for (int i=0;i<dim;i++)
    {
        bin_size[i]=cut_ph_s[i];
        tot_bin_grid[i]=static_cast<int>
        ((2.0*cut_ph_s[i]+s_hi[i]-s_lo[i])/bin_size[i])+1;
        tot_bin*=tot_bin_grid[i];
    }
    
    CREATE1D(bin_neigh_list,tot_bin);
    CREATE1D(bin_neigh_list_size,tot_bin);
    for(int i=0;i<tot_bin;i++)
        bin_neigh_list_size[i]=0;
    for(int i=0;i<dim;i++)
    {
        int no=1;
        for(int j=0;j<i;j++)
            no*=tot_bin_grid[j];
        bin_denom_list[i]=no;
    }
    int* ibin_loc;
    CREATE1D(ibin_loc,dim);
    int* list_ii;
    CREATE1D(list_ii,dim);
    int tmp;
    for(int ibin=0;ibin<tot_bin;ibin++)
    {
        tmp=ibin;
        for(int i=dim-1;i>-1;i--)
        {
            ibin_loc[i]=tmp/bin_denom_list[i];
            tmp-=ibin_loc[i]*bin_denom_list[i];
        }

        for(int i=0;i<dim;i++)
            list_ii[i]=0;
        
        find_bin_no(dim,0,list_ii
        ,ibin,ibin_loc,bin_denom_list);
        
    }
    delete [] ibin_loc;
    delete [] list_ii;

 }
/*--------------------------------------------
 find the bin no based
 --------------------------------------------*/
void Neighbor::find_bin_no(int dim,int pos,int*& list_ii
,int ibin,int* ibin_loc,int* list)
{

    if(dim>0)
    {
        for(int i=0;i<atoms->dimension;i++)
        {
            list_ii[pos]=-1+i;
            find_bin_no(dim-1,pos+1,list_ii,ibin,ibin_loc,list);
        }
    }
    else
    {
        int no=0;
        int tmp;
        int chk=0;
        for(int i=0;i<pos;i++)
        {
            tmp=ibin_loc[i]+list_ii[i];
            if(tmp<0 || tmp>(tot_bin_grid[i]-1))
                chk=1;
            no+=tmp*list[i];
        }
        if(chk==0)
        {
            GROW(bin_neigh_list[ibin]
            ,bin_neigh_list_size[ibin]
            ,bin_neigh_list_size[ibin]+1);
            bin_neigh_list[ibin][bin_neigh_list_size[ibin]]=no;
            bin_neigh_list_size[ibin]++;
        }
    }    
}
/*--------------------------------------------
 bin the atoms
 --------------------------------------------*/
void Neighbor::bin_atoms()
{
    if(first_atom_bin_size)
        delete [] first_atom_bin;

    CREATE1D(first_atom_bin,tot_bin);
    first_atom_bin_size=tot_bin;
    for(int i=0;i<first_atom_bin_size;i++)
        first_atom_bin[i]=-1;
    
    int tot_natms=atoms->natms+atoms->natms_ph;

    if(next_atm_size)
        delete [] next_atm;
    CREATE1D(next_atm,tot_natms);
    next_atm_size=tot_natms;
    
    
    //TYPE0* x=(TYPE0*)atoms->vectors[0].ret_vec();
    TYPE0* x;
    atoms->vectors[0].ret(x);
    
    int x_dim=atoms->vectors[0].dim;
    
    if(atm_bin_size)
        delete [] atm_bin;
    CREATE1D(atm_bin,atoms->natms);
    atm_bin_size=atoms->natms;
    
    int bin;
    for(int i=tot_natms-1;i>-1;i--)
    {
        bin=x2bin(&x[x_dim*i]);
        if(i<atoms->natms)
            atm_bin[i]=bin;
        next_atm[i]=first_atom_bin[bin];
        first_atom_bin[bin]=i;
    }
}
/*--------------------------------------------
 x 2 bin no
 --------------------------------------------*/
int Neighbor::x2bin(TYPE0* x)
{
    TYPE0** B=atoms->B;
    
    int dim=atoms->dimension;
    for(int j=0;j<dim;j++)
    {
        s_tmp[j]=0.0;
        for(int k=j;k<dim;k++)
            s_tmp[j]+=B[k][j]*x[k];
    }
    
    TYPE0* cut_ph_s=atoms->cut_ph_s;
    TYPE0* s_lo=atoms->s_lo;
    
    int no=0;
    for(int i=0;i<dim;i++)
    {
        no+=static_cast<int>
        ((s_tmp[i]+cut_ph_s[i]-s_lo[i])/bin_size[i])
        *bin_denom_list[i];
    }
    
    return no;
}
/*--------------------------------------------
 bin the atoms
 --------------------------------------------*/
void Neighbor::bin_atoms_s()
{
    if(first_atom_bin_size)
        delete [] first_atom_bin;
    CREATE1D(first_atom_bin,tot_bin);
    first_atom_bin_size=tot_bin;
    for(int i=0;i<first_atom_bin_size;i++)
        first_atom_bin[i]=-1;
    
    int tot_natms=atoms->natms+atoms->natms_ph;
    
    if(next_atm_size)
        delete [] next_atm;
    CREATE1D(next_atm,tot_natms);
    next_atm_size=tot_natms;
    
    if(atm_bin_size)
        delete [] atm_bin;
    CREATE1D(atm_bin,atoms->natms);
    atm_bin_size=atoms->natms;
    
    int bin;
    //TYPE0* s=(TYPE0*)atoms->vectors[0].ret_vec();
    TYPE0* s;
    atoms->vectors[0].ret(s);

    int s_dim=atoms->vectors[0].dim;
    for(int i=tot_natms-1;i>-1;i--)
    {
        bin=s2bin(&s[s_dim*i]);
        if(i<atoms->natms)
            atm_bin[i]=bin;
        next_atm[i]=first_atom_bin[bin];
        first_atom_bin[bin]=i;
    }
    atoms->s2x(tot_natms);
}
/*--------------------------------------------
 s 2 bin no
 --------------------------------------------*/
int Neighbor::s2bin(TYPE0* s)
{
    int dim=atoms->dimension;
    
    TYPE0* cut_ph_s=atoms->cut_ph_s;
    TYPE0* s_lo=atoms->s_lo;
    
    int no=0;
    for(int i=0;i<dim;i++)
    {
        no+=static_cast<int>
        ((s[i]+cut_ph_s[i]-s_lo[i])/bin_size[i])
        *bin_denom_list[i];
        
    }
    
    return no;
    
}


