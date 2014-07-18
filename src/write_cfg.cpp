#include "atom_types.h"
#include "write_cfg.h"
#include "error.h"
#include "memory.h"
using namespace std;
using namespace MAPP_NS;

/*--------------------------------------------
 constructor
 write cfg 1000 dump x
 --------------------------------------------*/
WriteCFG::WriteCFG(MAPP* mapp,int narg
,char** arg) : Write(mapp)
{
    sorting=0;
    if(narg<=4)
        error->abort("wrong write cfg command");

    write_step_tally=atoi(arg[2]);
    int lngth= static_cast<int>(strlen(arg[3]))+1;
    CREATE1D(file_name,lngth);
    for(int i=0;i<lngth;i++)
        file_name[i]=arg[3][i];
    

    id_xst=0;
    id_cmp=-1;
    int iarg=4;
    no_vecs=0;
    int tmp;
    while(iarg<narg)
    {
        if(strcmp(arg[iarg],"sort")==0)
        {
            sorting=1;
            iarg++;
        }
        else
        {
            if(strcmp(arg[iarg],"id")==0 ||
               strcmp(arg[iarg],"type")==0)
                error->abort("%s cannot be included in auxilary",arg[iarg]);
            
            tmp=atoms->find_exist(arg[iarg]);
            if(tmp<0)
                error->abort("cannot find vector %s",arg[iarg]);
            GROW(vec_list,no_vecs,no_vecs+1);
            vec_list[no_vecs]=tmp;
            no_vecs++;
            iarg++;
        }
    }
    
    
    
    for(int i=0;i<no_vecs;i++)
    {
        for(int j=i+1;j<no_vecs;j++)
        {
            if(vec_list[i]>vec_list[j])
            {
                if(vec_list[i]==vec_list[j])
                    error->abort("two components cannot be the same");
                
                tmp=vec_list[i];
                vec_list[i]=vec_list[j];
                vec_list[j]=tmp;
            }
        }
    }
    
    tot_dim=0;
    for(int i=0;i<no_vecs;i++)
        tot_dim+=atoms->vectors[vec_list[i]].dim;
    
    if(vec_list[0]!=0)
        error->abort("vector x should be included");
    
    type_cmp=atoms->find("type");
    GROW(vec_list,no_vecs,no_vecs+1);
    vec_list[no_vecs]=type_cmp;
    no_vecs++;
    
    if(sorting)
    {
        id_cmp=atoms->find("id");
        GROW(vec_list,no_vecs,no_vecs+1);
        vec_list[no_vecs]=id_cmp;
        no_vecs++;
    }
    
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
WriteCFG::~WriteCFG()
{
    delete [] file_name;
    
    if(no_vecs)
        delete [] vec_list;
    
}
/*--------------------------------------------
 write file
 --------------------------------------------*/
void WriteCFG::write_file(int stp)
{

    atoms->x2s_no_correction(atoms->natms);
    
    for(int i=0;i<no_vecs;i++)
        atoms->vectors[vec_list[i]].gather_dump();
    
    
    FILE* fp=NULL;
    
    if(atoms->my_p_no==0)
    {
        char* filename;
        CREATE1D(filename,MAXCHAR);
        sprintf (filename, "%s.%08d.cfg",file_name,stp);
        fp=fopen(filename,"w");
        delete [] filename;
        
        // write the header
        fprintf(fp,"Number of particles = %d\n",atoms->tot_natms);
        fprintf(fp,"A = 1 Angstrom (basic length-scale)\n");
        
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                fprintf(fp,"H0(%d,%d) = %lf A\n",i+1,j+1,atoms->H[i][j]);
        
        fprintf(fp,".NO_VELOCITY.\n");

        fprintf(fp,"entry_count = %d\n",tot_dim);
        
        
        // write the body
        if(sorting)
        {
            int icomp=0;
            for(int i=1;i<no_vecs-2;i++)
            {
                if(atoms->vectors[vec_list[i]].dim>1)
                {
                    for(int idim=0;idim<atoms->vectors[vec_list[i]].dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                                ,icomp++,atoms->vectors[vec_list[i]].name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,atoms->vectors[vec_list[i]].name);
                }
                
            }
            
            int* id;
            int tot_natms=atoms->tot_natms;
            atoms->vectors[id_cmp].ret_dump(id);
            int* sort;
            CREATE1D(sort,tot_natms);
            for(int i=0;i<tot_natms;i++)
                sort[id[i]]=i;
            
            int* type;
            atoms->vectors[type_cmp].ret_dump(type);
            
            for(int itype=0;itype<atom_types->no_types;itype++)
            {
                fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                       ,atom_types->atom_names[itype]);
                for(int i=0;i<tot_natms;i++)
                {
                    if(type[sort[i]]==itype)
                    {
                        for(int j=0;j<no_vecs-2;j++)
                            atoms->vectors[vec_list[j]].print_dump(fp,sort[i]);
                        fprintf(fp,"\n");
                    }
                }
            }
            
            if(tot_natms)
                delete [] sort;
        }
        else
        {
            int icomp=0;
            for(int i=1;i<no_vecs-1;i++)
            {
                if(atoms->vectors[vec_list[i]].dim>1)
                {
                    for(int idim=0;idim<atoms->vectors[vec_list[i]].dim;idim++)
                        fprintf(fp,"auxiliary[%d] = %s_%d [reduced unit]\n"
                                ,icomp++,atoms->vectors[vec_list[i]].name,idim);
                }
                else
                {
                    fprintf(fp,"auxiliary[%d] = %s [reduced unit]\n"
                            ,icomp++,atoms->vectors[vec_list[i]].name);
                }
                
            }
            int tot_natms=atoms->tot_natms;
            
            int* type;
            atoms->vectors[type_cmp].ret_dump(type);
            
            for(int itype=0;itype<atom_types->no_types;itype++)
            {
                fprintf(fp,"%lf \n%s \n",atom_types->mass[itype]
                       ,atom_types->atom_names[itype]);
                for(int i=0;i<tot_natms;i++)
                {
                    if(type[i]==itype)
                    {
                        for(int j=0;j<no_vecs-1;j++)
                            atoms->vectors[vec_list[j]].print_dump(fp,i);
                        fprintf(fp,"\n");
                    }
                }
            }
        }
        
        fclose(fp);
    }
    
    for(int i=0;i<no_vecs;i++)
        atoms->vectors[vec_list[i]].del_dump();
    
    atoms->s2x(atoms->natms);
}

