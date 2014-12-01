/*--------------------------------------------
 Created by Sina on 05/11/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "MAPP.h"
#include "memory.h"
#include "error.h"
#include "atoms.h"
#include "neighbor.h"
#include "atom_types.h"
#include "min.h"
#include "min_styles.h"
#include "clock.h"
#include "clock_styles.h"
#include "ff.h"
#include "ff_styles.h"
#include "md.h"
#include "md_styles.h"
#include "read.h"
#include "read_styles.h"
#include "write.h"
#include "write_styles.h"
#include "command_styles.h"
#include <cmath>
#include <stdio.h>
#define MAPP_VERSION "2.0.0"
using namespace MAPP_NS;
using namespace std;
/*--------------------------------------------
 constructor of the main executer
 --------------------------------------------*/
MAPP::
MAPP(int narg,char** args,MPI_Comm communicator)
{

    step_no=0;
    step_tally=1000;
    output = stdout;
    world=communicator;
    
    error = new Error(this);
    memory = new Memory(this);
    atoms = new Atoms(this,world);
    atom_types = new AtomTypes(this);
    neighbor = new Neighbor(this);
    
    //print the version
    if (atoms->my_p_no==0)
        fprintf(output,"MAPP Version: %s\n\n",
            (char*) MAPP_VERSION);
    
    
    forcefield=NULL;
    md=NULL;
    min=NULL;
    clock=NULL;
    write=NULL;
    
    mode=MD;
    atoms->add<TYPE0>(1, 3,"x");
    atoms->add<int>(1, 1,"type");
    atoms->add<int>(0, 1,"id");
    
    input_file=NULL;
    input_file=stdin;
    
    int iarg=1;
    while(iarg<narg)
    {
        
        if(strcmp(args[iarg],"-i")==0)
        {
            iarg++;
            if(iarg==narg)
                error->abort("no input file");
            input_file=fopen(args[iarg],"r");
            iarg++;
        }
        else if(strcmp(args[iarg],"-o")==0)
        {
            iarg++;
            if(iarg==narg)
                error->abort("no output file");
            output=fopen(args[iarg],"w");
            iarg++;
        }
        else
            error->abort("unknown postfix: %s"
            ,args[iarg]);
    }
    
    if(input_file==NULL)
        error->abort("input file not found");
    
    read_file();
    
    if(input_file!=stdin)
        fclose(input_file);
    
    
    
    //test();
    
}
/*--------------------------------------------
 destructor of the main executer
 --------------------------------------------*/
MAPP::~MAPP()
{
    if (atoms->my_p_no==0)
        fprintf(output,"Finito\n");
    
    if(output!=stdout)
        fclose(output);
    
    if(write!=NULL)
        delete write;

    delete neighbor;
    delete atom_types;
    delete atoms;
    delete error;
    delete memory;
}
/*--------------------------------------------
 test field
 --------------------------------------------*/
void MAPP::read_file()
{

    int input_file_chk=1;
    char* line;
    CREATE1D(line,MAXCHAR);
    while (input_file_chk)
    {
        if(atoms->my_p_no==0)
        {
            fgets(line,MAXCHAR,input_file);
            if(feof(input_file))
                input_file_chk=0;
        }
        MPI_Bcast(&input_file_chk,1,MPI_INT,0,world);
        if(input_file_chk==0)
            continue;
        MPI_Bcast(line,MAXCHAR,MPI_CHAR,0,world);
        command(line);
    }
    
    delete [] line;
}
/*--------------------------------------------
 analysing the commands
 --------------------------------------------*/
void MAPP::command(char* command)
{
    char** args;
    int narg;
    narg=parse_line(command,args);
    
    if(narg==0)
    {
        return;
    }
    else if(strcmp(args[0],"skin")==0)
    {
        atoms->add_skin(narg,args);
    }
    else if(strcmp(args[0],"mode")==0)
    {
        change_mode(narg,args);
    }
    else if(strcmp(args[0],"step_tally")==0)
    {
        if(narg!=2)
            error->abort("unknown command: %s",command);
        step_tally=atoi(args[1]);
        if(step_tally<=0)
            error->abort("step tally cannot "
            "be equal or less than zero");
    }
    else if(strcmp(args[0],"reset")==0)
    {
        if(narg!=1)
            error->abort("unknown command: %s",command);
        step_no=0;
    }

    else if(strcmp(args[0],"ff_coef")==0)
    {
        if(forcefield==NULL)
            error->abort("cannot add the coefficients"
            " before the forcefield is initiated");
        forcefield->coef(narg,args);
    }
    else if(strcmp(args[0],"read")==0) read_style(narg,args);
    else if(strcmp(args[0],"ff")==0) ff_style(narg,args);
    else if(strcmp(args[0],"min")==0) min_style(narg,args);
    else if(strcmp(args[0],"clock")==0) clock_style(narg,args);
    else if(strcmp(args[0],"md")==0) md_style(narg,args);
    else if(strcmp(args[0],"write")==0) write_style(narg,args);
    else if(strcmp(args[0],"del_write")==0)
    {
        if(write!=NULL)
            delete write;
        write=NULL;
    }
    else if(strcmp(args[0],"time_step")==0)
    {
        if(md==NULL)
            error->abort("before adjusting the "
            "time_step, ensemble should be initialized");
        md->add_dt(narg,args);
    }
    else if(strcmp(args[0],"boltzmann")==0)
    {
        if(md==NULL)
            error->abort("before adjusting the boltzmann"
            " constant, ensemble should be initialized");
        md->add_boltzmann(narg,args);
    }
    else if(strcmp(args[0],"run")==0)
    {
        if(md==NULL)
            error->abort("before run, ensemble should be initialized");
        md->run(narg,args);
    }
    else
        command_style(narg,args);
    
    
    for(int i=0;i<narg;i++)
        delete [] args[i];
    if(narg)
        delete [] args;
}
/*--------------------------------------------
 differnt forcefield styles
 --------------------------------------------*/
void MAPP::ff_style(int narg,char** args)
{
    if(narg!=2)
        error->abort("wrong command: %s",args[0]);
    
    if(forcefield!=NULL)
        delete forcefield;
    
    #define FF_Style
    #define FFStyle(class_name,style_name)     \
    else if(strcmp(args[1],#style_name)==0)    \
    forcefield=new class_name(this);
    
    //different forcefileds
    if(0){}
    #include "ff_styles.h"
    else
        error->abort("unknown forcefield: %s"
                     ,args[1]);
    
    #undef FF_Style
}
/*--------------------------------------------
 differnt minimization styles
 --------------------------------------------*/
void MAPP::min_style(int narg,char** args)
{
    if(narg<2)
        error->abort("wrong command: %s",args[0]);
    if(min!=NULL)
        delete min;
    
    #define Min_Style
    #define MinStyle(class_name,style_name)     \
    else if(strcmp(args[1],#style_name)==0)     \
        {if(min!=NULL)delete min;               \
        min= new class_name(this,narg,args);    \
        min->init();min->run();min->fin();      \
        delete min;min=NULL;}
    
    if(0){}
    #include "min_styles.h"
    else
        error->abort("wrong style of minimization"
        ": %s",args[1]);
    
    #undef Min_Style
}
/*--------------------------------------------
 differnt minimization styles
 --------------------------------------------*/
void MAPP::clock_style(int narg,char** args)
{
    if(narg<2)
        error->abort("wrong command: %s",args[0]);
    if(clock!=NULL)
        delete clock;
    
    #define Clock_Style
    #define ClockStyle(class_name,style_name)   \
    else if(strcmp(args[1],#style_name)==0)     \
        {if(clock!=NULL)delete clock;           \
        clock= new class_name(this,narg,args);  \
        clock->init();clock->run();clock->fin();\
        delete clock;clock=NULL;}
    
    if(0){}
    #include "clock_styles.h"
    else
        error->abort("wrong style of clock"
                     ": %s",args[1]);
    
    #undef Clock_Style
}
/*--------------------------------------------
 differnt MD styles
 --------------------------------------------*/
void MAPP::md_style(int narg,char** args)
{
    if(narg<2)
        error->abort("wrong command: %s",args[0]);
    
    int nh_xist=0;
    TYPE0 t_step = 0.0,boltz=0.0;
    if(md!=NULL)
    {
        nh_xist=1;
        t_step=md->dt;
        boltz=md->boltz;
        delete md;
    }
    #define MD_Style
    #define MDStyle(class_name,style_name)      \
    else if(strcmp(args[1],#style_name)==0)     \
    md= new class_name(this,narg,args);

    if(0){}
    #include "md_styles.h"
    else
        error->abort("wrong style of md: %s",args[1]);
    #undef MD_Style
    if(nh_xist)
    {
        md->dt=t_step;
        md->boltz=boltz;
    }
}
/*--------------------------------------------
 differnt read styles
 --------------------------------------------*/
void MAPP::read_style(int narg,char** args)
{
    if(narg<2)
        error->abort("wrong command: %s",args[0]);
    
    Read* read;
    #define Read_Style
    #define ReadStyle(class_name,style_name)    \
    else if(strcmp(args[1],#style_name)==0)     \
    read= new class_name(this,narg,args);
    
    if(0){}
    #include "read_styles.h"
    else
        error->abort("wrong style of md: %s",args[1]);
    #undef Read_Style
    delete read;
    
}
/*--------------------------------------------
 differnt read styles
 --------------------------------------------*/
void MAPP::write_style(int narg,char** args)
{
    if(narg<2)
        error->abort("wrong command: %s",args[0]);

    if(write!=NULL)
        delete write;

    #define Write_Style
    #define WriteStyle(class_name,style_name)   \
    else if(strcmp(args[1],#style_name)==0)     \
    write= new class_name(this,narg,args);
    
    if(0){}
    #include "write_styles.h"
    else
        error->abort("wrong style of write:"
            " %s",args[1]);
    #undef Write_Style

    
}
/*--------------------------------------------
 differnt command styles
 --------------------------------------------*/
void MAPP::command_style(int narg,char** args)
{
    if(narg<2)
        error->abort("wrong command: %s",args[0]);
    
    
    #define Command_Style
    #define CommandStyle(class_name,style_name) \
    else if(strcmp(args[0],#style_name)==0){    \
    class class_name* command =                 \
    new class_name(this,narg,args);             \
    delete command;}
    
    if(0){}
    #include "command_styles.h"
    else
        error->abort("unknown command:"
                     " %s",args[0]);
    #undef Command_Style
    
    
}
/*--------------------------------------------
 differnt command styles
 --------------------------------------------*/
void MAPP::change_mode(int narg,char** args)
{    
    if(narg!=2)
        error->abort("wrong command: %s",args[0]);
    int new_mode;
    if(strcmp(args[1],"md")==0)
    {
        new_mode=MD;
    }
    else if(strcmp(args[1],"dmd")==0)
    {
        new_mode=DMD;
    }
    else
        error->abort("unknown mode: %s",args[1]);
    if(new_mode==mode)
        return;
    mode=new_mode;
    
    /*
    int x_d_n;
    int f_n;
    if(new_mode==DMD && mode==MD)
    {
        x_d_n=atoms->find_exist("x_d");
        if(x_d_n>=0)
        {
            if(atoms->my_p_no==0)
                fprintf(output,"chnaging the mode to dmd,"
                "all the velocity values will be lost\n");
            atoms->del(x_d_n);
        }
        atoms->vectors[0].change_dimension(3+atom_types->no_types);
        TYPE0* x;
        atoms->vectors[0].ret(x);
        for(int i=0;i<atoms->natms;i++)
            x[4*i+3]=0.0;
        
        f_n=atoms->find_exist("f");
        if(f_n>=0)
        {
            atoms->vectors[f_n].change_dimension(4);
            TYPE0* f;
            atoms->vectors[f_n].ret(f);
        }
        else
        {
            atoms->add<TYPE0>(0,4,"f");
        }
    }
    else if(new_mode==MD && mode==DMD)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"chnaging the mode to dmd,"
            "all the alpha values will be lost\n");
        atoms->vectors[0].change_dimension(3);
        atoms->vectors[atoms->find("f")].change_dimension(3);
    }
    
    mode=new_mode;
 
    
    
    
    
    
    if(new_mode==VG && mode==MD)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"chnaging the mode from md to vg,"
            "all the velocity values will be lost\n");
    }
    else if(new_mode==MD && mode==VG)
    {
        if(atoms->my_p_no==0)
            fprintf(output,"chnaging the mode from vg to md,"
            "all the alpha values will be lost\n");
        atoms->vectors[0].change_dimension(3);
        atoms->vectors[atoms->find("f")].change_dimension(3);
        atoms->add<TYPE0>(0,4,"f");
    }
     */
}
/*--------------------------------------------
 parse a command line:
 chops up a 1d array line into 2d array of
 char, also returns the number of arguments
 if the line starts with # it returns 0
 for example:
 line=
 "  this  is a    test, #here is comment";
 will be converted to:
 arg[0]="this";
 arg[1]="is";
 arg[2]="a";
 arg[3]="test,";
 narg=4;
 --------------------------------------------*/
int MAPP::parse_line(char* line,char**& arg)
{
    arg = NULL;
    int narg = 0;
    int cursor = 0;
    int length= static_cast<int>(strlen(line));
    int ll[length];
    
    while (cursor < length
           && line[cursor]!='#'){
        if (isspace(line[cursor]))
            cursor++;
        else
        {
            if (narg==0 && line[cursor]=='#')
                return 0;
            int i = 0;
            while(!isspace(line[cursor])
                  && cursor < length)
            {
                cursor++;
                i++;
            }
            ll[narg]=i+1;
            narg++;
        }
    }
    
    if (narg==0) return 0;
    
    CREATE1D(arg,narg);
    
    for(int i=0;i<narg;i++)
        CREATE1D(arg[i],ll[i]);
    
    narg = 0;
    cursor = 0;
    
    while (cursor<length&&line[cursor]!='#')
    {
        if ( isspace(line[cursor]))
            cursor++;
        else
        {
            int i=0;
            while(!isspace(line[cursor])
                  && cursor < strlen(line))
                arg[narg][i++]=line[cursor++];
            arg[narg][i] = '\0';
            narg++;
        }
    }
    return narg;
}
/*--------------------------------------------
 hash sign remover:
 it removes the statement after hash sign,
 replaces the withe spaces with space;
 for example:
 line=
 "  this  is a    test, #here is comment";
 will be converted to:
 newline=
 "this is a test,";
 narg=4;
 --------------------------------------------*/
int MAPP::hash_remover(char* line,char*& newline)
{
    newline = NULL;
    CREATE1D(newline,MAXCHAR);
    
    int narg = 0;
    int cursor = 0;
    int icursor = 0;
    while (cursor < strlen(line)&& line[cursor]!='#')
    {
        if (isspace(line[cursor]))
            cursor++;
        
        else
        {   if (narg!=0)
            newline[icursor++]=' ';
            while(!isspace(line[cursor])
                  && cursor < strlen(line))
                newline[icursor++]=line[cursor++];
            narg++;
        }
    }
    newline[icursor]='\0';
    return narg;
}
/*--------------------------------------------
 test
 --------------------------------------------*/
void MAPP::test()
{
}

