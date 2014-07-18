/*--------------------------------------------
 Created by Sina on 05/12/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/

#ifndef MAPP_init_h
#define MAPP_init_h
#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include "MAPP.h"

namespace MAPP_NS{
    class InitPtrs{
    private:
    public:
        InitPtrs(MAPP* ptr):
        mapp(ptr),
        memory(ptr->memory),
        error(ptr->error),
        forcefield(ptr->forcefield),
        world(ptr->world),
        atoms(ptr->atoms),
        neighbor(ptr->neighbor),
        atom_types(ptr->atom_types),
        output(ptr->output),
        write(ptr->write),
        step_no(ptr->step_no),
        step_tally(ptr->step_tally)
        {}
        virtual ~InitPtrs(){}
    protected:
        MAPP* mapp;
        Memory*& memory;
        Error*& error;
        ForceField*& forcefield;
        MPI_Comm &world;
        Neighbor*& neighbor;
        Atoms*& atoms;
        AtomTypes*& atom_types;
        FILE*& output;
        Write*& write;
        int& step_no;
        int& step_tally;
    };
    
}
#endif

/* some useful macros */
#define AMU2ePA 0.000103642691516828
/* some usefule functions */
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MEM_ERR(vec) printf("%s:%s:vec\n",__FILE__,__FUNCTION__)

/* some importante constants */
#define MAXCHAR 1024
#define MAXSIZE 512
#define TOLERANCE 1.0e-10

/* some vector & matrix functions */
#define V3ZERO(V) (V[0]=0, V[1]=0, V[2]=0)
#define M3DET(A) A[0][0]*A[1][1]*A[2][2]\
+A[1][0]*A[2][1]*A[0][2]\
+A[2][0]*A[0][1]*A[1][2]\
-A[0][2]*A[1][1]*A[2][0]\
-A[1][2]*A[2][1]*A[0][0]\
-A[2][2]*A[0][1]*A[1][0]
#define M3ZERO(A) ( A[0][0]=0, A[0][1]=0, A[0][2]=0, A[1][0]=0, \
A[1][1]=0, A[1][2]=0, A[2][0]=0, A[2][1]=0, A[2][2]=0 )
#define M3INV(A,B,determinant) ( \
B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1], \
B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2], \
B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0], \
B[1][0] = A[1][2]*A[2][0]-A[1][0]*A[2][2], \
B[2][1] = A[2][0]*A[0][1]-A[2][1]*A[0][0], \
B[0][2] = A[0][1]*A[1][2]-A[0][2]*A[1][1], \
B[2][0] = A[1][0]*A[2][1]-A[2][0]*A[1][1], \
B[0][1] = A[2][1]*A[0][2]-A[0][1]*A[2][2], \
B[1][2] = A[0][2]*A[1][0]-A[1][2]*A[0][0], \
(determinant) = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0], \
B[0][0] /= (determinant),B[1][1] /= (determinant),B[2][2] /= (determinant), \
B[1][0] /= (determinant),B[2][1] /= (determinant),B[0][2] /= (determinant), \
B[2][0] /= (determinant),B[0][1] /= (determinant),B[1][2] /= (determinant) )
#define M3COFAC(A,B) ( \
B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1], \
B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2], \
B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0], \
B[1][0] = A[1][2]*A[2][0]-A[1][0]*A[2][2], \
B[2][1] = A[2][0]*A[0][1]-A[2][1]*A[0][0], \
B[0][2] = A[0][1]*A[1][2]-A[0][2]*A[1][1], \
B[2][0] = A[1][0]*A[2][1]-A[2][0]*A[1][1], \
B[0][1] = A[2][1]*A[0][2]-A[0][1]*A[2][2], \
B[1][2] = A[0][2]*A[1][0]-A[1][2]*A[0][0])
#define THICKNESS(B,d) (\
d[0]=1.0/sqrt(B[0][0]*B[0][0]+B[1][0]*B[1][0]+B[2][0]*B[2][0]),\
d[1]=1.0/sqrt(B[0][1]*B[0][1]+B[1][1]*B[1][1]+B[2][1]*B[2][1]),\
d[2]=1.0/sqrt(B[0][2]*B[0][2]+B[1][2]*B[1][2]+B[2][2]*B[2][2]))
#define M3EQV(A,B) ( B[0][0] = A[0][0], B[0][1] = A[0][1], \
B[0][2] = A[0][2], B[1][0] = A[1][0], B[1][1] = A[1][1], \
B[1][2] = A[1][2], B[2][0] = A[2][0], B[2][1] = A[2][1], \
B[2][2] = A[2][2] )
#define M3IDENTITY(A) ( A[0][0]=1, A[0][1]=0, A[0][2]=0, \
A[1][0]=0, A[1][1]=1, A[1][2]=0, A[2][0]=0, A[2][1]=0, A[2][2]=1 )
#define M3MUL_TRI_LOWER(A,B,C) (C[0][0]=A[0][0]*B[0][0],C[0][1]=0,C[0][2]=0,\
C[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0],C[1][1]=A[1][1]*B[1][1],C[1][2]=0,\
C[2][0]=A[2][0]*B[0][0]+A[2][1]*B[1][0]+A[2][2]*B[2][0],\
C[2][1]=A[2][1]*B[1][1]+A[2][2]*B[2][1],C[2][2]=A[2][2]*B[2][2])
#define M3INV_TRI_LOWER(A,B) (B[0][0]=1.0/A[0][0],B[0][1]=0,B[0][2]=0,\
B[1][0]=-A[1][0]/(A[0][0]*A[1][1]),B[1][1]=1.0/A[1][1],B[1][2]=0,\
B[2][0]=(A[1][0]*A[2][1]-A[1][1]*A[2][0])/(A[0][0]*A[1][1]*A[2][2]),\
B[2][1]=-A[2][1]/(A[1][1]*A[2][2]),B[2][2]=1.0/A[2][2])
#define COMP(ityp,jtyp) ((ityp) < (jtyp) ? ((jtyp+1)*jtyp/2+ityp) : ((ityp+1)*ityp/2+jtyp))









