/*--------------------------------------------
 Created by Sina on 07/30/13.
 Copyright (c) 2013 MIT. All rights reserved.
 --------------------------------------------*/
#include "rand_engine.h"
#include "error.h"
#include <cmath>
#define RNG_M 2147483647
#define RNG_MI 1.0/2147483647
#define RNG_A 16807
#define RNG_Q 127773
#define RNG_R 2836
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Random::Random(MAPP* mapp,int num)
: InitPtrs(mapp)
{
    if(num<=0)
        error->abort("random seed cannot be"
        " equal or less than zero");
    seed=num;
    gauss_chk=0;
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Random::~Random()
{
    
}
/*--------------------------------------------
 uniform random number generator
 uses Park Miller algorith to generate a
 random number (0,1)
 --------------------------------------------*/
double Random::uniform()
{
    int first,second,rnd;
    
    first = seed/RNG_Q;
    second = seed%RNG_Q;
    rnd = RNG_A*second-RNG_R*first;
    if (rnd > 0)
        seed=rnd;
    else
        seed=rnd+RNG_M;
    return seed*RNG_MI;
}
/*--------------------------------------------
 gaussian random number generator
 --------------------------------------------*/
double Random::gaussian()
{
    if(gauss_chk==0)
    {
        double u=0.0,v=0.0,s,r;
        s=1.0;
        while(s>=1.0)
        {
            u=2.0*uniform()-1.0;
            v=2.0*uniform()-1.0;
            s=u*u+v*v;
        }
        
        r=sqrt(-2.0*log(s)/s);
        gauss_chk=1;
        reserved=v*r;
        return u*r;
    }
    else
    {
        gauss_chk=0;
        return reserved;
    }  
}
/*--------------------------------------------
 gaussian random number generator
 --------------------------------------------*/
double Random::
gaussian(double sigma,double mu)
{
    return (mu+sigma*gaussian());
}


