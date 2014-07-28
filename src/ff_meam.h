
#ifdef FF_Style
    FFStyle(ForceField_meam,MEAM)
#else
#ifndef __MAPP__ff_meam__
#define __MAPP__ff_meam__

#include <iostream>
#include "ff.h"
#include "atoms.h"
namespace MAPP_NS
{
    class ForceField_meam : public ForceField
    {
    private:
    protected:
        
        /*---------------------------------*/
        int* v2d;
        int* v3d;
        int*** vind3d;
        int** vind2d;
        /*---------------------------------*/
        
        
        
        /*---------------------------------*/
        int emb_lin_neg;//not set yet
        int ialloy;//not set yet
        int mix_ref_t;//not set yet
        int bkgd_dyn;//not set yet
        TYPE0 delr_meam;//not set yet
        TYPE0 delr_meam_inv;//not set yet
        TYPE0 rc_meam;//not set yet
        TYPE0 gsmooth_factor;//not set yet
        
        TYPE0*** c_min;//not set yet
        TYPE0*** c_max;//not set yet
        
        TYPE0** re_meam;//not set yet
        TYPE0** ebound_meam;//not set yet
        TYPE0** Ec_meam;//not set yet
        int** lattice;//not set yet
        
        TYPE0* rho0_meam;//not set yet
        TYPE0* beta0_meam;//not set yet
        TYPE0* beta1_meam;//not set yet
        TYPE0* beta2_meam;//not set yet
        TYPE0* beta3_meam;//not set yet
        TYPE0* t1_meam;//not set yet
        TYPE0* t2_meam;//not set yet
        TYPE0* t3_meam;//not set yet
        TYPE0* Z_meam;//not set yet
        TYPE0* rho_ref_meam;//not set yet
        TYPE0* A_meam;//not set yet
        int* ibar_meam;//not set yet
        /*---------------------------------*/
        
        
        
        /*---------------------------------*/
        TYPE0* scrfcn;//not set yet
        TYPE0* dscrfcn;//not set yet
        TYPE0* fcpair;//not set yet
        /*---------------------------------*/
        
        
        /*---------------------------------*/
        int rho_dim;
        TYPE0* rho;//not set yet
        TYPE0* rho_vec;//not set yet it does not need to be communicated
        TYPE0* x;
        int* type;
        /*---------------------------------*/
        
        void screen();
        void ff();

        void fcut(TYPE0,TYPE0&);
        void dfcut(TYPE0,TYPE0&,TYPE0&);
        void dCfunc(TYPE0,TYPE0,TYPE0,TYPE0&);
        void dCfunc2(TYPE0,TYPE0,TYPE0,TYPE0&,TYPE0&);
        
        void G_gam(TYPE0,int,TYPE0,TYPE0&);
        void dG_gam(TYPE0,int,TYPE0,TYPE0&,TYPE0&);
        void get_shpfcn(TYPE0*,int);
    public:
        ForceField_meam(MAPP *);
        ~ForceField_meam();
        void force_calc(int,double*);
        TYPE0 energy_calc();
        void init();
        void fin();
        void coef(int,char**){}

    };
}


#endif
#endif
