
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
        
        TYPE0* delij;
        TYPE0* delik;
        TYPE0* deljk;
        TYPE0* s;
        TYPE0* si;
        TYPE0* sj;
        TYPE0* drho1drm1;
        TYPE0* drho1drm2;
        TYPE0* drho2drm1;
        TYPE0* drho2drm2;
        TYPE0* drho3drm1;
        TYPE0* drho3drm2;
        TYPE0* drhodrm1;
        TYPE0* drhodrm2;
        TYPE0* dUdrijm;
        TYPE0* fi;
        TYPE0* fj;
        TYPE0* v;
        
    protected:
        /*---------------------------------*/
        int rho_dim;
        int rho_n,rho_vec_n,type_n,x_n,f_n;
        /*---------------------------------*/
        
        
        
        /*---------------------------------*/
        int* v2d;
        int* v3d;
        int*** vind3d;
        int** vind2d;
        /*---------------------------------*/
        
        /*
         rho dimension is 38!!!!
         **** rho functions (4)
         00: rho0
         01: rho1
         02: rho2
         03: rho3
         
         **** rho_1 terms (3)
         04: Arho1_1
         05: Arho1_2
         06: Arho1_3
         
         **** rho_2 terms (7)
         ******** 1st terms (6)
         07: Arho2_1
         08: Arho2_2
         09: Arho2_3
         10: Arho2_4
         11: Arho2_5
         12: Arho2_6
         ******** 2nd term (1)
         13: Arho2b
         
         **** rho_3 terms
         ******** 1st terms (10)
         14: Arho3_1
         15: Arho3_2
         16: Arho3_3
         17: Arho3_4
         18: Arho3_5
         19: Arho3_6
         20: Arho3_7
         21: Arho3_8
         22: Arho3_9
         23: Arho3_10
         ******** 2nd terms (3)
         24: Arho3b_1
         25: Arho3b_2
         26: Arho3b_3
         
         **** t_ave terms (3)
         27: t_ave_1
         28: t_ave_2
         29: t_ave_3
         
         **** tsq_ave terms (3)
         30: tsq_ave_1
         31: tsq_ave_2
         32: tsq_ave_3
         
         **** Gamma function (1)
         33: gamma
         
         **** dGamma functions (3)
         34: dgamma1
         35: dgamma2
         36: dgamma3
         
         **** fhop term  (1)
         37: fhop
         */
        
        /*---------------------------------*/
        int emb_lin_neg;//not set yet
        int ialloy;//not set yet
        int mix_ref_t;//not set yet
        int bkgd_dyn;//not set yet
        int augt1;//not set yet
        int erose_form;//not set yet
        int nr;
        TYPE0 delr_meam;//not set yet
        TYPE0 delr_meam_inv;//not set yet
        TYPE0 rc_meam;//not set yet
        TYPE0 gsmooth_factor;//not set yet
        TYPE0 dr;
        TYPE0 dr_inv;
        
        TYPE0*** c_min;//not set yet
        TYPE0*** c_max;//not set yet
        TYPE0*** phirar;
        
        TYPE0** re_meam;//not set yet
        TYPE0** ebound_meam;//not set yet
        TYPE0** Ec_meam;//not set yet
        TYPE0** alpha_meam;//not set yet
        TYPE0** delta_meam;//not set yet
        TYPE0** attrac_meam;//not set yet
        TYPE0** repuls_meam;//not set yet
        int** lattice;//not set yet
        int** nn2_meam;//not set yet
        int** zbl_meam;//not set yet
        
        TYPE0* rho0_meam;//not set yet
        TYPE0* beta0_meam;//not set yet
        TYPE0* beta1_meam;//not set yet
        TYPE0* beta2_meam;//not set yet
        TYPE0* beta3_meam;//not set yet
        TYPE0* t0_meam;//not set yet
        TYPE0* t1_meam;//not set yet
        TYPE0* t2_meam;//not set yet
        TYPE0* t3_meam;//not set yet
        TYPE0* Z_meam;//not set yet
        TYPE0* rho_ref_meam;//not set yet
        TYPE0* A_meam;//not set yet
        int* ibar_meam;//not set yet
        int* ielt_meam;//not set yet
        int* type_ref;//not set yet
        
        void get_dens_ref(TYPE0,int,int,TYPE0&,TYPE0&,TYPE0&,TYPE0&,TYPE0&,TYPE0&,TYPE0&,TYPE0&);
        void get_sijk(TYPE0,int,int,int,TYPE0&);
        void get_tavref(TYPE0&,TYPE0&,TYPE0&,TYPE0&,TYPE0&,TYPE0&,TYPE0,TYPE0,TYPE0,TYPE0,TYPE0,TYPE0,TYPE0,int,int,int);
        TYPE0 phi_meam(TYPE0,int,int);
        /*---------------------------------*/
        
        
        /*---------------------------------*/
        int max_pairs;
        TYPE0* scrfcn;//not set yet
        TYPE0* dscrfcn;//not set yet
        TYPE0* fcpair;//not set yet
        /*---------------------------------*/
                

        void reset();

        void fcut(TYPE0,TYPE0&);
        void dfcut(TYPE0,TYPE0&,TYPE0&);
        void dCfunc(TYPE0,TYPE0,TYPE0,TYPE0&);
        void dCfunc2(TYPE0,TYPE0,TYPE0,TYPE0&,TYPE0&);
        
        void G_gam(TYPE0,int,TYPE0,TYPE0&);
        void dG_gam(TYPE0,int,TYPE0,TYPE0&,TYPE0&);
        void get_shpfcn(TYPE0*,int);
        
        void read_global(char*);
        void read_local(char*);
        
        void get_Zij(TYPE0&,int);
        void get_Zij2(TYPE0&,TYPE0&,TYPE0&,int,TYPE0,TYPE0);
        TYPE0 zbl(TYPE0,TYPE0,TYPE0);
        TYPE0 erose(TYPE0,TYPE0,TYPE0,TYPE0,TYPE0,TYPE0,int);
        void alloy_params();
        void compute_reference_density();
        void compute_pair_meam();
        void compute_phi(TYPE0,int,int,TYPE0&);
        void compute_phi_dphi(TYPE0,int,int,TYPE0&,TYPE0&);
        
        TYPE0* nrgy_strss;
        TYPE0 third,sixth;
        
    public:
        ForceField_meam(MAPP *);
        ~ForceField_meam();
        void force_calc1(int,double*);
        void force_calc(int,double*);        
        TYPE0 energy_calc();
        void init();
        void fin();
        void coef(int,char**);

    };
}


#endif
#endif
