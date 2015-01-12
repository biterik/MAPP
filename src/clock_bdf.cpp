#include "clock_bdf.h"
#include "neighbor.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_BDF::Clock_BDF(MAPP* mapp,int narg
,char** arg):Clock(mapp)
{
    min_gamma=1.0e-10;
    gamma_red=0.8;
    slope=0.4;
    max_iter=50;
    max_order=6;
    m_tol=1.0e-10;
    a_tol=1.0e-6;
    e_tol=0.0;
    min_del_t=1.0e-12;
    max_del_t=1.0e4;
    initial_del_t=-1.0;
    
    if(narg<3)
        error->abort("invalid arguments");
    no_steps=atoi(arg[2]);
    if(no_steps<=0)
        error->abort("no of steps should "
        "be greater than zero");
    
    if(narg>3)
    {
        if((narg-3)%2!=0)
            error->abort("wrong number of inputs");
        int iarg=3;
        while(iarg<narg)
        {
            if(strcmp(arg[iarg],"min_gamma")==0)
            {
                iarg++;
                min_gamma=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"red_gamma")==0)
            {
                iarg++;
                gamma_red=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"slope")==0)
            {
                iarg++;
                slope=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_iter")==0)
            {
                iarg++;
                max_iter=atoi(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_order")==0)
            {
                iarg++;
                max_order=atoi(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"m_tol")==0)
            {
                iarg++;
                m_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"a_tol")==0)
            {
                iarg++;
                a_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"e_tol")==0)
            {
                iarg++;
                e_tol=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"min_del_t")==0)
            {
                iarg++;
                min_del_t=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"max_del_t")==0)
            {
                iarg++;
                max_del_t=atof(arg[iarg]);
                iarg++;
            }
            else if(strcmp(arg[iarg],"initial_del_t")==0)
            {
                iarg++;
                initial_del_t=atof(arg[iarg]);
                if(initial_del_t<=0.0)
                    error->abort("initial_del_t should be larger than zero");
                iarg++;
            }
            else
                error->abort("unknown keyword: %s",arg[iarg]);
        }
    }
    
    if(min_gamma<=0.0)
        error->abort("min_gamma should be larger than zero");
    if(gamma_red<=0.0 || gamma_red>=1.0)
        error->abort("gamma_red should be between zero and one");
    if(slope<=0.0 || slope>=1.0)
        error->abort("slope should be between zero and one");
    if(max_iter<=0)
        error->abort("max_iter should be larger than zero");
    if(max_order<=0)
        error->abort("max_order should be larger than zero");
    if(m_tol<=0.0)
        error->abort("m_tol should be larger than zero");
    if(a_tol<=0.0)
        error->abort("a_tol should be larger than zero");
    if(e_tol<0.0)
        error->abort("e_tol should be larger than zero");
    if(min_del_t<=0.0)
        error->abort("min_del_t should be larger than zero");
    if(max_del_t<=0.0)
        error->abort("max_del_t should be larger than zero");
    if(min_del_t>=max_del_t)
        error->abort("max_del_t should be larger than min_del_t");
    
    
    c_n=atoms->find("c");
    c_d_n=atoms->find("c_d");
    
    int c_dim=atoms->vectors[c_n].dim;
    dof_lcl=atoms->natms*c_dim;
    MPI_Allreduce(&dof_lcl,&dof_tot,1,MPI_INT,MPI_SUM,world);
    
    CREATE1D(y_0,dof_lcl);
    CREATE1D(e_n,dof_lcl);
    CREATE1D(a,dof_lcl);
    CREATE1D(g,dof_lcl);
    CREATE1D(g0,dof_lcl);
    CREATE1D(c0,dof_lcl);
    CREATE1D(h,dof_lcl);
    CREATE1D(dy,dof_lcl);
    
    CREATE1D(t,max_order);
    CREATE1D(y,max_order);
    for(int i=0;i<max_order;i++)
        CREATE1D(y[i],dof_lcl);
    
    CREATE1D(alpha_y,max_order);
    CREATE1D(dalpha_y,max_order);
    CREATE1D(coef,max_order);
    CREATE1D(lwr_alpha,max_order);
    
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_BDF::~Clock_BDF()
{
    
    
    delete thermo;
    
    for(int i=0;i<max_order;i++)
        if(dof_lcl)
            delete [] y[i];
    if(max_order)
    {
        delete [] y;
        delete [] t;
        delete [] alpha_y;
        delete [] dalpha_y;
        delete [] coef;
        delete [] lwr_alpha;
    }
    
    if(dof_lcl)
    {
        delete [] y_0;
        delete [] e_n;
        delete [] a;
        delete [] g;
        delete [] g0;
        delete [] c0;
        delete [] h;
        delete [] dy;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_BDF::init()
{
    TYPE0* c;
    TYPE0* c_d;
    
    int f_n=atoms->find_exist("f");
    if(f_n<0)
        f_n=atoms->add<TYPE0>(0,atoms->vectors[0].dim,"f");
    
    if(cdof_n>-1)
    {
        int dof_n=atoms->find("dof");
        vecs_comm=new VecLst(mapp,5,0,c_n,c_d_n,cdof_n,dof_n);
        
    }
    else
    {
        vecs_comm=new VecLst(mapp,3,0,c_n,c_d_n);
    }
    
    vecs_comm->add_update(0);
    atoms->reset_comm(vecs_comm);
    forcefield->init();
    atoms->ph_setup(1,vecs_comm);
    
    neighbor->init();
    neighbor->create_list(0,1);
    atoms->store_0();
    forcefield->create_2nd_neigh_lst();
    
    
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    
    forcefield->force_calc(1,energy_stress);
    
    thermo->update(fe_idx,energy_stress[0]);
    thermo->update(stress_idx,6,&energy_stress[1]);
    thermo->update(time_idx,0.0);
    
    if(write!=NULL)
        write->init();
    thermo->init();
    
    delete [] energy_stress;
    
    
    t[0]=0.0;
    thermo->start_force_time();
    forcefield->c_d_calc();
    thermo->stop_force_time();
    atoms->vectors[c_n].ret(c);
    atoms->vectors[c_d_n].ret(c_d);
    
    rectify(c_d);
    if(initial_del_t<0.0)
    {
        TYPE0 sum,sum_lcl;
        sum_lcl=0.0;
        for(int i=0;i<dof_lcl;i++)
            sum_lcl+=c_d[i]*c_d[i];
        sum=0.0;
        MPI_Allreduce(&sum_lcl,&sum,1,MPI_TYPE0,MPI_SUM,world);
        sum=sqrt(sum/static_cast<TYPE0>(dof_tot));
        initial_del_t=MAX(a_tol/sum,min_del_t);
        initial_del_t=MIN(initial_del_t,max_del_t);
    }
    
    memcpy(y[0],c,dof_lcl*sizeof(TYPE0));
    memcpy(dy,c_d,dof_lcl*sizeof(TYPE0));
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_BDF::fin()
{
    
    if(write!=NULL)
        write->fin();
    
    forcefield->fin();
    neighbor->fin();
    
    delete vecs_comm;
    thermo->fin();
    
    
    atoms->x2s(atoms->natms);
}
/*--------------------------------------------
 calculate the coefficients
 --------------------------------------------*/
int Clock_BDF::interpolate(TYPE0 del_t,int q)
{
    TYPE0 tmp0,tmp1,c0,c1,c2,c3,c4,t_new;
    t_new=t[0]+del_t;
    c0=1.0;
    c1=0.0;
    c2=0.0;
    c4=1.0;
    for(int i=1;i<q;i++)
    {
        c4*=1.0/(t[0]-t[i]);
        tmp0=(t_new-t[0])/(t[i]-t[0]);
        tmp1=tmp0*tmp0;
        for(int j=1;j<q;j++)
        {
            if(i!=j)
                tmp1*=(t_new-t[j])/(t[i]-t[j]);
        }
        alpha_y[i]=tmp1;
        
        c0*=(t_new-t[i])/(t[0]-t[i]);
        c1+=1.0/(t[0]-t[i]);
        c2+=1.0/(t_new-t[i]);
    }
    
    alpha_y[0]=(1.0-(t_new-t[0])*c1)*c0;
    alpha_dy_0=c0*(t_new-t[0]);
    
    c3=2.0/(t_new-t[0]);
    
    
    dalpha_dy_0=alpha_dy_0*(c2+0.5*c3);
    dalpha_y[0]=c0*(c2-c1+(t_new-t[0])*c1*c2);
    
    for(int i=1;i<q;i++)
        dalpha_y[i]=alpha_y[i]*(c2+c3-1.0/(t_new-t[i]));
    
    
    
    tmp0=0.0;
    for(int i=0;i<q;i++)
        tmp0+=1.0/static_cast<TYPE0>(i+1);
    beta=del_t/tmp0;
    
    
    int ret_val=1;
    int tot_ret_val;
    int idof=0;
    while(idof<dof_lcl && ret_val==1)
    {
        
        a[idof]=beta*dalpha_dy_0*dy[idof];
        y_0[idof]=alpha_dy_0*dy[idof];
        for(int j=0;j<q;j++)
        {
            a[idof]+=beta*dalpha_y[j]*y[j][idof];
            y_0[idof]+=alpha_y[j]*y[j][idof];
        }
        
        a[idof]-=y_0[idof];
        if(y_0[idof]<0.0 || y_0[idof]>1.0)
        {
            ret_val=-1;
        }
        
        
        idof++;
    }
    
    MPI_Allreduce(&ret_val,&tot_ret_val,1,MPI_INT,MPI_MIN,world);
    
    return tot_ret_val;
}

/*--------------------------------------------
 lower order error estimate
 --------------------------------------------*/
void Clock_BDF::hi_lo_err_est(TYPE0 del_t,int q)
{
    
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    
    TYPE0 curr_t=t[0]+del_t;

    
    TYPE0 c0=0.0,c1=1.0,tmp0;
    TYPE0 lo_err_lcl,hi_err_lcl;
    
    
    TYPE0 a_0,a_1,l1_q;
    TYPE0 hi_a_0,hi_a_1,hi_l1_q,Q_n=1.0,hi_err_pre_fac=1.0;
    TYPE0 lo_l1_q,lo_err_pre_fac;
    
    a_0=a_1=1.0;
    l1_q=0.0;
    
    for(int i=0;i<q;i++)
    {
        l1_q+=1.0/(curr_t-t[i]);
        a_0*=curr_t-t[i];
        if(i!=0)
            a_1*=(curr_t-t[i])/(t[0]-t[i]);
    }
    
    if(hi_ord_avail)
    {
        hi_a_0=hi_a_1=1.0;
        for(int i=0;i<q;i++)
        {
            hi_a_0*=t[0]-t[i+1];
            if(i!=0)
            {
                hi_a_1*=(t[0]-t[i+1])/(t[1]-t[i+1]);
            }
        }
        Q_n=(a_0/hi_a_0)*(1.0+a_1)/(1.0+hi_a_1)*(del_t/(t[0]-t[1]));
        hi_l1_q=l1_q+1.0/(curr_t-t[q]);
        hi_err_pre_fac=-(curr_t-t[q])/(q+2);
        hi_err_pre_fac*=1.0/(hi_l1_q*(1.0+a_1)*del_t*del_t);
    }
    
    if(lo_ord_avail)
    {
        lo_l1_q=l1_q-1.0/(curr_t-t[q-1]);
        lo_err_pre_fac=-del_t*a_0/(curr_t-t[q-1]);
        c0=0.0;
        c1=1.0;
        for(int i=0;i<q-1;i++)
        {
            tmp0=1.0/(curr_t-t[i]);
            c0+=tmp0;
            c1*=tmp0;
            lwr_alpha[i]=tmp0*tmp0;
            
            for(int j=0;j<q-1;j++)
                if(i!=j)
                    lwr_alpha[i]*=1.0/(t[i]-t[j]);
            
            lwr_alpha[i]*=lo_err_pre_fac;
        }
        lwr_alpha_dy=c1*lo_err_pre_fac;
        lwr_alpha_y=-c0*c1*lo_err_pre_fac;
    }
    
    lo_err_lcl=hi_err_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        
        if(hi_ord_avail)
        {
            tmp0=(c[i]-y_0[i])-Q_n*e_n[i];
            tmp0*=hi_err_pre_fac;
            hi_err_lcl+=tmp0*tmp0;
        }
        
        if(lo_ord_avail)
        {
            tmp0=lwr_alpha_y*c[i];
            tmp0+=lwr_alpha_dy*c_d[i];
            for(int j=0;j<q-1;j++)
                tmp0+=lwr_alpha[j]*y[j][i];
            
            tmp0*=c1/c0;
            lo_err_lcl+=tmp0*tmp0;
        }
        
        e_n[i]=c[i]-y_0[i];
        
    }
    
    lo_err=hi_err=0.0;
    MPI_Allreduce(&lo_err_lcl,&lo_err,1,MPI_TYPE0,MPI_SUM,world);
    MPI_Allreduce(&hi_err_lcl,&hi_err,1,MPI_TYPE0,MPI_SUM,world);
    
    lo_err=sqrt(lo_err/(static_cast<TYPE0>(dof_tot)*a_tol));
    hi_err=sqrt(hi_err/(static_cast<TYPE0>(dof_tot)*a_tol));
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_BDF::run()
{
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    TYPE0* tmp_y;
    
    TYPE0 del_t=initial_del_t,del_t_tmp,err1;
    TYPE0 ratio=1.0,cost,lo_ratio=1.0,hi_ratio=1.0;
    int ord=1;
    int chk;
    int const_stp=0;
    int initial_phase=0;
    
    eq_ratio=1.0;
    
    int istep=0;
    while (eq_ratio>=1.0 && istep <no_steps)
    {
        chk=interpolate(del_t,ord);
        while(chk==-1)
        {
            if(initial_phase) initial_phase=0;
            del_t*=0.8;
            del_t=MAX(del_t,min_del_t);
            chk=interpolate(del_t,ord);
        }
        cost=solve(del_t,ord);
        //cost=0.0;
        while(err>=1.0 || cost>=1.0)
        {
            if(initial_phase) initial_phase=0;
            
            err1=MAX(err,cost);
            ratio=pow(0.5/err1,1.0/static_cast<TYPE0>(ord+1));
            
            if(ratio<0.5)
                ratio=0.5;
            else if(ratio>0.9)
                ratio=0.9;
            
            del_t_tmp=del_t*ratio;
            if(del_t_tmp<min_del_t)
                del_t=min_del_t;
            else
                del_t=del_t_tmp;
            
            
            chk=interpolate(del_t,ord);
            
            while(chk==-1)
            {
                del_t*=0.8;
                chk=interpolate(del_t,ord);
            }
            
            cost=solve(del_t,ord);
            //cost=0.0;
        }
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step())
        {
            thermo->start_force_time();
            forcefield->force_calc(1,energy_stress);
            thermo->stop_force_time();
            thermo->update(fe_idx,energy_stress[0]);
            thermo->update(stress_idx,6,&energy_stress[1]);
            thermo->update(time_idx,t[0]+del_t);
        }
        
        
        
        tmp_y=y[max_order-1];
        for(int i=max_order-1;i>0;i--)
        {
            t[i]=t[i-1];
            y[i]=y[i-1];
        }
        
        y[0]=tmp_y;
        t[0]+=del_t;
        
        memcpy(y[0],c,dof_lcl*sizeof(TYPE0));
        memcpy(dy,c_d,dof_lcl*sizeof(TYPE0));
        
        
        
        
        
        if(initial_phase)
        {
            del_t=MIN(del_t*2.0,max_del_t);
            if(ord<max_order && istep>=ord+3)
                ord++;
        }
        else
        {
            if(ord!=1 && const_stp>ord+1)
                lo_ord_avail=1;
            else
                lo_ord_avail=0;
            
            if((ord<max_order && istep!=0) && const_stp>ord+1)
                hi_ord_avail=1;
            else
                hi_ord_avail=0;
            
            hi_lo_err_est(del_t,ord);
            
            if(hi_ord_avail && lo_ord_avail)
            {
                lo_ratio=pow(0.5/lo_err,1.0/static_cast<TYPE0>(ord));
                hi_ratio=pow(0.5/hi_err,1.0/static_cast<TYPE0>(ord+2));
                ratio=pow(0.5/err,1.0/static_cast<TYPE0>(ord+1));
                
                if(hi_ratio>lo_ratio && hi_ratio>ratio)
                {
                    ord++;
                    ratio=hi_ratio;
                }
                else if(lo_ratio>hi_ratio && lo_ratio>ratio)
                {
                    ord--;
                    ratio=lo_ratio;
                }
                
            }
            else if(hi_ord_avail==1&& lo_ord_avail==0)
            {
                hi_ratio=pow(0.5/hi_err,1.0/static_cast<TYPE0>(ord+2));
                ratio=pow(0.5/err,1.0/static_cast<TYPE0>(ord+1));
                
                if(hi_ratio>ratio)
                {
                    ord++;
                    ratio=hi_ratio;
                }
            }
            else if(hi_ord_avail==0 && lo_ord_avail==1)
            {
                lo_ratio=pow(0.5/lo_err,1.0/static_cast<TYPE0>(ord));
                ratio=pow(0.5/err,1.0/static_cast<TYPE0>(ord+1));
                if(lo_ratio>ratio)
                {
                    ord--;
                    ratio=lo_ratio;
                }
            }
            else
                ratio=pow(0.5/err,1.0/static_cast<TYPE0>(ord+1));
            
            if(ratio>=2.0)
            {
                ratio=2.0;
                const_stp=0;
            }
            else if(ratio<=0.5)
            {
                ratio=0.5;
                const_stp=0;
            }
            else
            {
                ratio=1.0;
                const_stp++;
            }
            
        }
        
        
        del_t_tmp=ratio*del_t;
        if(del_t_tmp>max_del_t)
            del_t=max_del_t;
        else if (del_t_tmp<min_del_t)
            del_t=min_del_t;
        else
            del_t=del_t_tmp;
        

        step_no++;
        istep++;
    }
    

    delete [] energy_stress;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
TYPE0 Clock_BDF::solve(TYPE0 del_t,int q)
{
    
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    TYPE0 gamma,max_gamma=1.0;
    TYPE0 inner,tmp,tmp1;
    TYPE0 ratio;
    TYPE0 g0_g0,g_g,g_g0,g_h;
    TYPE0 curr_cost,ideal_cost,cost;
    int chk;
    

    /*
    memcpy(c,y_0,dof_lcl*sizeof(TYPE0));
    thermo->start_comm_time();
    atoms->update(c_n);
    thermo->stop_comm_time();
    */
   
    
    TYPE0 tot_ratio;
    ratio=1.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp=dy[i]*del_t+y[0][i];
        if(tmp>1.0)
            ratio=MIN((1.0-y[0][i])/(tmp-y[0][i]),ratio);
        else if(tmp<0.0)
            ratio=MIN((y[0][i]-0.0)/(y[0][i]-tmp),ratio);
    }
    MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
    for(int i=0;i<dof_lcl;i++)
        c[i]=y[0][i]+dy[i]*del_t*tot_ratio;

    thermo->start_comm_time();
    atoms->update(c_n);
    thermo->stop_comm_time();
    
    
    
    thermo->start_force_time();
    curr_cost=forcefield->calc_g(0,beta,a,g);
    rectify(g);
    thermo->stop_force_time();
    memcpy(h,g,dof_lcl*sizeof(TYPE0));
    
    inner=0.0;
    for(int i=0;i<dof_lcl;i++)
        inner+=g[i]*g[i];
    g0_g0=0.0;
    MPI_Allreduce(&inner,&g0_g0,1,MPI_TYPE0,MPI_SUM,world);
    g_h=g0_g0;
    
    int iter=0;
    //printf("1st cost: %e \n",curr_cost);
    
    while(curr_cost>m_tol*static_cast<TYPE0>(dof_tot)
          && iter<max_iter && max_gamma>min_gamma)
    {
        memcpy(g0,g,dof_lcl*sizeof(TYPE0));
        memcpy(c0,c,dof_lcl*sizeof(TYPE0));
        
        gamma=0.99;
        
        for(int i=0;i<dof_lcl;i++)
        {
            tmp=c0[i]+h[i];
            
            if(tmp>1.0)
            {
                gamma=MIN((1.0-c0[i])/(tmp-c0[i]),gamma);
            }
            else if(tmp<0.0)
            {
                gamma=MIN((c0[i]-0.0)/(c0[i]-tmp),gamma);
            }
        }
        MPI_Allreduce(&gamma,&max_gamma,1,MPI_TYPE0,MPI_MIN,world);

        max_gamma*=0.9;
        chk=1;
        
        cost=curr_cost;
        while(chk && max_gamma>min_gamma)
        {
            for(int i=0;i<dof_lcl;i++)
                c[i]=c0[i]+max_gamma*h[i];
            
            thermo->start_comm_time();
            atoms->update(c_n);
            thermo->stop_comm_time();
            
            thermo->start_force_time();
            curr_cost=forcefield->calc_g(0,beta,a,g);
            rectify(g);
            thermo->stop_force_time();

            ideal_cost=cost-slope*max_gamma*g_h;
            if(curr_cost<ideal_cost)
                chk=0;
            max_gamma*=gamma_red;
            if(max_gamma<=min_gamma)
            {
                memcpy(c,c0,dof_lcl*sizeof(TYPE0));
                
                thermo->start_comm_time();
                atoms->update(c_n);
                thermo->stop_comm_time();
                
                thermo->start_force_time();
                curr_cost=forcefield->calc_g(1,beta,a,g);
                thermo->stop_force_time();
            }
        }
        
        if(chk==0)
        {
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                inner+=g[i]*g0[i];
            g_g0=0.0;
            MPI_Allreduce(&inner,&g_g0,1,MPI_TYPE0,MPI_SUM,world);
            
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
                inner+=g[i]*g[i];
            g_g=0.0;
            MPI_Allreduce(&inner,&g_g,1,MPI_TYPE0,MPI_SUM,world);
            
            ratio=(g_g-g_g0)/g0_g0;
            
            g0_g0=g_g;
            
            inner=0.0;
            for(int i=0;i<dof_lcl;i++)
            {
                h[i]*=ratio;
                h[i]+=g[i];
                inner+=h[i]*g[i];
            }
            g_h=0.0;
            MPI_Allreduce(&inner,&g_h,1,MPI_TYPE0,MPI_SUM,world);
            
            if(g_h<0.0)
            {
                memcpy(h,g,dof_lcl*sizeof(TYPE0));
                g_h=g_g;
            }
        }
        
        
        iter++;
    }
    
    rectify(c_d);
    //printf("2nd cost: %e \n\n",curr_cost);
    
    TYPE0 a_0,a_1,curr_t,l1_q,err_pre_fac,tmp0,err_lcl;
    a_0=a_1=1.0;
    l1_q=0.0;
    
    for(int i=0;i<q;i++)
    {
        l1_q+=1.0/(curr_t-t[i]);
        a_0*=curr_t-t[i];
        if(i!=0)
        {
            a_1*=(curr_t-t[i])/(t[0]-t[i]);
        }
        
        
    }
    
    err_pre_fac=-1.0/(l1_q*(1.0+a_1));
    
    err_lcl=0.0;
    err=0.0;
    tmp1=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp0=(y_0[i]-c[i]);
        tmp1+=c_d[i]*c_d[i];
        err_lcl+=tmp0*tmp0;
    }
    
    MPI_Allreduce(&err_lcl,&err,1,MPI_TYPE0,MPI_SUM,world);
    MPI_Allreduce(&tmp1,&eq_ratio,1,MPI_TYPE0,MPI_SUM,world);
    eq_ratio=sqrt(eq_ratio/static_cast<TYPE0>(dof_tot))/e_tol;
    
    err=sqrt(err/(static_cast<TYPE0>(dof_tot)*a_tol));
    err*=fabs(err_pre_fac);

    
    
    for(int i=0;i<dof_lcl;i++)
        if(c[i]<0.0 || c[i]>1.0)
            error->abort("c exceeded the domain");
    
    
    return curr_cost/(m_tol*static_cast<TYPE0>(dof_tot));
    
}
/*--------------------------------------------
 for error
 --------------------------------------------*/
void Clock_BDF::err_coef(TYPE0* coef,int q)
{
    
    TYPE0 c0=1.0,tmp0;
    coef[0]=0.0;
    for(int i=1;i<q;i++)
    {
        coef[i]=0.0;
        coef[0]+=1.0/(t[0]-t[i]);
        c0*=t[0]-t[i];
    }
    
    for(int i=1;i<q;i++)
    {
        tmp0=1.0;
        for(int j=0;j<q;j++)
        {
            if(i!=j)
                tmp0*=t[i]-t[j];
        }
        coef[i]=c0/(tmp0*(t[0]-t[i]));
    }
    
}
