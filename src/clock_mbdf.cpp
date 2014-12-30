#include "clock_mbdf.h"
#include "neighbor.h"
#include "ff.h"
#include "write.h"
#include "thermo_dynamics.h"
using namespace MAPP_NS;
/*--------------------------------------------
 constructor
 --------------------------------------------*/
Clock_MBDF::Clock_MBDF(MAPP* mapp,int narg
                   ,char** arg):Clock(mapp)
{
    
    min_gamma=1.0e-30;
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
            else if(strcmp(arg[iarg],"gamma_red")==0)
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
    
    CREATE1D(dy,dof_lcl);
    CREATE1D(y_0,dof_lcl);
    CREATE1D(a,dof_lcl);
    CREATE1D(g,dof_lcl);
    CREATE1D(g0,dof_lcl);
    CREATE1D(c0,dof_lcl);
    CREATE1D(h,dof_lcl);
    
    CREATE1D(t,max_order+2);
    CREATE1D(mod_alpha,max_order+2);
    CREATE1D(mod_d_alpha,max_order+2);
    CREATE1D(y,max_order+2);
    for(int i=0;i<max_order+2;i++)
        CREATE1D(y[i],dof_lcl);
    
    CREATE1D(alph_err,max_order+2);
}
/*--------------------------------------------
 destructor
 --------------------------------------------*/
Clock_MBDF::~Clock_MBDF()
{
    
    delete thermo;
    
    for(int i=0;i<max_order+2;i++)
        if(dof_lcl) delete [] y[i];
    
    delete [] alph_err;
    delete [] y;
    delete [] t;
    delete [] mod_alpha;
    delete [] mod_d_alpha;
    
    if(dof_lcl)
    {
        delete [] dy;
        delete [] y_0;
        delete [] a;
        delete [] g0;
        delete [] c0;
        delete [] g;
        delete [] h;
    }
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_MBDF::init()
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
    
    t[0]=0.0;
    t[1]=-initial_del_t;
    memcpy(y[0],c,dof_lcl*sizeof(TYPE0));
    memcpy(y[1],c,dof_lcl*sizeof(TYPE0));
    memcpy(dy,c_d,dof_lcl*sizeof(TYPE0));
    
    
}
/*--------------------------------------------
 init
 --------------------------------------------*/
void Clock_MBDF::fin()
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
 init
 --------------------------------------------*/
int Clock_MBDF::interpolate(TYPE0 del_t,int ord)
{
    
    TYPE0 tmp0,tmp1,tmp2;
    TYPE0 curr_t=t[0]+del_t;
    TYPE0 tmp4=0.0;
    TYPE0 c0=0.0;
    for(int i=0;i<ord+1;i++)
    {
        c0+=1.0/(curr_t-t[i]);
        
        tmp0=tmp1=1.0;
        tmp2=0.0;
        for(int j=0;j<ord+1;j++)
            if(i!=j)
                tmp0*=(curr_t-t[j])/(t[i]-t[j]);
        
        
        mod_alpha[i]=tmp0;
        mod_d_alpha[i]=tmp0;
        
        if(i!=0)
            tmp4+=1.0/(static_cast<TYPE0>(i));
    }
    
    beta=del_t/tmp4;
    
    
    for(int i=0;i<ord+1;i++)
    {
        tmp0=c0-1.0/(curr_t-t[i]);
        mod_d_alpha[i]*=tmp0;
    }
    
    int idof=0;
    int ret_val=1;
    int tot_ret_val;
    while(idof<dof_lcl && ret_val==1)
    {
        y_0[idof]=0.0;
        a[idof]=0.0;
        for(int j=0;j<ord+1;j++)
        {
            a[idof]+=beta*mod_d_alpha[j]*y[j][idof];
            y_0[idof]+=mod_alpha[j]*y[j][idof];
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
 run
 --------------------------------------------*/
void Clock_MBDF::run()
{
    int dim=atoms->dimension;
    TYPE0* energy_stress;
    CREATE1D(energy_stress,dim*(dim+1)/2+1);
    
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0* c_d;
    atoms->vectors[c_d_n].ret(c_d);
    
    TYPE0* tmp_y;
    
    TYPE0 del_t=initial_del_t,del_t_tmp;
    TYPE0 err1,cost;
    TYPE0 ratio;
    int ord=1;
    int chk;
    int initial_phase=1;
    int const_stps=0,ord_const_stps=0;
    int del_ord;
    
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
                del_t=MAX(del_t,min_del_t);
                chk=interpolate(del_t,ord);
            }
            
            cost=solve(del_t,ord);
            
        }
        
        
        tmp_y=y[max_order+1];
        for(int i=max_order+1;i>0;i--)
        {
            t[i]=t[i-1];
            y[i]=y[i-1];
        }
        
        y[0]=tmp_y;
        t[0]+=del_t;
        
        memcpy(y[0],c,dof_lcl*sizeof(TYPE0));
        memcpy(dy,c_d,dof_lcl*sizeof(TYPE0));
        
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
            thermo->update(time_idx,t[0]);
        }
        
        if(initial_phase)
        {
            del_t=MIN(del_t*2.0,max_del_t);
            if(ord<max_order)
                ord++;
        }
        else
        {
            del_ord=err_calc(ord,const_stps,del_t,err);
            
            ord+=del_ord;
            
            if(del_ord==0)
                ord_const_stps++;
            else
                ord_const_stps=0;
            
            ratio=step_size(del_t,ord);
            
            if(del_ord==0)
            {
                if(ratio==1.0)
                    const_stps++;
                else
                    const_stps=0;
            }
            else
                const_stps=0;
            
            
            del_t*=ratio;
        }
        
        step_no++;

        istep++;
    }
    
    /*
    for(int istep=0;istep<no_steps;istep++)
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
                del_t=MAX(del_t,min_del_t);
                chk=interpolate(del_t,ord);
            }
            
            cost=solve(del_t,ord);
            
        }
        
        
        tmp_y=y[max_order+1];
        for(int i=max_order+1;i>0;i--)
        {
            t[i]=t[i-1];
            y[i]=y[i-1];
        }
        
        y[0]=tmp_y;
        t[0]+=del_t;
        
        memcpy(y[0],c,dof_lcl*sizeof(TYPE0));
        memcpy(dy,c_d,dof_lcl*sizeof(TYPE0));
        
        if(write!=NULL)
            write->write();
        thermo->thermo_print();
        
        if(thermo->test_prev_step())
        {
            forcefield->force_calc(1,energy_stress);
            thermo->update(fe_idx,energy_stress[0]);
            thermo->update(stress_idx,6,&energy_stress[1]);
            thermo->update(time_idx,t[0]);
        }
        
        if(initial_phase)
        {
            del_t=MIN(del_t*2.0,max_del_t);
            if(ord<max_order)
                ord++;
        }
        else
        {
            del_ord=err_calc(ord,const_stps,del_t,err);
            
            ord+=del_ord;
            
            if(del_ord==0)
                ord_const_stps++;
            else
                ord_const_stps=0;
            
            ratio=step_size(del_t,ord);
            
            if(del_ord==0)
            {
                if(ratio==1.0)
                    const_stps++;
                else
                    const_stps=0;
            }
            else
                const_stps=0;
            
            
            del_t*=ratio;
            
        }
        
        step_no++;
    }
    */
    delete [] energy_stress;
}
/*--------------------------------------------
 init
 --------------------------------------------*/
TYPE0 Clock_MBDF::solve(TYPE0 del_t,int ord)
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
    curr_cost=forcefield->calc_g(1,beta,a,g);
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
    
    
    while(curr_cost>m_tol*static_cast<TYPE0>(dof_tot)
          && iter<max_iter && max_gamma>min_gamma)
    {
        
        memcpy(g0,g,dof_lcl*sizeof(TYPE0));
        memcpy(c0,c,dof_lcl*sizeof(TYPE0));
        
        gamma=1.0;
        
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
            curr_cost=forcefield->calc_g(1,beta,a,g);
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
                curr_cost=forcefield->calc_g(0,beta,a,g);
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
    
    
    TYPE0 tmp0,sum0,tot_sum0,norm;
    sum0=0.0;
    tmp1=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp1+=c_d[i]*c_d[i];
        tmp0=y_0[i]-c[i];
        sum0+=tmp0*tmp0;
        
    }
    
    MPI_Allreduce(&sum0,&tot_sum0,1,MPI_TYPE0,MPI_SUM,world);
    norm=sqrt(tot_sum0/static_cast<TYPE0>(dof_tot));
    MPI_Allreduce(&tmp1,&eq_ratio,1,MPI_TYPE0,MPI_SUM,world);
    eq_ratio=sqrt(eq_ratio/static_cast<TYPE0>(dof_tot))/e_tol;
    
    TYPE0 curr_t=t[0]+del_t;
    TYPE0 alphak=del_t/(curr_t-t[ord-1]);
    TYPE0 alpha0=0.0;
    TYPE0 alphas=0.0;
    for(int i=0;i<ord+1;i++)
    {
        alpha0+=1.0/(curr_t-t[i]);
        if(i!=0)
            alphas+=1.0/(static_cast<TYPE0>(i));
    }
    alpha0*=-del_t;
    
    err=MAX(alphak,fabs(alphak-alphas+alpha0))*norm/a_tol;
    
    return curr_cost/(m_tol*static_cast<TYPE0>(dof_tot));
}
/*--------------------------------------------
 factorial
 --------------------------------------------*/
int Clock_MBDF::fac(int no)
{
    int ans=1;
    for(int i=1;i<no+1;i++)
        ans*=i;
    return ans;
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
TYPE0 Clock_MBDF::step_size(TYPE0 del_t,int ord)
{
    
    
    TYPE0* c;
    atoms->vectors[c_n].ret(c);
    TYPE0 tmp,ratio,tot_ratio;
    TYPE0 max_h;
    
    
    ratio=INFINITY;
    
    
    for(int i=0;i<dof_lcl;i++)
    {
        if(dy[i]>0.0)
        {
            tmp=(1.0-c[i])/dy[i];
            ratio=MIN(ratio,tmp);
        }
        else if(dy[i]<0.0)
        {
            tmp=(0.0-c[i])/dy[i];
            ratio=MIN(ratio,tmp);
        }
    }
    MPI_Allreduce(&ratio,&tot_ratio,1,MPI_TYPE0,MPI_MIN,world);
    max_h=0.99*tot_ratio;
    
    TYPE0 opt_h=del_t*pow(0.5/err,1.0/static_cast<TYPE0>(ord+1));
    
    if(max_h<opt_h)
        opt_h=max_h;
    
    opt_h=MAX(min_del_t,opt_h);
    
    ratio=opt_h/del_t;
    
    if(ratio>2.0)
        ratio=1.0;
    else if (ratio<0.5)
        ratio=0.5;
    
    opt_h*=ratio;
    
    if(opt_h>=max_del_t)
    {
        opt_h=max_del_t;
        
    }
    else if(opt_h<=min_del_t)
    {
        opt_h=min_del_t;
    }
    ratio=opt_h/del_t;
    
    
    return ratio;
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
int Clock_MBDF::err_calc(int q,int const_stps,TYPE0 del_t,TYPE0& error)
{
    
    TYPE0 tmp0,tmp1;
    TYPE0 terkm2=0.0,terkm1=0.0,terk=0.0,terkp1=0.0;
    int terkm2_flag,terkm1_flag,terk_flag,terkp1_flag;
    
    terk_flag=1;
    terkp1_flag=0;
    terkm1_flag=0;
    terkm2_flag=0;
    
    if(q>2)
        terkm2_flag=1;
    if(q>1)
        terkm1_flag=1;
    
    
    if(const_stps==q+2 && q<max_order)
        terkp1_flag=1;
    
    if(terkm2_flag)
    {
        tmp1=fac(q-1)*pow(del_t,q-1);
        for(int i=0;i<q;i++)
        {
            tmp0=1.0;
            for(int j=0;j<q;j++)
            {
                if(i!=j)
                    tmp0*=t[i]-t[j];
            }
            alph_err[i]=tmp1/tmp0;
        }
        terkm2=err_est(q);
    }
    if(terkm1_flag)
    {
        tmp1=fac(q)*pow(del_t,q);
        for(int i=0;i<q+1;i++)
        {
            tmp0=1.0;
            for(int j=0;j<q+1;j++)
            {
                if(i!=j)
                    tmp0*=t[i]-t[j];
            }
            alph_err[i]=tmp1/tmp0;
        }
        terkm1=err_est(q+1);
    }
    if(terk_flag)
    {
        tmp1=fac(q+1)*pow(del_t,q+1);
        for(int i=0;i<q+2;i++)
        {
            tmp0=1.0;
            for(int j=0;j<q+2;j++)
            {
                if(i!=j)
                    tmp0*=t[i]-t[j];
            }
            alph_err[i]=tmp1/tmp0;
        }
        terk=err_est(q+2);
    }
    if(terkp1_flag)
    {
        tmp1=1.0;
        for(int i=0;i<q+3;i++)
        {
            tmp0=1.0;
            for(int j=0;j<q+3;j++)
            {
                if(i!=j)
                    tmp0*=t[i]-t[j];
            }
            alph_err[i]=tmp1/tmp0;
        }
        terkp1=err_est(q+3);
    }
    
    
    //printf("terkm2 %e terkm1 %e terk %e \n",terkm2,terkm1,terk);
    
    
    
    //lower order
    //possible if q>1
    int lower_flag=0;
    int upper_flag=0;
    
    if(q>1)
    {
        if(q>2)
        {
            if(MAX(terkm1,terkm2)<=terk)
                lower_flag=1;
        }
        else
        {
            if(terkm1<=0.5*terk)
                lower_flag=1;
            
        }
    }
    
    if(terkp1_flag)
    {
        if(q>2)
        {
            if(terkp1<terk && terk<MAX(terkm1,terkm2))
                upper_flag=1;
        }
        else if(1<q && q<3)
        {
            if(terkp1<terk && terk<terkm1)
                upper_flag=1;
        }
        else
        {
            if(terkp1<0.5*terk)
                upper_flag=1;
        }
    }
    
    if(upper_flag)
    {
        error=terkp1;
        return 1;
    }
    if(lower_flag)
    {
        error=terkm1;
        return -1;
    }
    
    error=terk;
    return 0;
    
    
}
/*--------------------------------------------
 max step size
 --------------------------------------------*/
TYPE0 Clock_MBDF::err_est(int q)
{
    TYPE0 tmp0,err_lcl,err_tot;
    err_lcl=0.0;
    for(int i=0;i<dof_lcl;i++)
    {
        tmp0=0.0;
        for(int j=0;j<q;j++)
            tmp0+=alph_err[j]*y[j][i];
        err_lcl+=tmp0*tmp0;
    }
    MPI_Allreduce(&err_lcl,&err_tot,1,MPI_TYPE0,MPI_SUM,world);
    err_tot=sqrt(err_tot/(static_cast<TYPE0>(dof_tot)*a_tol));
    return err_tot;
}


