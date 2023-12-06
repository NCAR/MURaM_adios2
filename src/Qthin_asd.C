#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "run.H"
#include "grid.H"
#include "comm_split.H"
#include <iostream>
#include <algorithm>

#define YZ_LOOP(G,j,k) \
for((k)=(G).lbeg[2];(k)<=(G).lend[2];(k)++) \
for((j)=(G).lbeg[1];(j)<=(G).lend[1];(j)++)

using std::min;
using std::max;
using std::cout;
using std::endl;

void IntegrateColumn(const GridData &Grid, const RunData &Run, const double * Array, double * result, int zl, int zh);
double spline_1d(double *X, double *Y, double X0, int Nx);
double lin_1d(double *X, double *Y, double X0, int Nx);

// N are the table sizes
// the axes temperature (tax), column depth (cmassax), hydrogen column depth (tauhax)
// And the tables,
// Neutral fractions (_nf), excitaiton (_exc) and Escape probabilities (_esc)
// for Hydrogen, Calcium and Magnesium
int N_T;
int N_Hesc, N_Caesc, N_Mgesc;
double *tax;
double *C_thin;
double *H_nf;
double *H_nf_exc;
double *Ca_nf_exc;
double *Mg_nf_exc;

double *tauhax_Hesc;
double *H_esc;

double *cmassax_Caesc;
double *Ca_esc;

double *cmassax_Mgesc;
double *Mg_esc;

double dT_Tab;

// Abundances relative to hydrgoen from Asplund 2009
const double mp   = 1.67262158E-24;
/* Number fraction Abundance relative to hydrogen from
 * a 15 element mix using Asplund 2009 abundances */

const double A_Mg = 3.9810e-5*1.00784/24.305;
const double A_Ca = 2.1877e-6*1.00784/40.078;

// Electron and hydrogen number density from simple H/He mix.
const double XX = 0.7;
const double nh = XX/mp; 

void ELTE_init () {

  // Read in all the required tables. For now, I think this is enough.

  std::string strbuffer;

  std::ifstream fptr("./radtab.dat",std::ios::in);

  if(fptr){
    fptr >> strbuffer;
    fptr >> N_T;

    tax = new double [N_T]();
    C_thin = new double [N_T]();
    H_nf = new double [N_T]();
    H_nf_exc = new double [N_T]();
    Ca_nf_exc = new double [N_T]();
    Mg_nf_exc = new double [N_T]();

    double temp = 0;

    for(int nn = 0; nn<N_T;nn++){
      fptr >> tax[nn];
    }

    dT_Tab = 1.0/(tax[1]-tax[0]);

    fptr >> strbuffer;
    for(int nn = 0; nn<N_T;nn++){
      fptr >> C_thin[nn];
    }
    fptr >> strbuffer;
    for(int nn = 0; nn<N_T;nn++){
      fptr >> H_nf[nn];
      H_nf_exc[nn] = H_nf[nn];
    }
    fptr >> strbuffer;
    for(int nn = 0; nn<N_T;nn++){
      fptr >> temp;
      H_nf_exc[nn]*=temp;
    }
    fptr >> strbuffer;
    for(int nn = 0; nn<N_T;nn++){
      fptr >> Mg_nf_exc[nn];
    }
    fptr >> strbuffer;
    for(int nn = 0; nn<N_T;nn++){
      fptr >> temp;
      Mg_nf_exc[nn]*=temp;
    }
    fptr >> strbuffer;
    for(int nn = 0; nn<N_T;nn++){
      fptr >> Ca_nf_exc[nn];
    }
    fptr >> strbuffer;
    for(int nn = 0; nn<N_T;nn++){
      fptr >> temp;
      Ca_nf_exc[nn]*=temp;
    }

    fptr >> strbuffer;
    fptr >> N_Hesc;

    tauhax_Hesc = new double [N_Hesc]();
    H_esc = new double [N_Hesc]();

    for(int nn = 0; nn<N_Hesc;nn++){
      fptr >> tauhax_Hesc[nn];
    }

    fptr >> strbuffer;

    for(int nn = 0; nn<N_Hesc;nn++){
      fptr >> H_esc[nn];
    }

    fptr >> strbuffer;
    fptr >> N_Mgesc;

    cmassax_Mgesc = new double [N_Mgesc]();
    Mg_esc = new double [N_Mgesc]();

    for(int nn = 0; nn<N_Mgesc;nn++){
      fptr >> cmassax_Mgesc[nn];
    }

    fptr >> strbuffer;
    for(int nn = 0; nn<N_Mgesc;nn++){
      fptr >> Mg_esc[nn];
    }

    fptr >> strbuffer;
    fptr >> N_Caesc;

    cmassax_Caesc = new double [N_Caesc]();
    Ca_esc = new double [N_Caesc]();

    for(int nn = 0; nn<N_Caesc;nn++){
      fptr >> cmassax_Caesc[nn];
    }

    fptr >> strbuffer;
    for(int nn = 0; nn<N_Caesc;nn++){
      fptr >> Ca_esc[nn];
    }
  } else {
    fprintf(stdout,"radtab.dat not found. Aborting \n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

}

/*****************************************************************************/
void Get_Radloss(const RunData&  Run, GridData& Grid,const PhysicsData& Physics){

  static int radloss_ini_flag = 1;
  
  register int i,j,k,off,node,n,n1,n2;

  const int i_beg    = Grid.lbeg[0];
  const int i_end    = Grid.lend[0];
/*
  int next[3];
  for(i=0;i<3;i++)
    next[i] = Grid.stride[i];
*/

  //int next0 = Grid.stride[0];
  int next1 = Grid.stride[1];
  int next2 = Grid.stride[2];

  double t1,t2,r1,r2,ff,pt,pr,ts,rr,qq,qloss,t_a,t_b,tmin,tmax;
    
  // Chianti assumes Q=-n_e*n_H*L(T)
  const double X_H    = 0.726;  // Asplund09
  const double mp     = 1.67262158E-24;
  const double ne_ref = 1e9;
  const double ne_par = (sqrt(0.5*(1.0+X_H)*X_H)/mp)/ne_ref;

  const double inv_pmax = 1.0/Physics.rt[i_rt_pre_cut];
  const double p_thresh = Physics.rt[i_rt_tr_pre];
  const double T_thresh = Physics.rt[i_rt_tr_tem];


  static int ntab=0;

  static double* T_tab;
  static double* Q_tab;
  
  static double lgT0,del0;
 
  if(radloss_ini_flag){
    std::ifstream fptr("Radloss_Chianti.dat",std::ios::in);
    if(fptr){
      fptr.precision(16);
      fptr >> ntab;
      T_tab = new double[ntab];
      Q_tab = new double[ntab];
      for(i=0;i<ntab;i++)
	fptr >> T_tab[i] >> Q_tab[i];
      fptr.close();
    } else {
      cout << "Radloss Table not found .... abort" << endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    lgT0=log(T_tab[0]);
    del0=1.0/(log(T_tab[1])-log(T_tab[0]));

    for(i=0;i<ntab;i++){
      T_tab[i]  = log(T_tab[i]);
      Q_tab[i] *= ne_ref*ne_ref;
     }

    if(Run.rank == 0){
      cout << "Radloss: use log(T), log(rho)" << endl;
      cout << "Radloss: use X_H = " << X_H << endl;
      //cout << "Radloss: pressure_cutoff = " << 1./inv_pmax << endl;
      cout << "Radloss: use tr_switch" << endl;
    }
   
    radloss_ini_flag = 0;

#pragma acc enter data copyin(T_tab[:ntab])
#pragma acc enter data copyin(Q_tab[:ntab])
  }

  const int ibeg = Grid.lbeg[0];
  const int iend = Grid.lend[0];
  const int jbeg = Grid.lbeg[1];
  const int jend = Grid.lend[1];
  const int kbeg = Grid.lbeg[2];
  const int kend = Grid.lend[2];
  const int bufsize = Grid.bufsize;
  double qthin[2][i_end+2];
  double pswitch,tswitch;

#pragma acc parallel loop collapse(2) gang \
 present(Grid[:1], Grid.Qthin[:bufsize], Grid.temp[:bufsize], \
         Grid.pres[:bufsize], Grid.U[:bufsize], \
         T_tab[:ntab], Q_tab[:ntab]) \
 private(off, qthin[:2][:i_end+2],tswitch,pswitch)
  for(k=kbeg; k<=kend; k++)
  for(j=jbeg; j<=jend; j++) {
//#pragma acc cache(qthin[:2][:i_end+2])
    off = j*next1+k*next2;

#pragma acc loop vector
    for(i=i_beg-1;i<=i_end+1;i++) {
      qthin[0][i] = 0.0;
      qthin[1][i] = 0.0;
      //Grid.Qthin[off+i]=0.0;
    }

    //node = off+i_beg-1;
    //t1 = log(Grid.temp[node]);
    //r1 = log(Grid.U[node].d*ne_par);
    
#pragma acc loop vector private(node, t1, t2, r1, r2, n1, n2, \
 t_a, t_b, ff, ts, pr, pt, rr, qq, qloss)
    for(i=i_beg;i<=i_end+1;i++){
      node = off+i;
     
      t1 = log(Grid.temp[node-1]);
      r1 = log(Grid.U[node-1].d*ne_par);

      t2 = log(Grid.temp[node]);
      r2 = log(Grid.U[node].d*ne_par);

      tmin = min(t1,t2);
      tmax = max(t1,t2);

      n1 = (int) ( (tmin-lgT0)*del0 );
      n2 = (int) ( (tmax-lgT0)*del0 );

      n1 = max(0,n1);
      n2 = min(ntab-2,n2);
      
#pragma acc loop seq
      for(n=n1;n<=n2;n++){
	
	t_a=max(tmin,T_tab[n]);
	t_b=min(tmax,T_tab[n+1]);
	
	if(t_b-t_a > 1e-6 ){
	  ff = (t_b-t_a)/(tmax-tmin);
	  ts = 0.5*(t_a+t_b);
	  pr = (t2-ts)/(t2-t1);
	  pt = (T_tab[n+1]-ts)/(T_tab[n+1]-T_tab[n]);
	} else if (t_b-t_a >= 0.) {
	  ff = 1.0;
	  ts = 0.5*(t_a+t_b);
	  pr = 0.5;
	  pt = (T_tab[n+1]-ts)/(T_tab[n+1]-T_tab[n]);
	} else {
	  ff=0.0;
	  pr=0.5;
	  pt=0.5;
	}
	
	rr = exp(pr*r1+(1.0-pr)*r2);
	qq = pt*Q_tab[n]+(1.0-pt)*Q_tab[n+1];
	
	qloss = -rr*rr*qq*ff;

	//Grid.Qthin[node-1] += qloss*pr;
	//Grid.Qthin[node]   += qloss*(1.0-pr);
	qthin[0][i-1] += qloss*pr;
        qthin[1][i] += qloss*(1.0-pr);
      }
      //t1 = t2;
      //r1 = r2;
    }

#pragma acc loop vector
    for(i=i_beg-1; i<=i_end+1; i++)
      Grid.Qthin[off+i] = qthin[0][i] + qthin[1][i];

    // remove radiative loss in high pressure regions 
#pragma acc loop vector private(pswitch,tswitch)
    for(i=i_beg;i<=i_end;i++){
      //Grid.Qthin[off+i] *= max(0.0,1.-pow(Grid.pres[off+i]*inv_pmax,2));
      // *** new ***
      pswitch = min(max(p_thresh-Grid.pres[off+i],0.0),1.0);
      tswitch = min(max(Grid.temp[off+i]-T_thresh,0.0),1.0);
      Grid.Qthin[off+i] *= min(pswitch,tswitch);
      // *** new ***
    }
  }

}

void ELTE_decon(){

  delete [] C_thin;
  delete [] tax;
  delete [] H_nf;
  delete [] H_nf_exc;
  delete [] tauhax_Hesc;
  delete [] H_esc;
  delete [] Ca_nf_exc;
  delete [] cmassax_Caesc;
  delete [] Ca_esc;
  delete [] Mg_nf_exc;
  delete [] cmassax_Mgesc;
  delete [] Mg_esc;
}

void ELTE_Qx(const RunData&  Run, GridData& Grid, const PhysicsData& Physics){

// Calculate Chromospheric heating/cooling from NLTE using empirical tables of
// C&L2012

  const int i_beg = Grid.lbeg[0];
  const int i_end = Grid.lend[0];
  const int v_length = i_end-i_beg+1; 
  double H_neutral[v_length], density[v_length], cmass[v_length], Hcol[v_length];
  double ff, ts, pr, pt;
  int i,j,k;

  const double inv_pmax = 1.0/Physics.rt[i_rt_pre_cut];

  YZ_LOOP(Grid,j,k){
    for(i=0;i<v_length;i++){
      int inode = Grid.node(i+i_beg,j,k);
      density[i] = Grid.U[inode].d;
      double nf = lin_1d(tax,H_nf,log(Grid.temp[inode]),N_T);
      H_neutral[i] = nh*density[i]*nf;
    }

    for(i=i_beg-1;i<=i_end+1;i++){
      int inode = Grid.node(i,j,k);
      Grid.QH[inode] =  0.0;
      Grid.QMg[inode] = 0.0;
      Grid.QCa[inode] =  0.0;
      Grid.Qthin[inode] =  0.0;
    }
    
    // Integrate for column mass and neutral hydrogen column mass
    IntegrateColumn(Grid, Run, H_neutral,Hcol,0,v_length);
    IntegrateColumn(Grid, Run, density,cmass,0,v_length);
      
    // Overlap Interval for NLTE line cooling (all on same Tax)
    
    double lgT0 = tax[0];
 
    int inode = Grid.node(i_beg-1,j,k);
    double lgT1 = log(Grid.temp[inode]);
    double r1 = log(sqrt(Grid.ne[inode]*Grid.U[inode].d*nh));
    
    for (i = i_beg ; i <= i_end+1 ; i++){
      inode = Grid.node(i-1,j,k);
      
      double lgT2 = log(Grid.temp[inode+1]);
      double r2 = log(sqrt(Grid.ne[inode+1]*Grid.U[inode+1].d*nh)); 
      
      double tmin = min(lgT1,lgT2);
      double tmax = max(lgT1,lgT2);
      
      int n1 = (int) ( (tmin-lgT0)*dT_Tab );
      int n2 = (int) ( (tmax-lgT0)*dT_Tab );
      
      n1 = max(0,n1);
      n2 = min(N_T-2,n2);
      
      for(int n=n1;n<=n2;n++){

        double t_a=max(tmin,tax[n]);
        double t_b=min(tmax,tax[n+1]);

        if(t_b-t_a > 1e-6){
          ff = (t_b-t_a)/(tmax-tmin);
          ts = 0.5*(t_a+t_b);
          pr = (lgT2-ts)/(lgT2-lgT1);
          pt = (tax[n+1]-ts)/(tax[n+1]-tax[n]);
        } else if (t_b-t_a >= 0.0) {
          ff = 1.0;
          ts = t_a;
          pr = 0.5;
          pt = (tax[n+1]-ts)/(tax[n+1]-tax[n]);
        } else {
          ff=0.0;
          pr=0.5;
          pt=0.5;
        }

        double rr = exp(pr*r1+(1.0-pr)*r2);

        double qq = pt*H_nf_exc[n]+(1.0-pt)*H_nf_exc[n+1];
        double qloss = -rr*rr*qq*ff;

        Grid.QH[inode] += qloss*pr;
        Grid.QH[inode+1]+= qloss*(1.0-pr);

        qq = pt*Ca_nf_exc[n]+(1.0-pt)*Ca_nf_exc[n+1];
        qloss = -rr*rr*qq*ff;

        Grid.QCa[inode] += qloss*pr;
        Grid.QCa[inode+1] += qloss*(1.0-pr);

        qq = pt*Mg_nf_exc[n]+(1.0-pt)*Mg_nf_exc[n+1];
        qloss = -rr*rr*qq*ff;

        Grid.QMg[inode] += qloss*pr;
        Grid.QMg[inode+1] += qloss*(1.0-pr);

        qq = pt*C_thin[n]+(1.0-pt)*C_thin[n+1];
        qloss = -rr*rr*qq*ff;

        Grid.Qthin[inode] += qloss*pr;
        Grid.Qthin[inode+1] += qloss*(1.0-pr);
      }
      lgT1 = lgT2;
      r1 = r2;
    }

    for (i = i_beg ; i <= i_end ; i++){

      inode = Grid.node(i,j,k);

      // pressure cut off 
      //Grid.Qthin[inode] *= max(0.0,1.-pow(Grid.pres[inode]*inv_pmax,2));
      // *** new ***
      double pswitch = min(max(Physics.rt[i_rt_tr_pre]-Grid.pres[inode],0.0),1.0);
      double tswitch = min(max(Grid.temp[inode]-Physics.rt[i_rt_tr_tem],0.0),1.0);
      Grid.Qthin[inode] *= min(pswitch,tswitch);
      // *** new ***
      
      // Multiple NLTE line losses by escape probability
      double escape = lin_1d(tauhax_Hesc,H_esc,log(Hcol[i-i_beg]*4.0e-14),N_Hesc);
      Grid.QH[inode]  *= escape;

      escape = lin_1d(cmassax_Caesc,Ca_esc,log(cmass[i-i_beg]),N_Caesc);
      Grid.QCa[inode] *= escape*A_Ca;

      escape = lin_1d( cmassax_Mgesc,Mg_esc,log(cmass[i-i_beg]),N_Mgesc);
      Grid.QMg[inode] *= escape*A_Mg;
    }
  }
}

void IntegrateColumn(const GridData &Grid, const RunData &Run, const double * Array, double * result, int zl, int zh) {
  
  double dx = Grid.dx[0];
  double temp_int;
  
  double * temp_int_all = new double [Grid.procs[0]]();
  double * temp_int_sum = new double [Grid.procs[0]]();
  
  //top point
  //Set to something very low but not numerically 0.

  if (Grid.is_gend[0]){
    result[zh-1] = 1.0e-20;
    temp_int = 0.5*dx*Array[zh-1];
  }  else  {
    result[zh-1] = 0.5*dx*Array[zh-1];
    temp_int = dx*Array[zh-1];
  }
  // integrate
  for (int kk = zh-2; kk >= zl ; --kk){
    result[kk] = temp_int + Array[kk]*0.5*dx;
    temp_int += dx*Array[kk];
  }

  // To share in only one function call put temp_int in another array
  // if I am located above that core
  for (int nn=0;nn < Grid.procs[0];nn++){
    if (nn < Run.zrank){
      temp_int_all[nn] = temp_int;
    } 
  }

  // SUM contributions for different cores
  MPI_Allreduce(&temp_int_all[0],&temp_int_sum[0], Grid.procs[0],MPI_DOUBLE,MPI_SUM,XCOL_COMM);

  // Compute final sum by adding contribution of all cores above me
  for (int kk = zl; kk < zh ; ++kk){
    result[kk] += temp_int_sum[Run.zrank];
  }

  delete [] temp_int_all;
  delete [] temp_int_sum;

  return;
}

double spline_1d(double *X, double *Y, double X0, int Nx) {
  
  double Y0;

  if (X0 <= X[0]){
    // if less than the table axis set to first point
    Y0 = Y[0];
  } else if (X0 >= X[Nx-1]){
    // if greater than the table axis set to last point
    Y0 = Y[Nx-1];
  } else if (X0 < X[1]) {
    // if within the first cell, interpolate

  double c0 = Y[0];
  double c1 = Y[1];
  double c2 = Y[2];
  double c3 = Y[3];
  double d1 = X[0] - X[1];
  double d2 = X[2] - X[1];
  double d3 = X[3] - X[1];

  double b = ((c0-c1+(c1-c3)*d1/d3)/(d1*d1*d1-d3*d3*d1)-(c3-c1+(c1-c2)*d3/d2)/(d3*d3*d3-d2*d2*d3))*(d1*d2+d1*d3+d3*d2+d3*d3)/(d2-d1);
  double a = ((c0-c1+(c1-c3)*d1/d3)-b*(d1*d1-d3*d1))/(d1*d1*d1-d3*d3*d1);
  double c = (c0-c1)/d1 -a*d1*d1-b*d1;
  double d = c1;

  double delx = X0-X[1];
  Y0 = a*delx*delx*delx+b*delx*delx+c*delx+d;

  } else if (X0 > X[Nx-2]){
   // if within the last cell, interpolate  

  double c0 = Y[Nx-4];
  double c1 = Y[Nx-3];
  double c2 = Y[Nx-2];
  double c3 = Y[Nx-1];
  double d1 = X[Nx-4] - X[Nx-3];
  double d2 = X[Nx-2] - X[Nx-3];
  double d3 = X[Nx-1] - X[Nx-3];

  double b = ((c0-c1+(c1-c3)*d1/d3)/(d1*d1*d1-d3*d3*d1)-(c3-c1+(c1-c2)*d3/d2)/(d3*d3*d3-d2*d2*d3))*(d1*d2+d1*d3+d3*d2+d3*d3)/(d2-d1);
  double a = ((c0-c1+(c1-c3)*d1/d3)-b*d1*(d1-d3))/(d1*d1*d1-d3*d3*d1);
  double c = (c0-c1)/d1 -a*d1*d1-b*d1;
  double d = c1;

  double delx = X0-X[Nx-3];
  Y0 = a*delx*delx*delx+b*delx*delx+c*delx+d;

  } else {
   // if anywhere else, find and interpolate  

    int nn=0;
    while(X[nn] < X0){
      nn++;
    }
    if (X[nn] == X0){
    // if exactly equal to an X-index
    Y0 = Y[nn];
    } else {
    // else set up interpolant
    double c0 = Y[nn-2];
    double c1 = Y[nn-1];
    double c2 = Y[nn];
    double c3 = Y[nn+1];
    double d1 = X[nn-2] - X[nn-1];
    double d2 = X[nn] - X[nn-1];
    double d3 = X[nn+1] - X[nn-1];

    double b = ((c0-c1+(c1-c3)*d1/d3)/(d1*d1*d1-d3*d3*d1)-(c3-c1+(c1-c2)*d3/d2)/(d3*d3*d3-d2*d2*d3))*(d1*d2+d1*d3+d3*d2+d3*d3)/(d2-d1);
    double a = ((c0-c1+(c1-c3)*d1/d3)-b*(d1*d1-d3*d1))/(d1*d1*d1-d3*d3*d1);
    double c = (c0-c1)/d1 -a*d1*d1-b*d1;
    double d = c1;

    double delx = X0-X[nn-1];
    Y0 = a*delx*delx*delx+b*delx*delx+c*delx+d; 
    }
  }


  return Y0;
  }

double lin_1d(double *X, double *Y, double X0, int Nx) {
  
  double Y0;

  if (X0 <= X[0]){
    // if less than the table axis set to first point
    Y0 = Y[0];
  } else if (X0 >= X[Nx-1]){
    // if greater than the table axis set to last point
    Y0 = Y[Nx-1];
  } else {
   // if anywhere else, find and interpolate  
    int nn=0;
    while(X[nn] < X0){
      nn++;
    }
    if (X[nn] == X0){
    // if exactly equal to an X-index
    Y0 = Y[nn];
    } else {
    // else set up interpolant
    Y0 = Y[nn-1]+(Y[nn]-Y[nn-1])/(X[nn]-X[nn-1])*(X0-X[nn-1]);
    }
  }
  return Y0;
}

