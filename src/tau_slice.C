#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include "rt/rt.h"
#include <iostream>

using std::cout;
using std::endl;

static adios2_variable *varList1[72], *varTime;
typedef double realtype;

extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);

//======================================================================
void tau_slice(const RunData&  Run, const GridData& Grid, 
	      const PhysicsData& Physics,RTS *rts) {

  static int ini_flag = 1;

  const int iroot = 0;
  int enable_adios2 = 1;
  register int i, j, k, ind, nsl, v, ind1, node1, node2;

  int ibeg = Grid.lbeg[0], iend = Grid.lend[0];
  int jbeg = Grid.lbeg[1], jend = Grid.lend[1];
  int kbeg = Grid.lbeg[2], kend = Grid.lend[2];

  int localsize  = Grid.lsize[1]*Grid.lsize[2];
   
  float* iobuf;
  float* iosum;
  double* iobuf_double;
  double* iosum_double;

  char filename[128];

  double q1,q2;

  static int nslice; 
  static double* tau_lev;
  static int nslvar;

  FILE* fhandle=NULL;

  int define_var=1;
  //MPI_File fhandle_mpi;
  //int offset;
 
  if (ini_flag ==1){
    nslice=Physics.slice[i_sl_tau];
    tau_lev = (realtype*) malloc(nslice*sizeof(realtype));
    for (i=0;i<nslice;i++){
      tau_lev[i] = Physics.tau_lev[i];
    }
    if (Run.rank == 0) {
      cout << "tau_slice: " << nslice << endl;
      for (i=0;i<nslice;i++){
	cout << "tau_slice: " << tau_lev[i]<< endl;
      } 
    }   

    nslvar = 0;
    for (v=0;v<14;v++){
      if (Physics.tau_var[v] == 1) nslvar+=1;
    }

    ini_flag = 0;
  }

  iobuf = (float*) malloc(nslvar*localsize*sizeof(float));
  iosum = (float*) malloc(nslvar*localsize*sizeof(float));
  if (Run.enable_adios2_double){
    iobuf_double = (double*) malloc(nslvar*localsize*sizeof(double));
    iosum_double = (double*) malloc(nslvar*localsize*sizeof(double));
  }
  adios2_step_status err; 
  if ((Run.enable_adios2)&& (Run.globiter > Run.begin_iter)   && (Run.globiter > 0)){
    cout<<"tau begin step"<<endl;
    adios2_begin_step(Run.engineH, adios2_step_mode_append, 10.0, &err);
  }

  if (varList1[0] == NULL)
      define_var=0;
  for (nsl = 0; nsl<nslice; nsl++){
 
    for(v=0;v<nslvar*localsize;v++){
      iobuf[v] = 0.0;
      iosum[v] = 0.0;
      if (Run.enable_adios2_double){
        iobuf_double[v] = 0.0;
        iosum_double[v] = 0.0;
      }
    }

    for (k=kbeg; k<=kend; k++){
      for (j=jbeg; j<=jend; j++){
	ind  = j-jbeg + (k-kbeg)*Grid.lsize[1];
        for (i=ibeg; i<=iend; i++){
	  node1 = Grid.node(i,j,k);
	  node2 = Grid.node(i+1,j,k);
	  
	  if( (Grid.Tau[node1] >= tau_lev[nsl]) && (tau_lev[nsl] > Grid.Tau[node2]) ){
	    q1 = (tau_lev[nsl]-Grid.Tau[node2])/(Grid.Tau[node1]-Grid.Tau[node2]);
	    //q1 = (log(tau_lev[nsl])-log(Grid.Tau[node2]))/(log(Grid.Tau[node1])-log(Grid.Tau[node2]));
	    q2 = 1.0-q1;
	    
	    ind1 = ind;
	    
	    if (Physics.tau_var[0] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].d*q1+Grid.U[node2].d*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (Grid.U[node1].d*q1+Grid.U[node2].d*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[1] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.x*q1+Grid.U[node2].M.x*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (Grid.U[node1].M.x*q1+Grid.U[node2].M.x*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[2] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.y*q1+Grid.U[node2].M.y*q2); 
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (Grid.U[node1].M.y*q1+Grid.U[node2].M.y*q2); 
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[3] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].M.z*q1+Grid.U[node2].M.z*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (Grid.U[node1].M.z*q1+Grid.U[node2].M.z*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[4] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].e/Grid.U[node1].d*q1+Grid.U[node2].e/Grid.U[node2].d*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (Grid.U[node1].e/Grid.U[node1].d*q1+Grid.U[node2].e/Grid.U[node2].d*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[5] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.x*q1+Grid.U[node2].B.x*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (Grid.U[node1].B.x*q1+Grid.U[node2].B.x*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[6] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.y*q1+Grid.U[node2].B.y*q2);  
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (Grid.U[node1].B.y*q1+Grid.U[node2].B.y*q2);  
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[7] == 1){
	      iobuf[ind1] = (float) (Grid.U[node1].B.z*q1+Grid.U[node2].B.z*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (Grid.U[node1].B.z*q1+Grid.U[node2].B.z*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[8] == 1){
	      iobuf[ind1] = (float) (sqrt(Grid.U[node1].M.sqr())*q1+sqrt(Grid.U[node2].M.sqr())*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (sqrt(Grid.U[node1].M.sqr())*q1+sqrt(Grid.U[node2].M.sqr())*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[9] == 1){
	      iobuf[ind1] = (float) (sqrt(Grid.U[node1].B.sqr())*q1+sqrt(Grid.U[node2].B.sqr())*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = (sqrt(Grid.U[node1].B.sqr())*q1+sqrt(Grid.U[node2].B.sqr())*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[10] == 1){
	      iobuf[ind1] = (float) (Grid.temp[node1]*q1+Grid.temp[node2]*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] =  (Grid.temp[node1]*q1+Grid.temp[node2]*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[11] == 1){
	      iobuf[ind1] = (float) (Grid.pres[node1]*q1+Grid.pres[node2]*q2);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] =  (Grid.pres[node1]*q1+Grid.pres[node2]*q2);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[12] == 1){
	      iobuf[ind1] = (float) rts->Iout(j,k);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] =  rts->Iout(j,k);
	      ind1 += localsize;
	    }
	    if (Physics.tau_var[13] == 1){
	      iobuf[ind1] = float(Grid.coord(i,0)*q1+Grid.coord(i+1,0)*q2)/float(Grid.gxmax[0]);
	      if (Run.enable_adios2_double)
	        iobuf_double[ind1] = Grid.coord(i,0)*q1+Grid.coord(i+1,0)*q2/Grid.gxmax[0];
	    }
	  }
	}
      }
    }
    //cout<<"before MPI_reduce"<<endl;
    MPI_Reduce(iobuf,iosum,nslvar*localsize,MPI_FLOAT,MPI_SUM,iroot,
		  XCOL_COMM);
    if (Run.enable_adios2_double){
      MPI_Reduce(iobuf_double,iosum_double,nslvar*localsize,MPI_DOUBLE,MPI_SUM,iroot,
		  XCOL_COMM);
    }
    //cout<<"after MPI_reduce"<<endl;

    if (xcol_rank == iroot){

      /*
      sprintf(filename,"%s_%.3f.%06d","tau_slice",tau_lev[nsl],
	      Run.globiter);
      MPI_File_open(YZ_COMM,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
		    MPI_INFO_NULL,&fhandle_mpi);
      
      if(yz_rank == 0) {
	float header[4];            
	header[0] = (float) nslvar;
	header[1] = (float) Grid.gsize[1];
	header[2] = (float) Grid.gsize[2];
	header[3] = (float) Run.time;
	MPI_File_write(fhandle_mpi,header,4,MPI_FLOAT,MPI_STATUS_IGNORE);
      }
      
      offset = 4*sizeof(float);
      
      MPI_File_set_view(fhandle_mpi,offset,MPI_FLOAT,yz_subarray,"native",
			MPI_INFO_NULL);
      MPI_File_write_all(fhandle_mpi,&(iosum[0]),nslvar*localsize,
			 MPI_FLOAT,MPI_STATUS_IGNORE);
      
      MPI_File_close(&fhandle_mpi);   
      */

      char adios2_var_name[128];
      int curr=0;
      if(Physics.slice[i_sl_collect] == 0) {	
        if(Run.enable_adios2){
          if ((Run.globiter > Run.begin_iter && define_var==0) || Run.globiter < std::min(Run.resfreq,Run.slicefreq)){
            sprintf(adios2_var_name,"%s_%.3f.%06d","tau_slice",
                  tau_lev[nsl],curr);
            uint64_t gdim[2],start[2],count[2];
            uint64_t fixed_shape[1], fixed_start[1], fixed_count[1];
            int vi, ci, v_index;
            fixed_shape[0]=1;
            fixed_start[0]=0;
            fixed_count[0]=1;
            for(ci=0;ci<2; ci++){
              gdim[ci] = Grid.gsize[ci+1];
              start[ci] = Grid.beg[ci+1]-Grid.gbeg[ci+1];
              count[ci] = Grid.lsize[ci+1];
            }
            //cout<<"tau slice name "<<adios2_var_name<<" "<<start[0]<<" "<<start[1]<<" "<<count[0]<<" "<<count[1]<<endl;
            for (vi=0;vi<nslvar;vi++){
                v_index=vi+nslvar*nsl;
	        if(tau_lev[nsl] >= 1e-3){
                  sprintf(adios2_var_name,"%s_%.3f.%d.%06d","tau_slice",
                      tau_lev[nsl],vi,curr);
                }
                else{
                  sprintf(adios2_var_name,"%s_%.6f.%d.%06d","tau_slice",
                      tau_lev[nsl],vi,curr);
                }
                //cout<<"tau slice name "<<adios2_var_name<<" "<<tau_lev[nsl]<<" "<<nsl<<endl;
                if(Run.enable_adios2_double) {
                  varList1[v_index] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_double,2, gdim,
                               start, count, adios2_constant_dims_true);
                }
                else{
                  varList1[v_index] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_float,2, gdim,
                               start, count, adios2_constant_dims_true);
                }
            }
            if (varTime==NULL){
                varTime = adios2_define_variable(Run.ioH, "Run.time.2d",adios2_type_float,0,NULL,
                           NULL,NULL,adios2_constant_dims_true);
            }
          }
        } 
        if ((Run.enable_adios2)&& (Run.globiter > Run.begin_iter)  &&( Run.globiter >0)){
          int vi,v_index;
          if (nsl==0){
              float temp=(float)Run.time;
              adios2_put(Run.engineH,varTime,&temp,adios2_mode_sync);
          }
          //cout<<"put tau slice "<<Run.maxiter-Run.globiter<<" "<<std::min(Run.resfreq,Run.slicefreq)<<endl;
          //cout<<Run.rank<<" tau begin step "<<yz_rank<<" xcol "<<xcol_rank<<endl;
          for(vi=0;vi<nslvar;vi++){
              v_index=vi+nslvar*nsl;
              if (Run.enable_adios2_double)
                adios2_put(Run.engineH,varList1[v_index],&iosum_double[localsize*vi],adios2_mode_sync);
              else
                adios2_put(Run.engineH,varList1[v_index],&iosum[localsize*vi],adios2_mode_sync);
          }
        }
        if((Run.globiter == 0) || (Run.maxiter-Run.globiter)<std::min(Run.resfreq,Run.slicefreq) || !Run.enable_adios2) {
	if(yz_rank == 0) {
	  if(tau_lev[nsl] >= 1e-3)
	    sprintf(filename,"%s%s_%.3f.%06d",Run.path_2D,"tau_slice",
		    tau_lev[nsl],Run.globiter);
	  else
	    sprintf(filename,"%s%s_%.6f.%06d",Run.path_2D,"tau_slice",
		    tau_lev[nsl],Run.globiter);
	  fhandle=fopen(filename,"w");
	  
	  float header[4];            
	  header[0] = (float) nslvar;
	  header[1] = (float) Grid.gsize[1];
	  header[2] = (float) Grid.gsize[2];
	  header[3] = (float) Run.time;
	  fwrite(header,sizeof(float),4,fhandle);
	}
        //cout<<Run.rank<<" before slice_write"<<endl;	
	slice_write(Grid,0,&(iosum[0]),localsize,nslvar,1,2,fhandle);
	
          //std::cout<<localsize<<" "<<nslvar<<" "<<Run.rank<<" tau filename 1 "<<std::endl;
	if(yz_rank == 0)
	  fclose(fhandle);
      }	
      } else {
	if(yz_rank == 0){
	  sprintf(filename,"%s_%.3f.dat","tau_slice",tau_lev[nsl]);
	  fhandle=fopen(filename,"a");
	}
          //std::cout<<"filename 2 "<<filename<<std::endl;
	
	slice_write(Grid,0,&(iosum[0]),localsize,nslvar,1,2,fhandle);
	
	if(yz_rank == 0){ 
	  fclose(fhandle);
	  
	  sprintf(filename,"%s_%.3f.log","tau_slice",tau_lev[nsl]);
          //std::cout<<"filename 3 "<<filename<<std::endl;
	  
      std::fstream fptr;
	  int newfile = 0;
	  fptr.open(filename,std::ios::in);
	  if (!fptr) newfile = 1;
	  fptr.close();
	  
	  fptr.open(filename,std::ios::out|std::ios::app);
	  fptr.precision(10);
	  if (newfile) {      
	    fptr <<  nslvar << ' ' <<  Grid.gsize[1] << ' ' 
		 << Grid.gsize[2] << endl;
	    fptr << Physics.tau_var[0]  << ' ' 
		 << Physics.tau_var[1]  << ' ' 
		 << Physics.tau_var[2]  << ' ' 
		 << Physics.tau_var[3]  << ' ' 
		 << Physics.tau_var[4]  << ' ' 
		 << Physics.tau_var[5]  << ' ' 
		 << Physics.tau_var[6]  << ' ' 
		 << Physics.tau_var[7]  << ' ' 
		 << Physics.tau_var[8]  << ' ' 
		 << Physics.tau_var[9]  << ' ' 
		 << Physics.tau_var[10] << ' ' 
		 << Physics.tau_var[11] << ' '
		 << Physics.tau_var[12] << ' '        
		 << Physics.tau_var[13] << endl;
	  }
	  fptr << Run.globiter << ' ' << Run.time << endl;
	  fptr.close();
	}
      }
    }
  }
  if ((Run.enable_adios2) && (Run.globiter > Run.begin_iter)  && (Run.globiter >0 )){
    adios2_end_step(Run.engineH);
  }

  free(iobuf);
  free(iosum);
  if (Run.enable_adios2_double){
  free(iobuf_double);
  free(iosum_double);
  }
}

