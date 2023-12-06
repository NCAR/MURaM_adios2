#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include <iostream>

using std::cout;
using std::endl;

static adios2_variable *varList1[1];
extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);

void xy_slice(const RunData&  Run, const GridData& Grid, 
	      const PhysicsData& Physics) {

  static int ini_flag = 1;

  register int i, j, k, node, ind, nsl, v;

  int ibeg = Grid.lbeg[0], iend = Grid.lend[0];
  int jbeg = Grid.lbeg[1], jend = Grid.lend[1];

  int localsize  = Grid.lsize[0]*Grid.lsize[1];
   
  float* iobuf;
  double* iobuf_double;

  char filename[128];

  static int nslice; 
  static int* ixpos;
  static int nslvar;

  FILE* fhandle=NULL;
  
  int define_var=1;
  //MPI_File fhandle_mpi;
  //int offset;
  
  if (ini_flag ==1){
    nslice=Physics.slice[i_sl_xy];
    ixpos = (int*) malloc(nslice*sizeof(int));
    for (i=0;i<nslice;i++){
      ixpos[i] = Physics.xy_lev[i];
    }

    if (Run.rank == 0) {
      if(Run.verbose > 1)cout << "xy_slice: " << nslice << endl;
      for (i=0;i<nslice;i++){
      if(Run.verbose > 1) cout << "xy_slice: " << ixpos[i]<< endl;
      }     
    }

    nslvar = 0;
    for (v=0;v<12;v++){
      if (Physics.xy_var[v] == 1) nslvar+=1;
    }

    ini_flag = 0;
  }

  iobuf = (float*) malloc(nslvar*localsize*sizeof(float));
  if (Run.enable_adios2_double)
    iobuf_double = (double*) malloc(nslvar*localsize*sizeof(double));
 
  int var_count, adios2_ind;
  var_count = 0;
  adios2_ind = 0;
  adios2_step_status err;
  if ((Run.enable_adios2) && (Run.globiter > Run.begin_iter) && (Run.globiter > 0 )){
    cout<<"xy slice begin step"<<endl;
    adios2_begin_step(Run.engineH, adios2_step_mode_append, 10.0, &err);
  }

  if (varList1[0] == NULL)
      define_var=0;
  for (nsl = 0; nsl<nslice; nsl++){

    if ( (Grid.beg[2] <= ixpos[nsl]+Grid.gbeg[2] ) and 
         (Grid.end[2] >= ixpos[nsl]+Grid.gbeg[2] )){

      for (i=ibeg; i<=iend; i++){
        for (j=jbeg; j<=jend; j++){
          ind  = j-jbeg + (i-ibeg)*Grid.lsize[1];
          k    = Grid.lbeg[2]+ixpos[nsl]+Grid.gbeg[2]-Grid.beg[2];
          node = Grid.node(i,j,k);

	  if (Physics.xy_var[0] == 1){
	    iobuf[ind] = (float) Grid.U[node].d;
            if (Run.enable_adios2_double)
	      iobuf_double[ind] = Grid.U[node].d;
	    ind += localsize;
	  }
	  if (Physics.xy_var[1] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.x;
            if (Run.enable_adios2_double)
	      iobuf_double[ind] = Grid.U[node].M.x;
	    ind += localsize;
	  }
	  if (Physics.xy_var[2] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.y; 
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = Grid.U[node].M.y; 
	    ind += localsize;
	  }
	  if (Physics.xy_var[3] == 1){
	    iobuf[ind] = (float) Grid.U[node].M.z;
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = Grid.U[node].M.z;
	    ind += localsize;
	  }
	  if (Physics.xy_var[4] == 1){
	    iobuf[ind] = (float) (Grid.U[node].e/Grid.U[node].d);
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = (Grid.U[node].e/Grid.U[node].d);
	    ind += localsize;
	  }
	  if (Physics.xy_var[5] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.x;
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = Grid.U[node].B.x;
	    ind += localsize;
	  }
	  if (Physics.xy_var[6] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.y;  
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = Grid.U[node].B.y;  
	    ind += localsize;
	  }
	  if (Physics.xy_var[7] == 1){
	    iobuf[ind] = (float) Grid.U[node].B.z;
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = Grid.U[node].B.z;
	    ind += localsize;
	  }
	  if (Physics.xy_var[8] == 1){
	    iobuf[ind] = (float) sqrt(Grid.U[node].M.sqr());
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = sqrt(Grid.U[node].M.sqr());
	    ind += localsize;
	  }
	  if (Physics.xy_var[9] == 1){
	    iobuf[ind] = (float) sqrt(Grid.U[node].B.sqr());
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = sqrt(Grid.U[node].B.sqr());
	    ind += localsize;
	  }
	  if (Physics.xy_var[10] == 1){
	    iobuf[ind] = (float) Grid.temp[node];
            if (Run.enable_adios2_double)
	    iobuf_double[ind] = Grid.temp[node];
	    ind += localsize;
	  }
	  if (Physics.xy_var[11] == 1){
	    iobuf[ind] = (float) Grid.pres[node];
            if (Run.enable_adios2_double)
	    iobuf_double[ind] =  Grid.pres[node];
	  }
	}
      }

      char adios2_var_name[128]; 
      int curr=0;
      if(Physics.slice[i_sl_collect] == 0) {
        if(Run.enable_adios2){
          if (Run.globiter < std::min(Run.resfreq,Run.slicefreq) || (Run.globiter > Run.begin_iter && define_var==0) ){
            sprintf(adios2_var_name,"%s_%04d.%06d","xy_slice",
                  ixpos[nsl],curr);
            uint64_t gdim[3],start[3],count[3];
            gdim[0]=nslvar;
            start[0]=0;
            count[0]=nslvar;;
            for(int i=1;i<3; i++){
              gdim[i]=Grid.gsize[i-1];
              start[i] = Grid.beg[i-1]-Grid.gbeg[i-1];
              count[i] = Grid.lsize[i-1];
            }
            //cout<<"xy slice name "<<adios2_var_name<<" "<<localsize<<" "<<start[0]<<" "<<start[1]<<" "<<start[2]<<" "<<count[0]<<" "<<count[1]<<" "<<count[2]<<endl;

            if(Run.enable_adios2_double) {
              varList1[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_double,3, gdim,
                           start, count, adios2_constant_dims_true);
            }
            else{
              varList1[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_float,3, gdim,
                           start, count, adios2_constant_dims_true);
            }
            var_count++;
          }
        }
        //if (Run.enable_adios2){
        if (Run.enable_adios2 && (Run.globiter > Run.begin_iter) && Run.globiter >0){
          //cout<<" xy_slice put "<<Run.rank<<" "<<xy_rank<<endl;
          if (Run.enable_adios2_double)
            adios2_put(Run.engineH,varList1[adios2_ind],&iobuf_double[0],adios2_mode_sync);
          else
            adios2_put(Run.engineH,varList1[adios2_ind],&iobuf[0],adios2_mode_sync);
          adios2_ind++;
        }
        
        if(!Run.globiter || !Run.enable_adios2 ||(Run.maxiter-Run.globiter)<std::min(Run.resfreq,Run.slicefreq)){
	if(xy_rank == 0) { 
	  sprintf(filename,"%s%s_%04d.%06d",Run.path_2D,"xy_slice",ixpos[nsl],
		  Run.globiter);
	    fhandle=fopen(filename,"w");
	    
	    float header[4];            
	    header[0] = (float) nslvar;
	    header[1] = (float) Grid.gsize[1];
	    header[2] = (float) Grid.gsize[0];
	    header[3] = (float) Run.time;
	    fwrite(header,sizeof(float),4,fhandle);
	}
	
	slice_write(Grid,0,iobuf,localsize,nslvar,1,0,fhandle);
	
	if(xy_rank == 0) 
	  fclose(fhandle);
       }	
      } else {
	if(xy_rank == 0) { 
	  sprintf(filename,"%s_%04d.dat","xy_slice",ixpos[nsl]);
	  fhandle=fopen(filename,"a");
	}
	
	slice_write(Grid,0,iobuf,localsize,nslvar,1,0,fhandle);
	
	if(xy_rank == 0){ 
	  fclose(fhandle);
	  
	  sprintf(filename,"%s_%04d.log","xy_slice",ixpos[nsl]);
	  
      std::fstream fptr;
	  int newfile = 0;
	  fptr.open(filename,std::ios::in);
	  if (!fptr) newfile = 1;
	  fptr.close();
	  
	  fptr.open(filename,std::ios::out|std::ios::app);
	  fptr.precision(10);
	  if (newfile) {      
	    fptr <<  nslvar << ' ' <<  Grid.gsize[0] << ' ' 
		 << Grid.gsize[1] << endl;
	    fptr << Physics.xy_var[0]  << ' ' 
		 << Physics.xy_var[1]  << ' ' 
		 << Physics.xy_var[2]  << ' ' 
		 << Physics.xy_var[3]  << ' ' 
		 << Physics.xy_var[4]  << ' ' 
		 << Physics.xy_var[5]  << ' ' 
		 << Physics.xy_var[6]  << ' ' 
		 << Physics.xy_var[7]  << ' ' 
		 << Physics.xy_var[8]  << ' ' 
		 << Physics.xy_var[9]  << ' ' 
		 << Physics.xy_var[10] << ' ' 
		 << Physics.xy_var[11] << endl;
	  }
	  fptr << Run.globiter << ' ' << Run.time << endl;
	  fptr.close();
	}
      }
    }     
  }
  if ((Run.enable_adios2) && (Run.globiter > Run.begin_iter)  && (Run.globiter > 0 )){
    adios2_end_step(Run.engineH);
  }

  free(iobuf);
  if (Run.enable_adios2_double)
    free(iobuf_double);
}

