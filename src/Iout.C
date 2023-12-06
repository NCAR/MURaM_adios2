#include <mpi.h>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include "physics.H"
#include "grid.H"
#include "run.H"
#include "eos.H"
#include "comm_split.H"
#include "rt/rt.h"
#include <iostream>

using namespace std;
static adios2_variable *varList1[1];
extern double get_Iout(int,int);

extern void slice_write(const GridData&,const int,float*,int,int,const int,
			const int,FILE*);
//======================================================================
void Iout(const RunData&  Run, const GridData& Grid, 
	  const PhysicsData& Physics,RTS *rts) {
  
  register int ind;

  int jbeg = Grid.lbeg[1], jend = Grid.lend[1];
  int kbeg = Grid.lbeg[2], kend = Grid.lend[2];

  int localsize  = Grid.lsize[1]*Grid.lsize[2];

  float* icloc;  
  double* icloc_double;  

  char filename[128];

  FILE* fhandle=NULL;
  int enable_adios2=1;
  //MPI_File fhandle_mpi;
  //int offset;
  int var_count, adios2_ind;
  var_count = 0;
  adios2_ind = 0;
  adios2_step_status err;
  if ((Run.enable_adios2) && (Run.globiter > Run.begin_iter)  && (Run.globiter>0)){
    adios2_begin_step(Run.engineH, adios2_step_mode_append, 10.0, &err);
    cout <<"Iout begin step"<<endl;
  }
  if(Grid.is_gend[0]){  
    //std::cout<<"inside Iout"<<std::endl;
    icloc = (float*) malloc(localsize*sizeof(float));
    if(Run.enable_adios2_double) 
      icloc_double = (double*) malloc(localsize*sizeof(double));
    rts->UpdateIout();
    for (int k=kbeg; k<=kend; k++){
      for (int j=jbeg; j<=jend; j++){
	    ind = j-jbeg + (k-kbeg)*Grid.lsize[1];       
	    icloc[ind] = (float) rts->Iout(j,k);
            if(Run.enable_adios2_double) 
	      icloc_double[ind] = (double) rts->Iout(j,k);
      }
    }
    char adios2_var_name[128];
    int curr=0;

    if(Physics.slice[i_sl_collect] == 0) {
        if(Run.enable_adios2){
          if ((Run.globiter > Run.begin_iter && varList1[0]==NULL) ||Run.globiter < std::min(Run.resfreq,Run.slicefreq)){
            sprintf(adios2_var_name,"%s.%06d","I_out",curr);
            uint64_t gdim[2],start[2],count[2];
            for(int i=0;i<2; i++){
              gdim[i] = Grid.gsize[i+1];
              start[i] = Grid.beg[i+1]-Grid.gbeg[i+1];
              count[i] = Grid.lsize[i+1];
            }
            //std::cout<<"I_out name "<<adios2_var_name<<" "<<localsize<<" "<<start[0]<<" "<<start[1]<<" "<<start[2]<<" "<<count[0]<<" "<<count[1]<<" "<<count[2]<<endl;

            if(Run.enable_adios2_double) {
              varList1[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_double,2, gdim,
                           start, count, adios2_constant_dims_true);
            }
            else{
              varList1[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_float,2, gdim,
                           start, count, adios2_constant_dims_true);
            }
            var_count++;
          }
        }
        if ((Run.enable_adios2)&& (Run.globiter > Run.begin_iter)  && (Run.globiter>0)){
          //std::cout<<Run.rank<<" before Iout put "<<xy_rank<<endl;
          if (Run.enable_adios2_double)
            adios2_put(Run.engineH,varList1[adios2_ind],&icloc_double[0],adios2_mode_sync);
          else
            adios2_put(Run.engineH,varList1[adios2_ind],&icloc[0],adios2_mode_sync);
          adios2_ind++;
          //std::cout<<Run.rank<<" after Iout put "<<yz_rank<<endl;
        } 
      if (!Run.enable_adios2 || !Run.globiter||(Run.maxiter-Run.globiter)<std::min(Run.resfreq,Run.slicefreq)){
	//cout<<"Iout rank "<<Run.rank<<" "<<yz_rank<<endl;
      if(yz_rank == 0) {  
        sprintf(filename,"%s%s.%06d",Run.path_2D,"I_out",Run.globiter);
        fhandle=fopen(filename,"w");
        
        float header[4];            
        header[0] = (float) 1;
        header[1] = (float) Grid.gsize[1];
        header[2] = (float) Grid.gsize[2];
        header[3] = (float) Run.time;
        fwrite(header,sizeof(float),4,fhandle);
      }
      
      slice_write(Grid,0,icloc,localsize,1,1,2,fhandle);
      
      if(yz_rank == 0)
	fclose(fhandle);
      }
    } else {
      if(yz_rank == 0){
	sprintf(filename,"I_out.dat");
	fhandle=fopen(filename,"a");
      }
      
      slice_write(Grid,0,icloc,localsize,1,1,2,fhandle);
      
      if(yz_rank == 0){
	fclose(fhandle);
	
	fstream fptr;
	int newfile = 0;
	fptr.open("I_out.log",ios::in);
	if (!fptr) newfile = 1;
	fptr.close();
	
	fptr.open("I_out.log",ios::out|ios::app);
	fptr.precision(10);
	if (newfile)       
	  fptr << '1' << ' ' <<  Grid.gsize[1] << ' ' 
	       << Grid.gsize[2] << endl;
	fptr << Run.globiter << ' ' << Run.time << endl;
	fptr.close();
      }
    }
    free(icloc);
    if (Run.enable_adios2_double)
      free(icloc_double);
  }
  if ((Run.enable_adios2) && (Run.globiter > Run.begin_iter)  && (Run.globiter>0))
    adios2_end_step(Run.engineH);

}
