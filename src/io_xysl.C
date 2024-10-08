#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include "io.H"
#include "precision.h"
#include "grid.H"
#include "run.H"
#include <cmath>
#include "comm_split.H"
#include "rt/rt.h"
#include "eos.H"
#include "limit_va.H"
#include <iostream>
#include "adios2_c.h"

using std::cout;
using std::endl;


int mpi_io_in     = 1;
int mpi_io_out    = 1;
int blocksize     = 8; // has only effect for mpi_io = 0 !

int nblocks,blsz;

static  adios2_variable *varList[20];
static  adios2_variable *diag_varList[20];
static  adios2_variable *eos_varList[20];
//adios2_adios* adiosH;
//adios2_io* ioH;
//adios2_engine *engineH;
int io_rank;

MPI_Datatype io_subarray,io_subarray_write;
MPI_Info io_info,io_info_write;
MPI_Comm io_xy_comm,io_z_comm,io_comm;

extern void WriteBackupFile(const char*,const int,const double);
extern void ReadBackupFile(const char*,const int,int*,double*);

void z_gather_io(const GridData&,const int,float*,int,float*,int);
void z_scatter_io(const GridData&,const int,float*,int,float*,int);
void xy_slice_write(const GridData&,const int,float*,int,FILE*);
void xy_slice_read(const GridData&,const int,float*,int,FILE*);

inline int imin(int i, int k){ return i < k ? i : k; }

void WriteHeaderFile(const GridData& Grid,const RunData& Run) {
  FILE* fh;
  char file[128];
  sprintf(file,"%s%s.%06d",Run.path_3D,"Header",Run.globiter);
  fh=fopen(file,"w");
  fprintf(fh,"%d  %d  %d  %f  %f  %f  %f  %f  %f",
          Grid.gsize[0],Grid.gsize[1],Grid.gsize[2],
	  Grid.dx[0],Grid.dx[1],Grid.dx[2],Run.time,Run.dt,1.0/sqrt(inv_va2max));
  fclose(fh);
}

void IO_Init(const GridData& Grid,RunData& Run) {
  int i;
  int gsz[3]; int lsz[3]; int str[3];

  MPI_Comm_dup(MPI_COMM_WORLD,&io_comm);
  MPI_Comm_dup(XY_COMM,&io_xy_comm);
  MPI_Comm_dup(ZCOL_COMM,&io_z_comm);
  MPI_Comm_rank(io_comm, &io_rank);

 for (i=0;i<=2;i++){
    gsz[i]=Grid.gsize[i];
    lsz[i]=Grid.lsize[i];
    str[i]=Grid.beg[i]-Grid.gbeg[i];
  }
  MPI_Type_create_subarray(3,gsz,lsz,str,MPI_ORDER_FORTRAN,
               MPI_FLOAT,&io_subarray_write);
  MPI_Type_commit(&io_subarray_write); 
  MPI_Info_create(&io_info_write);
  gsz[2] = Grid.gsize[2];
  lsz[2] = Grid.gsize[2];
  str[2] = 0;
  
  MPI_Type_create_subarray(3,gsz,lsz,str,MPI_ORDER_FORTRAN,
               MPI_FLOAT,&io_subarray);
  MPI_Type_commit(&io_subarray); 

  MPI_Info_create(&io_info);

  // make sure MPI IO errors leed to program termination
  MPI_File_set_errhandler(MPI_FILE_NULL,MPI_ERRORS_ARE_FATAL); 

  if (Run.enable_adios2){
     char filename[128];
     Run.adiosH = adios2_init_mpi(io_comm);
     Run.ioH = adios2_declare_io(Run.adiosH, "io_xyzl");
     adios2_set_engine(Run.ioH, "Sst");
     //adios2_set_engine(Run.ioH, "BPFile");
     //adios2_set_parameter(Run.ioH, "ProfileUnits", "Microseconds");
     //adios2_set_parameter(Run.ioH, "NumAggregators", "4");
     adios2_set_parameter(Run.ioH, "Threads", "1");
     adios2_set_parameter(Run.ioH, "QueueLimit", "2");
     adios2_set_parameter(Run.ioH, "InitialBufferSize", "1Gb");
     //adios2_set_parameter(Run.ioH, "SstVerbose", "5");
     //adios2_set_parameter(Run.ioH,"MarshalMethod","FFS");
     //adios2_set_parameter(Run.ioH,"DataTransport","RDMA");
     adios2_set_parameter(Run.ioH,"DataTransport","WAN");
     //adios2_set_parameter(Run.ioH,"ControlTransport","enet");
     //adios2_set_parameter(Run.ioH,"NetworkInterface", "hsn0");
     
     sprintf(filename,"%s_%s","result","prim");
     Run.engineH = adios2_open(Run.ioH, filename, adios2_mode_write);
     sprintf(filename,"%s_%s","diag","iter");
     //engineH2 = adios2_open(Run.ioH, filename, adios2_mode_write);
  }

  if ( Grid.gsize[2]%blocksize == 0 ){
    nblocks = Grid.gsize[2]/blocksize;
    blsz    = blocksize;
  } else {
    nblocks = Grid.gsize[2];
    blsz    = 1;
  }

  if(xy_rank+zcol_rank == 0)
      cout << "xy_slice_io: " 
       << mpi_io_in << ' ' << mpi_io_out << ' ' 
       << nblocks << ' ' << blsz 
       << endl;    
}

void IO_Finalize(const RunData& Run) {
  if(Run.enable_adios2) {
    adios2_close(Run.engineH);
    adios2_finalize(Run.adiosH);
  }
  MPI_Type_free(&io_subarray);
  MPI_Type_free(&io_subarray_write);
  MPI_Info_free(&io_info);
  MPI_Info_free(&io_info_write);
  MPI_Comm_free(&io_comm);
  MPI_Comm_free(&io_xy_comm);
  MPI_Comm_free(&io_z_comm);
}
////////////////////// Output /////////////////////////////////
void OutputSolution(const RunData& Run,const GridData& Grid,const PhysicsData& Physics)
{
  char filename[128];

  register int i,j,k,loc,v1,v2,node,v1_max,v2_max,var;
  double* U0 = (double*)Grid.U0;

  float* iobuf_loc;
  double* iobuf_loc_double;
  float* iobuf_glo=NULL;

  int sizex=Grid.lend[0]-Grid.lbeg[0]+1;
  int sizey=Grid.lend[1]-Grid.lbeg[1]+1;
  int sizez=Grid.lend[2]-Grid.lbeg[2]+1;
  
  int lsize=sizex*sizey*sizez; 
  int gsize=sizex*sizey*Grid.gsize[2];

  FILE* fh=NULL;
  MPI_File mfh;

  if(Physics.params[i_param_spitzer] > 0.0){
    if(Grid.procs[2]<3){
      v1_max = 9;
      v2_max = 1;
    }
    else if(Grid.procs[2]<9){
      v1_max = 3;
      v2_max = 3;
    } else {
      v1_max = 1;
      v2_max = 9;
    }
  } else {
    if(Grid.procs[2]<2){
      v1_max = 8;
      v2_max = 1;
    }
    else if(Grid.procs[2]<4){
      v1_max = 4;
      v2_max = 2;
    }
    else if(Grid.procs[2]<8){
      v1_max = 2;
      v2_max = 4;
    } else {
      v1_max = 1;
      v2_max = 8;
    }
  }

  LOCAL_LOOP(Grid,i,j,k) {
    node = Grid.node(i,j,k);
    Grid.U0[node] = Grid.U[node];
    Grid.U0[node].M /= Grid.U0[node].d;
    Grid.U0[node].e -= 0.5*Grid.U0[node].d*Grid.U0[node].M.sqr();
  }

  iobuf_loc = (float*)malloc(lsize*sizeof(float));
  if (Run.enable_adios2_double)
    iobuf_loc_double = (double*)malloc(lsize*sizeof(double));
  if (Run.enable_adios2==0)
    for(v2=0;v2<v2_max;v2++)
      if (zcol_rank == v2) iobuf_glo = (float*)malloc(gsize*sizeof(float));

  char adios2_var_name[128];
  int var_count;
  double clk;

  if(Run.enable_adios2) {
    int curr=0;
    if (Run.globiter <= imin(Run.resfreq,Run.slicefreq)||(Run.globiter > Run.begin_iter && varList[0]==NULL)){
    var_count = 0;
    for(v1=0;v1<v1_max;v1++) {
    for(v2=0;v2<v2_max;v2++) {
      var = v1*v2_max + v2;
      //sprintf(filename,"%s%s_%s.%06d",Run.path_3D,Run.resfile,"prim",Run.globiter);
      sprintf(adios2_var_name, "%s_%d.%06d", "result_prim",var,curr);
      uint64_t gdim[3],start[3],count[3];
      for(int i=0; i<3; i++) {
        gdim[i] = Grid.gsize[i];
        start[i] = Grid.beg[i]-Grid.gbeg[i];
        count[i] = Grid.lsize[i];
        if (xy_rank == 0)
            std::cout<<"define prim_var_name "<<adios2_var_name<<" "<<gdim[i]<<" "<<start[i]<<" "<<count[i]<<std::endl;
      }
      if(Run.enable_adios2_double) 
          varList[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_double,3, gdim,
                               start, count, adios2_constant_dims_true);
      else
          varList[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_float,3, gdim,
                               start, count, adios2_constant_dims_true);
      var_count++;
    }
    }
    }
  }

  int adios2_ind=0;
  adios2_step_status err; 
  float ** iobuf_loc_list; 
  iobuf_loc_list=(float**)malloc(var_count*sizeof(float*));
  for (int i=0;i<var_count;i++){
     iobuf_loc_list[i]=(float*)malloc(lsize*sizeof(float));
  }
  if (Run.enable_adios2 && Run.NeedsOutput() && Run.NeedsBackup()){
	  //cout<<"begin prim step"<<endl;
    adios2_begin_step(Run.engineH, adios2_step_mode_append, 10.0, &err);
  }
  for(v1=0;v1<v1_max;v1++) {
    for(v2=0;v2<v2_max;v2++) {
      var = v1*v2_max + v2;
      for(k=0;k<sizez;k++){
	for(j=0;j<sizey;j++){
	  for(i=0;i<sizex;i++){
            if (var < 8){
	      loc=Grid.node(i+Grid.ghosts[0],j+Grid.ghosts[1],k+Grid.ghosts[2]);
	      //iobuf_loc[i+j*sizex+k*sizex*sizey] = (float)U0[loc*Physics.NVAR+var];
              iobuf_loc_list[var][i+j*sizex+k*sizex*sizey]=(float)U0[loc*Physics.NVAR+var];
              if (Run.enable_adios2_double)
                iobuf_loc_double[i+j*sizex+k*sizex*sizey] = U0[loc*Physics.NVAR+var];
            } else {
              loc=Grid.node(i+Grid.ghosts[0],j+Grid.ghosts[1],k+Grid.ghosts[2]);
              //iobuf_loc[i+j*sizex+k*sizex*sizey] = (float)Grid.sflx0[loc];
              iobuf_loc_list[var][i+j*sizex+k*sizex*sizey] = (float)Grid.sflx0[loc];
              if (Run.enable_adios2_double)
                iobuf_loc_double[i+j*sizex+k*sizex*sizey] = Grid.sflx0[loc];
            }
	  }
	}
      }
      //adios2_put
      if (Run.enable_adios2 && Run.NeedsOutput() && Run.NeedsBackup()){
        //sprintf(adios2_var_name, "%s_%d.%06d", "prim",var,Run.globiter);
        clk = MPI_Wtime();
        if (Run.enable_adios2_double)
            adios2_put(Run.engineH,varList[adios2_ind],&iobuf_loc_double[0],adios2_mode_deferred);
        else
            adios2_put(Run.engineH,varList[adios2_ind],&iobuf_loc_list[adios2_ind][0],adios2_mode_deferred);
        size_t name_size;
        adios2_variable_name(adios2_var_name,&name_size,varList[adios2_ind]);
        if (io_rank==0)
          std::cout<<adios2_var_name<<" put prim rank "<<io_rank<<" io time "<<MPI_Wtime()-clk<<std::endl;
        adios2_ind++;
      }
      //else if (Run.NeedsBackup() && !Run.NeedsOutput() || Run.enable_adios2==0 || (Run.maxiter-Run.globiter)<Run.resfreq)
      if(!Run.enable_adios2){
	//sprintf(filename,"%s%s_%s_%d.%06d",Run.path_3D,Run.resfile,"prim",var,Run.globiter);
	//if(xy_rank == 0) cout << "write " << filename << endl;
        //MPI_File_open(io_comm,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,io_info_write,&mfh);
        //MPI_File_set_view(mfh,0,MPI_FLOAT,io_subarray_write,(char *)"native",io_info_write);
        //MPI_File_write_all(mfh,iobuf_loc,lsize,MPI_FLOAT,MPI_STATUS_IGNORE);
        //MPI_File_close(&mfh);      
        z_gather_io(Grid,v2,iobuf_loc_list[var],lsize,iobuf_glo,gsize);
      }
    }
    
    for(v2=0;v2<v2_max;v2++) {
      var = v1*v2_max + v2;       
      if (zcol_rank == v2){
	sprintf(filename,"%s%s_%s_%d.%06d",Run.path_3D,Run.resfile,"prim",var,Run.globiter);
	if(xy_rank == 0) cout << "write " << filename << endl;
	if(mpi_io_out == 0){
	  if(xy_rank == 0){
	    fh=fopen(filename,"w");
	    if(fh == NULL){
	      cout << "Error opening " << filename << ". Aborting ... " << endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	    }
	  }
	  for (k=0;k<nblocks;k++) xy_slice_write(Grid,0,&(iobuf_glo[k*blsz*sizex*sizey]),sizex*sizey,fh);
	  if(xy_rank == 0) fclose(fh);
	}else {
	  //if (Run.NeedsBackup() && !Run.NeedsOutput() || Run.enable_adios2==0 || (Run.maxiter-Run.globiter)<Run.resfreq)
	  if ( Run.enable_adios2==0 ){
	    MPI_File_open(io_xy_comm,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,io_info,&mfh);
	    MPI_File_set_view(mfh,0,MPI_FLOAT,io_subarray,(char *)"native",io_info);
	    MPI_File_write_all(mfh,iobuf_glo,gsize,MPI_FLOAT,MPI_STATUS_IGNORE);
	    MPI_File_close(&mfh);      
	  }
	}
      }
    }
  }
  if (Run.enable_adios2 && (Run.globiter > Run.begin_iter)  && Run.NeedsOutput() && Run.NeedsBackup()){
    adios2_end_step(Run.engineH);
  }


  free(iobuf_loc);
  free(iobuf_loc_list);
  if (Run.enable_adios2_double)
    free(iobuf_loc_double);
  if ( Run.enable_adios2==0 ){
    for(v2=0;v2<v2_max;v2++)
      if (zcol_rank == v2) free(iobuf_glo); 
  }
}

////////////////////// Backup-Restore /////////////////////////
//  HasBackupFile -- Return 1 if backup file exists, 0 otherwise.

int HasBackupFile(const char *file) {
  FILE *fptr;
  int res=0;

  fptr = fopen(file,"r");
  if( fptr ) {
    res = 1;
    fclose(fptr);
  }
  return res;
}

//-----------------------------------------------------------------
void BackupSolution(const RunData& Run,const GridData& Grid,
            const PhysicsData& Physics) {

  int v;
  int erase_flag = 0;
  char resfile[128];
  int prevbackup;
  double time;

  int max_vars = Physics.NVAR;
  if(Physics.params[i_param_spitzer] > 0.0)
    max_vars = Physics.NVAR+1;

  if( Run.rank==0 ){
    ReadBackupFile(Run.backfile,0,&prevbackup,&time);
    int flag = 1;
    if (Run.outcad > 0){
      double nnn = time/Run.outcad;
      int ceil = (int) nnn;
      if (nnn == (double) ceil)
        flag= 0;
      else
        flag= (int) (ceil+1)*Run.outcad-time;
    }
    cout <<"Backup "<< time << " " << Run.outcad << " " << time/Run.outcad << " " << flag << endl;
    if((prevbackup%Run.resfreq)&&(flag))
      erase_flag = 1;
  }

  OutputSolution(Run,Grid,Physics);
  
  if( Run.rank==0 ) {
    WriteBackupFile(Run.backfile,Run.globiter,Run.time);  
    WriteHeaderFile(Grid,Run);
  }
  
  // erase old backup file after successful output of new file, use here a Barrier
  // to ensure all processes are done writing before erasing previous
  MPI_Barrier(MPI_COMM_WORLD);
  if( Run.rank==0 ) {
    if(erase_flag){
      cout << "erase: " << Run.globiter << " " << prevbackup << endl;
      for(v=0;v<max_vars;v++) {
	sprintf(resfile,"%s%s_%s_%d.%06d",Run.path_3D,Run.resfile,
		"prim",v,prevbackup);
	remove(resfile);
      }
      sprintf(resfile,"%s%s.%06d",Run.path_3D,"Header",prevbackup);
      remove(resfile);
    }
  }
}

//-----------------------------------------------------------------
void RestoreSolution(RunData& Run,GridData& Grid,const PhysicsData& Physics){

  register int i,j,k,node,loc;

  char filename[128];

  double* U = (double*) Grid.U;
  int v1,v2,var,v1_max,v2_max;

  float* iobuf_loc;
  float* iobuf_glo=NULL;

  int sizex=(Grid.end[0]-Grid.beg[0]+1);
  int sizey=(Grid.end[1]-Grid.beg[1]+1);
  int sizez=(Grid.end[2]-Grid.beg[2]+1);

  int lsize=sizex*sizey*sizez; 
  int gsize=sizex*sizey*Grid.gsize[2];

  FILE* fh=NULL;
  MPI_File mfh;

   if(Physics.params[i_param_spitzer] > 0.0){
    if(Grid.procs[2]<3){
      v1_max = 9;
      v2_max = 1;
    }
    else if(Grid.procs[2]<9){
      v1_max = 3;
      v2_max = 3;
    } else {
      v1_max = 1;
      v2_max = 9;
    }
  } else {
    if(Grid.procs[2]<2){
      v1_max = 8;
      v2_max = 1;
    }
    else if(Grid.procs[2]<4){
      v1_max = 4;
      v2_max = 2;
    }
    else if(Grid.procs[2]<8){
      v1_max = 2;
      v2_max = 4;
    } else {
      v1_max = 1;
      v2_max = 8;
    }
  }

  // call on all cpus -> initialize p_bc, s_bc
  ReadBackupFile(Run.backfile,1,&Run.globiter,&Run.time);
  Run.begin_iter = Run.globiter;
  /*****************************************************************/
  if(Run.rank == 0){
    sprintf(filename,"%s%s_%s_%d.%06d",Run.path_3D,Run.resfile,
    "prim",0,Run.globiter);
    fh=fopen(filename,"r");
    if(fh == NULL){
      cout << "Error opening " << filename 
           << ". Try conservative ... " << endl;
    } else {
      fclose(fh);
    }
  }

  /*****************************************************************/

  iobuf_loc = (float*)malloc(lsize*sizeof(float));
  for(v2=0;v2<v2_max;v2++)
    if (zcol_rank == v2) iobuf_glo = (float*)malloc(gsize*sizeof(float));
  
  for(v1=0;v1<v1_max;v1++) {
    for(v2=0;v2<v2_max;v2++) {
      var = v1*v2_max + v2;
      if (zcol_rank == v2){
	sprintf(filename,"%s%s_%s_%d.%06d",Run.path_3D,Run.resfile,
		"prim",var,Run.globiter);
	if(xy_rank == 0) cout << "restore " << filename << endl;
	if(mpi_io_in == 0) {
	  if(xy_rank == 0){
	    fh=fopen(filename,"r");
	    if(fh == NULL){
	      cout << "Error opening " << filename << ". Aborting ... " << endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	    }
	  }
	  for (k=0;k<nblocks;k++)
	    xy_slice_read(Grid,0,&(iobuf_glo[k*blsz*sizex*sizey]),sizex*sizey,
			  fh);
	  if(xy_rank == 0) fclose(fh);
	} else {
	 cout <<"Run.rank="<<Run.rank<< "restore " << filename << endl;
	  MPI_File_open(io_xy_comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&mfh);
	  MPI_File_set_view(mfh,0,MPI_FLOAT,io_subarray,(char *) "native",
			    io_info);
	  MPI_File_read_all(mfh,iobuf_glo,gsize,MPI_FLOAT,
			    MPI_STATUS_IGNORE);
	  MPI_File_close(&mfh);
	}
      }
    }
    
    for(v2=0;v2<v2_max;v2++) {
      var = v1*v2_max + v2;
      z_scatter_io(Grid,v2,iobuf_glo,gsize,iobuf_loc,lsize);       
      for(k=0;k<sizez;k++){
	for(j=0;j<sizey;j++){
	  for(i=0;i<sizex;i++){
            if (var < 8){ 
	      loc=Grid.node(i+Grid.ghosts[0],j+Grid.ghosts[1],k+Grid.ghosts[2]);
	      U[loc*Physics.NVAR+var]=(double)iobuf_loc[i+j*sizex+k*sizex*sizey] ;
            } else {
              loc=Grid.node(i+Grid.ghosts[0],j+Grid.ghosts[1],k+Grid.ghosts[2]);
              Grid.sflx[loc]=(double)iobuf_loc[i+j*sizex+k*sizex*sizey] ;
            }
	  }
	}
      }
    }
  }
  
  LOCAL_LOOP(Grid,i,j,k) {
    node = Grid.node(i,j,k);
    Grid.U[node].e += 0.5*Grid.U[node].M.sqr()*Grid.U[node].d;
    Grid.U[node].M *= Grid.U[node].d;
    Grid.U0[node] = Grid.U[node];
  }

  free(iobuf_loc);
  for(v2=0;v2<v2_max;v2++)
    if (zcol_rank == v2) free(iobuf_glo); 
}
//////////////////// tvar_output ///////////////////////
void diag_output(const RunData& Run, const GridData& Grid,const PhysicsData& Physics, RTS *rts) {

  char filename[128];

  register int i,j,k,loc;

  float* iobuf_loc;
  double* iobuf_loc_double;
  float* iobuf_glo=NULL;
  int enable_adios2=0;

  int sizex=Grid.lend[0]-Grid.lbeg[0]+1;
  int sizey=Grid.lend[1]-Grid.lbeg[1]+1;
  int sizez=Grid.lend[2]-Grid.lbeg[2]+1;
  
  int lsize=sizex*sizey*sizez; 
  int gsize=sizex*sizey*Grid.gsize[2];
    
  FILE* fh=NULL;
  MPI_File mfh;

  int v1_max,v2_max,v1,v2,var;

  int max_vars=11;
  
  double* tvars[max_vars];
  int var_index[max_vars];
  char diag_names[max_vars][128];

  int tot_vars = 0;
  // tvars are only set properly aftr first iteration
  if(Run.iteration > 0){
    for(var=0;var<max_vars-1;var++){
      if(Run.diag_output[var] == 1){
	var_index[tot_vars]=var;
	tot_vars +=1;
      }
    }
  
    if(Physics.params[i_param_ambipolar] > 0.0){
      var_index[tot_vars] = max_vars-1;
      tot_vars +=1;
    }
  }
    
  
  if(Grid.procs[2]<tot_vars){
    v1_max = tot_vars;
    v2_max = 1;
  } else {
    v1_max = 1;
    v2_max = tot_vars;
  }
  
  tvars[0] = Grid.tvar1;
  tvars[1] = Grid.tvar2;
  tvars[2] = Grid.tvar3;
  tvars[3] = Grid.tvar4;
  tvars[4] = Grid.tvar5;
  tvars[5] = Grid.tvar6;
  tvars[6] = Grid.tvar7;
  tvars[7] = Grid.tvar8;
  tvars[8] = Grid.Qres;
  tvars[9] = Grid.Qvis;
  tvars[10] = Grid.Qamb;

  sprintf(diag_names[0],"%s","tvar1");
  sprintf(diag_names[1],"%s","tvar2");
  sprintf(diag_names[2],"%s","tvar3");
  sprintf(diag_names[3],"%s","tvar4");
  sprintf(diag_names[4],"%s","tvar5");
  sprintf(diag_names[5],"%s","tvar6");
  sprintf(diag_names[6],"%s","tvar7");
  sprintf(diag_names[7],"%s","tvar8");
  sprintf(diag_names[8],"%s","Qres");
  sprintf(diag_names[9],"%s","Qvis");
  sprintf(diag_names[10],"%s","Qamb");
  
  
  iobuf_loc = (float*)malloc(lsize*sizeof(float));
  if(Run.enable_adios2_double) 
    iobuf_loc_double = (double*)malloc(lsize*sizeof(double));
  for(v2=0;v2<v2_max;v2++)
    if (zcol_rank == v2) iobuf_glo = (float*)malloc(gsize*sizeof(float));

        //std::cout<<"Run globiter "<<Run.globiter<<std::endl;
  char adios2_var_name[128];
  if(Run.enable_adios2) {
    int curr=0;
    int var_count;
    if ((Run.globiter <= Run.resfreq) ||(Run.globiter > Run.begin_iter && diag_varList[0]==NULL)){
    var_count = 0;
    for(v1=0;v1<v1_max;v1++) {
    for(v2=0;v2<v2_max;v2++) {
      var = var_index[v1*v2_max + v2];
      sprintf(adios2_var_name, "%s.%06d", diag_names[var],curr);
      uint64_t gdim[3],start[3],count[3];
      for(int i=0; i<3; i++) {
        gdim[i] = Grid.gsize[i];
        start[i] = Grid.beg[i]-Grid.gbeg[i];
        count[i] = Grid.lsize[i];
        if (xy_rank == 0)
            std::cout<<"define "<<adios2_var_name<<" "<<gdim[i]<<" "<<start[i]<<" "<<count[i]<<std::endl;
      }
      if(Run.enable_adios2_double) 
          diag_varList[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_double,3, gdim,
                               start, count, adios2_constant_dims_true);
      else
          diag_varList[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_float,3, gdim,
                               start, count, adios2_constant_dims_true);
      var_count++;
    }
    }
    }
  }

  int adios2_ind=0;
  if (Run.enable_adios2 && (Run.globiter > Run.begin_iter)  && (Run.globiter > 0)){
    adios2_step_status serr; 
    adios2_begin_step(Run.engineH, adios2_step_mode_append, 10.0, &serr);
  }

  for(v1=0;v1<v1_max;v1++){
    for(v2=0;v2<v2_max;v2++){
      var = var_index[v1*v2_max + v2];
      for(k=0;k<sizez;k++){
	for(j=0;j<sizey;j++){
	  for(i=0;i<sizex;i++){
	    loc=Grid.node(i+Grid.ghosts[0],j+Grid.ghosts[1],k+Grid.ghosts[2]);
	    iobuf_loc[i+j*sizex+k*sizex*sizey] = (float) tvars[var][loc];
            if(Run.enable_adios2_double) 
              iobuf_loc_double[i+j*sizex+k*sizex*sizey] = tvars[var][loc];
	  }
	}
      }
      //adios2_put
      if (Run.enable_adios2){
        size_t step;
        adios2_error err; 
        //sprintf(adios2_var_name, "%s_%d.%06d", diag_names[var],var,Run.globiter);
        if (Run.enable_adios2_double)
            err=adios2_put(Run.engineH,diag_varList[adios2_ind],&iobuf_loc_double[0],adios2_mode_sync);
        else
            err=adios2_put(Run.engineH,diag_varList[adios2_ind],&iobuf_loc[0],adios2_mode_sync);
        adios2_ind++;
        //std::cout<<adios2_var_name<<"put diag rank "<<io_rank<<" lsize "<<lsize<<std::endl;
      }
      if ( Run.enable_adios2==0 ){
	  //sprintf(filename,"%s%s.%06d",Run.path_3D,diag_names[var],Run.globiter);
	  //cout <<Run.rank<< " write " << filename << endl;
          //MPI_File_open(io_comm,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
          //              io_info_write,&mfh);
          //MPI_File_set_view(mfh,0,MPI_FLOAT,io_subarray_write,(char *) "native",
          //                  io_info_write);
          //MPI_File_write_all(mfh,iobuf_loc,lsize,MPI_FLOAT,MPI_STATUS_IGNORE);
          //MPI_File_close(&mfh);
      z_gather_io(Grid,v2,iobuf_loc,lsize,iobuf_glo,gsize);
      }
    }
    
    for(v2=0;v2<v2_max;v2++){
      var = var_index[v1*v2_max + v2];
      if (zcol_rank == v2){
	sprintf(filename,"%s%s.%06d",Run.path_3D,diag_names[var],Run.globiter);
	if(xy_rank == 0) cout << "write " << filename << endl;
	if(mpi_io_out == 0) { 
	  if(xy_rank == 0){
	    fh=fopen(filename,"w");
	    if(fh == NULL){
	      cout << "Error opening " << filename << ". Aborting ... " << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	    }
	  }
	  for (k=0;k<nblocks;k++)
	    xy_slice_write(Grid,0,&(iobuf_glo[k*blsz*sizex*sizey]),sizex*sizey,
			   fh);
	  if(xy_rank == 0) fclose(fh);
	} else{
	  //  if (!Run.globiter || Run.enable_adios2==0 || (Run.maxiter-Run.globiter)<Run.resfreq) 
	  if ( Run.enable_adios2==0 ){
	  MPI_File_open(io_xy_comm,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
			io_info,&mfh);
	  MPI_File_set_view(mfh,0,MPI_FLOAT,io_subarray,(char *) "native",
			    io_info);
	  MPI_File_write_all(mfh,iobuf_glo,gsize,MPI_FLOAT,MPI_STATUS_IGNORE);
	  MPI_File_close(&mfh);
  	  }
	}
      }
    }
  }
  
  if (Run.enable_adios2 && (Run.globiter > Run.begin_iter)  && (Run.globiter > 0)){
    adios2_end_step(Run.engineH);
  }

  free(iobuf_loc); 
  if(Run.enable_adios2_double)
    free(iobuf_loc_double);
  for(v2=0;v2<v2_max;v2++)
    if (zcol_rank == v2) free(iobuf_glo); 
  
}

//////////////////// eos_output ///////////////////////
void eos_output(const RunData& Run, const GridData& Grid,const PhysicsData& Physics, RTS *rts) {

  char filename[128];

  register int i,j,k,loc;

  float* iobuf_loc;
  double* iobuf_loc_double;
  float* iobuf_glo=NULL;

  int sizex=Grid.lend[0]-Grid.lbeg[0]+1;
  int sizey=Grid.lend[1]-Grid.lbeg[1]+1;
  int sizez=Grid.lend[2]-Grid.lbeg[2]+1;
 
  int enable_adios2=0; 
  int lsize=sizex*sizey*sizez; 
  int gsize=sizex*sizey*Grid.gsize[2];
    
  FILE* fh=NULL;
  MPI_File mfh;

  int v1_max,v2_max,v1,v2,var;

  int max_vars = 14;
  
  char eos_names[max_vars][128];
  double* eos_vars[max_vars];

  //Only write out variables that re allocated based on phyasics configuration
  int var_init[max_vars];
  for(var=0;var<max_vars;var++)
    var_init[var] = 0;
  
  var_init[0] = 1;
  var_init[1] = 1;
  var_init[2] = 1;
  var_init[3] = 1;
  var_init[4] = 1;
  var_init[5] = 1;
  var_init[6] = 1;
  var_init[7] = 1;
  var_init[8] = 1;

  if(Physics.rt_ext[i_ext_cor] >= 1)
    var_init[9] = 1;
  if(Physics.rt_ext[i_ext_cor] == 2){
    var_init[10] = 1;
    var_init[11] = 1;
    var_init[12] = 1;
    var_init[13] = 1;
  }
  
  int var_index[max_vars];

  int tot_vars = 0;
  for(var=0;var<max_vars;var++){
    if((Run.eos_output[var] == 1) && (var_init[var] == 1)){
      var_index[tot_vars]=var;
      tot_vars +=1;
    }
  }
      
  if(Grid.procs[2]<tot_vars){
    v1_max = tot_vars;
    v2_max = 1;
  } else {
    v1_max = 1;
    v2_max = tot_vars;
  }

  // This is ugly, but for now it works so I will stick with it -- DP 
  sprintf(eos_names[0],"%s","eosT");
  sprintf(eos_names[1],"%s","eosP");
  sprintf(eos_names[2],"%s","eosne");
  sprintf(eos_names[3],"%s","eosrhoi");
  sprintf(eos_names[4],"%s","eosamb");
  sprintf(eos_names[5],"%s","Qtot");
  sprintf(eos_names[6],"%s","tau");
  sprintf(eos_names[7],"%s","Jtot");
  sprintf(eos_names[8],"%s","Stot");
  sprintf(eos_names[9],"%s","QxCor");
  sprintf(eos_names[10],"%s","QxH");
  sprintf(eos_names[11],"%s","QxMg");
  sprintf(eos_names[12],"%s","QxCa");
  sprintf(eos_names[13],"%s","QxChr");

  eos_vars[0] = Grid.temp;
  eos_vars[1] = Grid.pres;
  eos_vars[2] = Grid.ne;
  eos_vars[3] = Grid.rhoi;
  eos_vars[4] = Grid.amb;
  eos_vars[5] = Grid.Qtot;
  eos_vars[6] = Grid.Tau;
  eos_vars[7] = Grid.Jtot;
  eos_vars[8] = Grid.Stot;
  eos_vars[9] = Grid.Qthin;
  eos_vars[10] = Grid.QH;
  eos_vars[11] = Grid.QMg;
  eos_vars[12] = Grid.QCa;
  eos_vars[13] = Grid.QChr;
 
  iobuf_loc = (float*)malloc(lsize*sizeof(float));
  if(Run.enable_adios2_double)
    iobuf_loc_double = (double*)malloc(lsize*sizeof(double));
  for(v2=0;v2<v2_max;v2++)
    if (zcol_rank == v2) iobuf_glo = (float*)malloc(gsize*sizeof(float));

  char adios2_var_name[128];
  int curr=0;
  int var_count;
  double clk;

  if(Run.enable_adios2) {
    if ((Run.globiter > Run.begin_iter) && eos_varList[0]==NULL || (Run.globiter < Run.resfreq)){
    var_count = 0;
    for(v1=0;v1<v1_max;v1++) {
    for(v2=0;v2<v2_max;v2++) {
      var = var_index[v1*v2_max + v2];
      sprintf(adios2_var_name, "%s.%06d", eos_names[var],curr);
      uint64_t gdim[3],start[3],count[3];
      for(int i=0; i<3; i++) {
        gdim[i] = Grid.gsize[i];
        start[i] = Grid.beg[i]-Grid.gbeg[i];
        count[i] = Grid.lsize[i];
        if (xy_rank == 0)
            std::cout<<"define "<< adios2_var_name<<" "<<gdim[i]<<" "<<start[i]<<" "<<count[i]<<std::endl;
      }
      if(Run.enable_adios2_double) 
          eos_varList[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_double,3, gdim,
                               start, count, adios2_constant_dims_true);
      else
          eos_varList[var_count] = adios2_define_variable(Run.ioH, adios2_var_name, adios2_type_float,3, gdim,
                               start, count, adios2_constant_dims_true);
      var_count++;
    }
    }
    }
  }

  int adios2_ind=0;
  if (Run.enable_adios2 && (Run.globiter > Run.begin_iter) && (Run.globiter > 0)){
    adios2_step_status serr; 
    //cout<<"begin eos"<<endl;
    adios2_begin_step(Run.engineH, adios2_step_mode_append, 10.0, &serr);
  }

  for(v1=0;v1<v1_max;v1++){
    for(v2=0;v2<v2_max;v2++){
      var = var_index[v1*v2_max + v2];
      for(k=0;k<sizez;k++){
	for(j=0;j<sizey;j++){
	  for(i=0;i<sizex;i++){
	    loc=Grid.node(i+Grid.ghosts[0],j+Grid.ghosts[1],k+Grid.ghosts[2]);
	    iobuf_loc[i+j*sizex+k*sizex*sizey] = (float) eos_vars[var][loc];
            if(Run.enable_adios2_double)
              iobuf_loc_double[i+j*sizex+k*sizex*sizey] = eos_vars[var][loc];
	  }
	}
      }
      //adios2_put
      if (Run.enable_adios2 && (Run.globiter > Run.begin_iter) && (Run.globiter > 0)){
        adios2_error err;
        size_t step;
        //sprintf(adios2_var_name, "%s.%06d", eos_names[var],Run.globiter);
        //cout<<"put eos data "<<endl;
        if (Run.enable_adios2_double)
            err=adios2_put(Run.engineH,eos_varList[adios2_ind],&iobuf_loc_double[0],adios2_mode_sync);
        else
            err=adios2_put(Run.engineH,eos_varList[adios2_ind],&iobuf_loc[0],adios2_mode_sync);
        adios2_ind++;
        //std::cout<<adios2_var_name<<"put eos rank "<<io_rank<<" lsize "<<lsize<<std::endl;
      }
      if ( Run.enable_adios2==0 ){
//	sprintf(filename,"%s%s.%06d",Run.path_3D,eos_names[var],Run.globiter);
//	if(xy_rank == 0) cout << "write " << filename << endl;
//      MPI_File_open(io_comm,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
//                    io_info_write,&mfh);
//      MPI_File_set_view(mfh,0,MPI_FLOAT,io_subarray_write,(char *) "native",
//                        io_info_write);
//      MPI_File_write_all(mfh,iobuf_loc,lsize,MPI_FLOAT,MPI_STATUS_IGNORE);
//      MPI_File_close(&mfh);
      z_gather_io(Grid,v2,iobuf_loc,lsize,iobuf_glo,gsize);
      }
    }
    
    for(v2=0;v2<v2_max;v2++){
      var = var_index[v1*v2_max + v2];
      if (zcol_rank == v2){
	sprintf(filename,"%s%s.%06d",Run.path_3D,eos_names[var],Run.globiter);
	if(xy_rank == 0) cout << "write " << filename << endl;
	if(mpi_io_out == 0) { 
	  if(xy_rank == 0){
	    fh=fopen(filename,"w");
	    if(fh == NULL){
	      cout << "Error opening " << filename << ". Aborting ... " << endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	    }
	  }
	  for (k=0;k<nblocks;k++)
	    xy_slice_write(Grid,0,&(iobuf_glo[k*blsz*sizex*sizey]),sizex*sizey,
			   fh);
	  if(xy_rank == 0) fclose(fh);
	} else{ 
	 //if (!Run.globiter || !Run.enable_adios2 || (Run.maxiter-Run.globiter)<Run.resfreq) 
	  if ( Run.enable_adios2==0 ){
	  MPI_File_open(io_xy_comm,filename,MPI_MODE_CREATE | MPI_MODE_WRONLY,
			io_info,&mfh);
	  MPI_File_set_view(mfh,0,MPI_FLOAT,io_subarray,(char *) "native",
			    io_info);
	  MPI_File_write_all(mfh,iobuf_glo,gsize,MPI_FLOAT,MPI_STATUS_IGNORE);
	  MPI_File_close(&mfh);
	  }
	}
      }
    }
  }
  
  if (Run.enable_adios2 && (Run.globiter > Run.begin_iter)  && (Run.globiter > 0)){
    adios2_end_step(Run.engineH);
  }

  free(iobuf_loc); 
  if (Run.enable_adios2_double)
    free(iobuf_loc_double);
  for(v2=0;v2<v2_max;v2++)
    if (zcol_rank == v2) free(iobuf_glo); 
  
}

//===========================================================================
void z_gather_io(const GridData& Grid, const int iroot,float* vloc,int nloc,
         float* vglo, int nglo){
  
  int np;
  const int nprocs = Grid.procs[2];
  
  static int* recvcounts;
  static int* offsets;
  static int recvbufsize;
  static int ini_flag = 1;
  if (ini_flag == 1){
    recvcounts  = (int*) malloc(nprocs*sizeof(int));
    offsets     = (int*) malloc(nprocs*sizeof(int)); 
 
    //std::cout<<"z_gather_io"<<std::endl;
    MPI_Allgather(&nloc,1,MPI_INT,recvcounts,1,MPI_INT,io_z_comm);
    
    offsets[0]=0;
    for (np=1;np<nprocs;np++){  
      offsets[np]=offsets[np-1]+recvcounts[np-1];
      //std::cout<<"offset "<<np<<"="<<offsets[np]<<std::endl;
      //std::cout<<"recvcounts "<<np<<"="<<recvcounts[np-1]<<std::endl;
    }
    recvbufsize=0;
    for (np=0;np<nprocs;np++){
      recvbufsize+=recvcounts[np];
      //std::cout<<nglo<<" recvcounts "<<np<<"="<<recvcounts[np]<<std::endl;
    }

    ini_flag = 0;
  }
    
  if(zcol_rank == iroot){
    if(nglo != recvbufsize){
      cout << "z_gather_io: nglo != recvbufsize " << zcol_rank << ' ' 
       << nglo << ' ' << recvbufsize << endl;
      MPI_Abort(io_comm,1);
    }
  }
  
  if(nloc != recvcounts[zcol_rank]){
    cout << "z_gather_io: nloc != recvcounts " << zcol_rank << ' ' 
         << nloc << ' ' << recvcounts[zcol_rank] << endl;
    MPI_Abort(io_comm,1);
  }
  
  MPI_Gatherv(vloc,nloc,MPI_FLOAT,vglo,recvcounts,offsets,
              MPI_FLOAT,iroot,io_z_comm);
}

//=======================================================================
void z_scatter_io(const GridData& Grid, const int iroot,float* vglo,int nglo,
          float* vloc, int nloc){
  
  int np;
  const int nprocs = Grid.procs[2];
  
  static int* sendcounts;
  static int* offsets;
  static int sendbufsize;
  static int ini_flag = 1;

  if (ini_flag == 1){
    sendcounts  = (int*) malloc(nprocs*sizeof(int));
    offsets     = (int*) malloc(nprocs*sizeof(int)); 
 
    MPI_Allgather(&nloc,1,MPI_INT,sendcounts,1,MPI_INT,io_z_comm);
    
    offsets[0]=0;
    for (np=1;np<nprocs;np++){  
      offsets[np]=offsets[np-1]+sendcounts[np-1];
    }
    sendbufsize=0;
    for (np=0;np<nprocs;np++){
      sendbufsize+=sendcounts[np];
    }

    ini_flag = 0;
  }
    
  if(zcol_rank == iroot){
    if(nglo != sendbufsize){
      cout << "z_scatter_io: nglo != sendbufsize " << zcol_rank << ' ' 
       << nglo << ' ' << sendbufsize << endl;
      MPI_Abort(io_comm,1);
    }
  }
  
  if(nloc != sendcounts[zcol_rank]){
    cout << "z_scatter_io: nloc != sendcounts " << zcol_rank << ' ' 
         << nloc << ' ' << sendcounts[zcol_rank] << endl;
    MPI_Abort(io_comm,1);
  }
  
  MPI_Scatterv(vglo,sendcounts,offsets,MPI_FLOAT,vloc,nloc,MPI_FLOAT,
           iroot,io_z_comm);
}
//=========================================================================
void xy_slice_write(const GridData& Grid,const int iroot,float* vloc,int nloc,
            FILE* fhandle){
  
  int np,i,k,iloc,iglo,v;
  int Gb[2], Ge[2], bounds[4];
  int isz,ksz,ioff,koff;

  const int nprocs = Grid.procs[0]*Grid.procs[1];
  const int xysize = Grid.gsize[0]*Grid.gsize[1];
  
  static int* proc_bounds;
  static int ini_flag = 1;

  int* recvcounts;
  int* offsets;
  float* recvbuf=NULL;
  float* iobuf=NULL;
  int recvbufsize;
 
  if(ini_flag == 1){
    proc_bounds = (int*) malloc(4*nprocs*sizeof(int));
    bounds[0]=Grid.beg[0];
    bounds[1]=Grid.end[0];
    bounds[2]=Grid.beg[1];
    bounds[3]=Grid.end[1];  
    MPI_Allgather(bounds,4,MPI_INT,proc_bounds,4,MPI_INT,io_xy_comm);
    ini_flag = 0;
  }

  recvcounts  = (int*) malloc(nprocs*sizeof(int));
  offsets     = (int*) malloc(nprocs*sizeof(int));

  for (np=0;np<nprocs;np++){
    Gb[0]=proc_bounds[4*np+0];
    Ge[0]=proc_bounds[4*np+1];
    Gb[1]=proc_bounds[4*np+2];
    Ge[1]=proc_bounds[4*np+3];  
    recvcounts[np]=(Ge[0]-Gb[0]+1)*(Ge[1]-Gb[1]+1)*blsz; 
  }
  offsets[0]=0;
  for (np=1;np<nprocs;np++){  
    offsets[np]=offsets[np-1]+recvcounts[np-1];
  }
  
  recvbufsize=0;
  for (np=0;np<nprocs;np++){   
    recvbufsize+=recvcounts[np];
  }
    
  if(xy_rank == iroot){
    if (Grid.gsize[0]*Grid.gsize[1]*blsz != recvbufsize){
      cout << "xy_gather_io: nglo != recvbufsize " << xy_rank << ' ' 
           << Grid.gsize[0]*Grid.gsize[1]*blsz << ' ' << recvbufsize << endl;
      MPI_Abort(io_comm,1);
    }
    
    recvbuf = (float*) malloc(recvbufsize*sizeof(float));
    iobuf   = (float*) malloc(xysize*sizeof(float));
  }
  
  if(nloc*blsz != recvcounts[xy_rank]){
    cout << "xy_gather_io: nloc != recvcounts " << xy_rank << ' ' 
         << nloc << ' ' << recvcounts[xy_rank] << endl;
    MPI_Abort(io_comm,1);
  }

  MPI_Gatherv(vloc,nloc*blsz,MPI_FLOAT,recvbuf,recvcounts,offsets,
              MPI_FLOAT,iroot,io_xy_comm);

  if (xy_rank == iroot){
    for (v=0;v<blsz;v++){
      for(np=0; np<nprocs; np++){
    Gb[0]=proc_bounds[4*np+0];
    Ge[0]=proc_bounds[4*np+1];
    Gb[1]=proc_bounds[4*np+2];
    Ge[1]=proc_bounds[4*np+3];

    isz  = Ge[0]-Gb[0]+1;
    ksz  = Ge[1]-Gb[1]+1;
    ioff = Gb[0]-Grid.gbeg[0];
    koff = Gb[1]-Grid.gbeg[1];
    
    for (k=0;k<ksz;k++){ 
      for (i=0;i<isz;i++){
        iloc = i + isz*(k+v*ksz) + offsets[np];
        iglo = (i+ioff) + (k+koff)*Grid.gsize[0];
        iobuf[iglo] = recvbuf[iloc];
      }
    }
      }
      int wcounts=fwrite(iobuf,sizeof(float),xysize,fhandle);
      if(wcounts != xysize){
    cout << "fwrite io error: " << wcounts << ' ' << xysize << endl;
    MPI_Abort(io_comm,1);
      }
    }
  }

  if (xy_rank == iroot){
    free(recvbuf); free(iobuf);
  }
  free(recvcounts); free(offsets);

}

//=============================================================================
void xy_slice_read(const GridData& Grid,const int iroot,float* vloc, 
           int nloc, FILE* fhandle){
  
  int np,i,k,iloc,iglo,v;
  int Gb[2], Ge[2], bounds[4];
  int isz,ksz,ioff,koff;

  const int nprocs = Grid.procs[0]*Grid.procs[1];
  const int xysize = Grid.gsize[0]*Grid.gsize[1];

  static int* proc_bounds;
  static int ini_flag = 1;

  int* sendcounts;
  int* offsets;
  float* sendbuf=NULL;
  float* iobuf=NULL;
  int sendbufsize;

  if (ini_flag == 1){
    proc_bounds = (int*) malloc(4*nprocs*sizeof(int));
    bounds[0]=Grid.beg[0];
    bounds[1]=Grid.end[0];
    bounds[2]=Grid.beg[1];
    bounds[3]=Grid.end[1];   
    MPI_Allgather(bounds,4,MPI_INT,proc_bounds,4,MPI_INT,io_xy_comm);
    ini_flag = 0;
  }

  sendcounts  = (int*) malloc(nprocs*sizeof(int));
  offsets     = (int*) malloc(nprocs*sizeof(int));

  for (np=0;np<nprocs;np++){
    Gb[0]=proc_bounds[4*np+0];
    Ge[0]=proc_bounds[4*np+1];
    Gb[1]=proc_bounds[4*np+2];
    Ge[1]=proc_bounds[4*np+3];  
    sendcounts[np]=(Ge[0]-Gb[0]+1)*(Ge[1]-Gb[1]+1)*blsz; 
  }
  offsets[0]=0;
  for (np=1;np<nprocs;np++){  
    offsets[np]=offsets[np-1]+sendcounts[np-1];
  }
  
  sendbufsize=0;
  for (np=0;np<nprocs;np++){   
    sendbufsize+=sendcounts[np];
  }

  if(xy_rank == iroot){
    if (xysize*blsz != sendbufsize){
      cout << "xy_scatter_io: nglo != sendbufsize " << xy_rank << ' ' 
           << xysize*blsz << ' ' << sendbufsize << endl;
      MPI_Abort(io_comm,1);
    }
    
    sendbuf = (float*) malloc(sendbufsize*sizeof(float));
    iobuf   = (float*) malloc(xysize*sizeof(float));
  }
  
  if(nloc*blsz != sendcounts[xy_rank]){
    cout << "xy_scatter_io: nloc != sendcounts " << xy_rank << ' ' 
         << nloc << ' ' << sendcounts[xy_rank] << endl;
    MPI_Abort(io_comm,1);
  }

  if (xy_rank == iroot){
    for (v=0;v<blsz;v++){
      int rcounts=fread(iobuf,sizeof(float),xysize,fhandle);
      if(rcounts != xysize){
    cout << "fread io error: " << rcounts << ' ' << xysize << endl;
    MPI_Abort(io_comm,1);
      }
      for(np=0; np<nprocs; np++){
    Gb[0]=proc_bounds[4*np+0];
    Ge[0]=proc_bounds[4*np+1];
    Gb[1]=proc_bounds[4*np+2];
    Ge[1]=proc_bounds[4*np+3];
    
    isz  = Ge[0]-Gb[0]+1;
    ksz  = Ge[1]-Gb[1]+1;
    ioff = Gb[0]-Grid.gbeg[0];
    koff = Gb[1]-Grid.gbeg[1];
    
    for (k=0;k<ksz;k++){ 
      for (i=0;i<isz;i++){
        iloc = i + isz*(k+v*ksz) + offsets[np];
        iglo = (i+ioff) + (k+koff)*Grid.gsize[0];
        sendbuf[iloc] = iobuf[iglo];
      }
    }
      } 
    }
  }

  MPI_Scatterv(sendbuf,sendcounts,offsets,MPI_FLOAT,vloc,nloc*blsz,MPI_FLOAT,
               iroot,io_xy_comm);

  if (xy_rank == iroot){
    free(sendbuf);
    free(iobuf);
  }

  free(sendcounts); free(offsets);

}
//==========================================================================
