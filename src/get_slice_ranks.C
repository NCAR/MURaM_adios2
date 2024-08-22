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

void get_slice_ranks(const RunData&  Run, const GridData& Grid, const PhysicsData& Physics) {

   FILE* fh;
   fh=fopen("slice_ranks.txt","w");
   fclose(fh);

  int xy_nranks=Grid.procs[0]*Grid.procs[1];
  int xy_ranks[xy_nranks];

  for (int nsl = 0; nsl<Physics.slice[i_sl_xy]; nsl++){

    if ( (Grid.beg[2] <= Physics.xy_lev[nsl]+Grid.gbeg[2] ) and 
         (Grid.end[2] >= Physics.xy_lev[nsl]+Grid.gbeg[2] )){
      
      MPI_Gather(&Run.rank,1,MPI_INT,xy_ranks,1,MPI_INT,0,XY_COMM);
      if (xy_rank == 0){
         fh=fopen("slice_ranks.txt","a");     
         if (nsl==0)
             fprintf(fh,"{\n");
         fprintf(fh,"\"%s_%04d\": {","xy_slice",Physics.xy_lev[nsl]);
         for (int i=0;i<xy_nranks;i++){
	     fprintf(fh,"\"%d\":1 ",xy_ranks[i]);
             if (i!=xy_nranks-1)
                 fprintf(fh,",");
         }
         fprintf(fh,"},\n");
	 fclose(fh);
      } 
    }
    MPI_Barrier(MPI_COMM_WORLD);

  }

  MPI_Barrier(MPI_COMM_WORLD);


  int xz_nranks=Grid.procs[0]*Grid.procs[2];
  int xz_ranks[xz_nranks];  

  for (int nsl = 0; nsl<Physics.slice[i_sl_xz]; nsl++){

    if ( (Grid.beg[1] <= Physics.xz_lev[nsl]+Grid.gbeg[1] ) and 
         (Grid.end[1] >= Physics.xz_lev[nsl]+Grid.gbeg[1] )){

      MPI_Gather(&Run.rank,1,MPI_INT,xz_ranks,1,MPI_INT,0,XZ_COMM);
      if (xz_rank == 0){
         fh=fopen("slice_ranks.txt","a");     
	 fprintf(fh,"\"%s_%04d\": {","xz_slice",Physics.xz_lev[nsl]);
         for (int i=0;i<xz_nranks;i++){
             fprintf(fh,"\"%d\":1 ",xz_ranks[i]);
             if (i!=xz_nranks-1)
                 fprintf(fh,",");
         }
         fprintf(fh,"},\n");
	 fclose(fh);
      } 
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int yz_nranks=Grid.procs[1]*Grid.procs[2];
  int yz_ranks[yz_nranks];

  for (int nsl = 0; nsl<Physics.slice[i_sl_yz]; nsl++){

    if ( (Grid.beg[0] <= Physics.yz_lev[nsl]+Grid.gbeg[0] ) and 
         (Grid.end[0] >= Physics.yz_lev[nsl]+Grid.gbeg[0] )){
	    
       MPI_Gather(&Run.rank,1,MPI_INT,yz_ranks,1,MPI_INT,0,YZ_COMM);
       if (yz_rank == 0){
         fh=fopen("slice_ranks.txt","a");      
	 fprintf(fh,"\"%s_%04d\": {","yz_slice",Physics.yz_lev[nsl]);
         for (int i=0;i<yz_nranks;i++){
             fprintf(fh,"\"%d\":1 ",yz_ranks[i]);
             if (i!=yz_nranks-1)
                 fprintf(fh,",");
         }
         if (nsl==Physics.slice[i_sl_yz]-1){
            fprintf(fh,"}\n");
            fprintf(fh,"}\n");
         }
         else
            fprintf(fh,"},\n");
	 fclose(fh);
       }

    }
    MPI_Barrier(MPI_COMM_WORLD);
  }  

}
