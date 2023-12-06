#if defined(MURAM_HEFFTE) || defined(MURAM_HEFFTE_CPU)

#include <mpi.h>
#include "grid.H"
#include "run.H"
#include "comm_split.H"
#include <vector>
#include "ACCH.h"
#include "heffte.h"

#include <iostream>
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;
/*
  Routine follows following naming conventions regardless of data layout
  z: vertical direction
  x: fast horizontal direction (second fftw dimension)
  y: slow horizontal direction (first fftw dimension, direction of fftw-mpi decomposition)

  Currently all operations are within the top layer of processors, i.e. it is sufficient to call
  the routine only in the boundary layer. 
*/

namespace {
  #ifndef MURAM_HEFFTE_CPU
    typedef heffte::fft3d<heffte::backend::cufft> heffte_fft3d;
    typedef heffte::gpu::vector<std::complex<double>> heffte_container;
    std::vector<std::complex<double>> cpu_input;
    std::vector<std::complex<double>> cpu_output;
  #else
    typedef heffte::fft3d<heffte::backend::fftw> heffte_fft3d;
    typedef std::vector<std::complex<double>> heffte_container;
    #define cpu_input input
    #define cpu_output output
  #endif

  void update_cpu(std::vector<std::complex<double>> &cpu, heffte_container &gpu)
  {
    #ifndef MURAM_HEFFTE_CPU
      cpu = heffte::gpu::transfer().unload(gpu);
    #endif
  }

  void update_gpu(std::vector<std::complex<double>> &cpu, heffte_container &gpu)
  {
    #ifndef MURAM_HEFFTE_CPU
      gpu = heffte::gpu::transfer().load(cpu);
    #endif
  }

  heffte_fft3d *fft;

  heffte_container input;
  heffte_container output;

  std::vector<std::vector<std::complex<double>>> kernel_fft(8);
  std::vector<std::complex<double>> bz_fft;

  void x_exchange(const GridData& Grid, double* b_ext, const int hdir_x, const int hdir_y);
  void y_exchange(const GridData&, double*, const int, const int);
}

//heffte::fft3d<heffte::backend::fftw> fft(heffte::box3d<>({0,0,0},{0,0,0}), heffte::box3d<>({0,0,0},{0,0,0}), MPI_COMM_WORLD);

void potential_ext_heffte(const GridData& Grid, double *bz0, double *b_ext){

  static int ini_flag = 1;

  static double *bz0_x, *b_ext_x;
  
  int i, j, k;

  // ==============================================================
  // ******** Change this for different dimensional layout ********
  int hdir_x = 1; // fast dimension
  int hdir_y = 2; // slow dimension, direction of fftw-mpi decomposition

  MPI_Comm x_comm = YCOL_COMM;
  MPI_Comm y_comm = ZCOL_COMM;
  MPI_Comm heffte_comm = YZ_COMM;

  int x_rank = ycol_rank;
  int y_rank = zcol_rank;
  //printf("(%d,%d)\n", x_rank, y_rank);
  //===============================================================

  int gx = Grid.gsize[hdir_x];
  int gy = Grid.gsize[hdir_y];

  int lx = Grid.lsize[hdir_x];
  int ly = Grid.lsize[hdir_y];

  int gh_x = Grid.ghosts[hdir_x];
  int gh_y = Grid.ghosts[hdir_y];
  
  int gx_gh = gx+2*gh_x;
  int lx_gh = lx+2*gh_x;
  int ly_gh = ly+2*gh_y;

  int begx = Grid.beg[hdir_x] - Grid.gbeg[hdir_x];
  int begy = Grid.beg[hdir_y] - Grid.gbeg[hdir_y];

  int buf_sz    = ly*lx;
  int buf_sz_gh = ly_gh*lx_gh;
  
  if( (Grid.NDIM == 3) and (Grid.dx[hdir_x] != Grid.dx[hdir_y]) ){
	  cout << "dx=dz required for potential field bnd: " 
		  << Grid.dx[hdir_x] << ' ' << Grid.dx[hdir_y] << endl;
	  MPI_Abort(MPI_COMM_WORLD,1);
  }

  if(ini_flag){
    for(i=0;i<8*buf_sz_gh;i++) b_ext[i] = 0.0;

    int me = heffte::mpi::comm_rank(heffte_comm);
    int num_ranks = heffte::mpi::comm_size(heffte_comm);

#ifndef MURAM_HEFFTE_CPU
    heffte::gpu::device_set(ACCH::GetGPU());
#endif

    heffte::box3d<> all_indexes({0,0,0}, {gx-1,gy-1,0});
    heffte::box3d<> inbox  = heffte::box3d<>({begx, begy, 0}, {begx+lx-1, begy+ly-1, 0});
    heffte::box3d<> outbox = heffte::box3d<>({begx, begy, 0}, {begx+lx-1, begy+ly-1, 0});

    fft = new heffte_fft3d(inbox, outbox, heffte_comm);

#ifndef MURAM_HEFFTE_CPU
    input = heffte_container(fft->size_inbox());
    output = heffte_container(fft->size_outbox());
#endif
    cpu_input.resize(fft->size_inbox());
    cpu_output.resize(fft->size_outbox());
    for(i = 0; i < 8; i++) kernel_fft[i].resize(fft->size_outbox());    
    bz_fft.resize(fft->size_outbox());

   
    //*************** Analytical kernels **********************************//	
    double kx,ky,k2,kabs,scale,hx,hy,hz,e0,e1,e2,norm;
 
    const double smooth=0.25;
    const double pi=2.0*asin(1.0);

    scale=Grid.dx[0]/Grid.dx[1];
    norm=double(gx*gy);
    
    for(j=0;j<ly;j++){
      for(i=0;i<lx;i++) {
	
	kx=Grid.beg[hdir_x] - gh_x+i;
	ky=Grid.beg[hdir_y] - gh_y+j;
	
	if (kx > Grid.gsize[hdir_x]/2) kx-=Grid.gsize[hdir_x];
	if (ky > Grid.gsize[hdir_y]/2) ky-=Grid.gsize[hdir_y];
	
	kx*=2.0*pi/Grid.gsize[hdir_x];
	ky*=2.0*pi/Grid.gsize[hdir_y];
	
	k2=kx*kx+ky*ky;
	kabs=sqrt(k2);
	
	if(kabs !=0.0){
	  hx=kx/kabs;
	  hy=ky/kabs;
	} else {
	  hx=1.0;
	  hy=1.0;
	}
	hz=1.0;
	
	e0=exp(-smooth*k2)/norm;
	e1=exp(-kabs*scale-smooth*k2)/norm;
	e2=exp(-2*kabs*scale-smooth*k2)/norm;
	
	//Bx components
	kernel_fft[0][j*lx+i] = std::complex<double>(hz*e1,0.0);
	kernel_fft[1][j*lx+i] = std::complex<double>(hz*e2,0.0);
	
	//By components
	kernel_fft[2][j*lx+i] = std::complex<double>(0.0,-hx*e0);
	kernel_fft[3][j*lx+i] = std::complex<double>(0.0,-hx*e1);
	kernel_fft[4][j*lx+i] = std::complex<double>(0.0,-hx*e2);
	
	//Bz components
	kernel_fft[5][j*lx+i] = std::complex<double>(0.0,-hy*e0);
	kernel_fft[6][j*lx+i] = std::complex<double>(0.0,-hy*e1);
	kernel_fft[7][j*lx+i] = std::complex<double>(0.0,-hy*e2);
      }
    }
    //*************** Analytical kernels **********************************//	
     
    ini_flag = 0;
  }
    
  // ======== forward fft of bz ========
  for (i=0;i<buf_sz;i++)
    cpu_input[i] = bz0[i];

  update_gpu(cpu_input, input);
  fft->forward(input.data(), output.data());
  update_cpu(cpu_output, output);

  for (i=0;i<buf_sz;i++)
    bz_fft[i] = cpu_output[i];

  //  ======== convolution with kernels ========
    
  for(k=0;k<8;k++){
    for (i=0;i<buf_sz;i++){
      cpu_output[i] = std::complex<double>(bz_fft[i].real()*kernel_fft[k][i].real()-bz_fft[i].imag()*kernel_fft[k][i].imag(),
					   bz_fft[i].real()*kernel_fft[k][i].imag()+bz_fft[i].imag()*kernel_fft[k][i].real());
    }

    update_gpu(cpu_output, output);
    fft->backward(output.data(), input.data());
    update_cpu(cpu_input, input);

    for(j=0;j<ly;j++)
      for(i=0;i<lx;i++) {
	b_ext[k*buf_sz_gh+(j+gh_y)*(lx+2*gh_x)+i+gh_x] = cpu_input[j*lx+i].real();
      }
  }
 
  x_exchange(Grid,b_ext,hdir_x,hdir_y); 
  y_exchange(Grid,b_ext,hdir_x,hdir_y);
 
}

namespace {
// ======================================================================================================================================
void y_exchange(const GridData& Grid, double* b_ext, const int hdir_x, const int hdir_y){

  static int ini_flag=1;
 
  register int i,j,k,buf;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int* leftr =Grid.leftr;
  const int* rightr=Grid.rightr;
  
  static int x_sz,y_sz,y_gh,y_beg,y_end,buf_sz;
  static double *sndbuf_l , *recbuf_l, *sndbuf_r, *recbuf_r;

  if(ini_flag){  
    y_gh  = Grid.ghosts[hdir_y];
    y_beg = Grid.ghosts[hdir_y];
    y_end = Grid.lsize[hdir_y]+Grid.ghosts[hdir_y];

    x_sz  = Grid.lsize[hdir_x]+2*Grid.ghosts[hdir_x];
    y_sz  = Grid.lsize[hdir_y]+2*Grid.ghosts[hdir_y];
    
    buf_sz = x_sz*y_gh;
    
    sndbuf_l = new double[8*buf_sz];
    recbuf_l = new double[8*buf_sz];
    sndbuf_r = new double[8*buf_sz];
    recbuf_r = new double[8*buf_sz];

    ini_flag=0;
  }
  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=y_beg;j<y_beg+y_gh;j++){
      for(i=0;i<x_sz;i++){
	sndbuf_l[buf++] = b_ext[k*x_sz*y_sz+j*x_sz+i];
      }  
    }
  }

  MPI_Irecv(recbuf_r,8*buf_sz,MPI_DOUBLE,rightr[hdir_y],0,MPI_COMM_WORLD,rr); 
  MPI_Isend(sndbuf_l,8*buf_sz,MPI_DOUBLE,leftr[hdir_y], 0,MPI_COMM_WORLD,rr+1);

  buf = 0;
  for(k=0;k<8;k++){
    for(j=y_end-y_gh;j<y_end;j++){
      for(i=0;i<x_sz;i++){
	sndbuf_r[buf++] = b_ext[k*x_sz*y_sz+j*x_sz+i];
      }  
    }
  }

  MPI_Irecv(recbuf_l,8*buf_sz,MPI_DOUBLE,leftr[hdir_y], 1,MPI_COMM_WORLD,rl);   
  MPI_Isend(sndbuf_r,8*buf_sz,MPI_DOUBLE,rightr[hdir_y],1,MPI_COMM_WORLD,rl+1);

  MPI_Waitall(2,rr,st);
  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=y_end;j<y_end+y_gh;j++){
      for(i=0;i<x_sz;i++){
	b_ext[k*x_sz*y_sz+j*x_sz+i] = recbuf_r[buf++];
      }  
    }
  }

  MPI_Waitall(2,rl,st);

  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=0;j<y_gh;j++){
      for(i=0;i<x_sz;i++){
	b_ext[k*x_sz*y_sz+j*x_sz+i] = recbuf_l[buf++];
      }  
    }
  } 
}

void x_exchange(const GridData& Grid, double* b_ext, const int hdir_x, const int hdir_y){

  static int ini_flag=1;
 
  register int i,j,k,buf;

  MPI_Status st[2];
  MPI_Request rr[2],rl[2];

  const int* leftr =Grid.leftr;
  const int* rightr=Grid.rightr;
  
  static int x_sz,y_sz,x_gh,x_beg,x_end,buf_sz;
  static double *sndbuf_l , *recbuf_l, *sndbuf_r, *recbuf_r;

  if(ini_flag){  
    x_gh  = Grid.ghosts[hdir_x];
    x_beg = Grid.ghosts[hdir_x];
    x_end = Grid.lsize[hdir_x]+Grid.ghosts[hdir_x];

    x_sz  = Grid.lsize[hdir_x]+2*Grid.ghosts[hdir_x];
    y_sz  = Grid.lsize[hdir_y]+2*Grid.ghosts[hdir_y];
    
    buf_sz = y_sz*x_gh;
    
    sndbuf_l = new double[8*buf_sz];
    recbuf_l = new double[8*buf_sz];
    sndbuf_r = new double[8*buf_sz];
    recbuf_r = new double[8*buf_sz];

    ini_flag=0;
  }
  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=0;j<y_sz;j++){
      for(i=x_beg;i<x_beg+x_gh;i++){
	sndbuf_l[buf++] = b_ext[k*x_sz*y_sz+j*x_sz+i];
      }  
    }
  }

  MPI_Irecv(recbuf_r,8*buf_sz,MPI_DOUBLE,rightr[hdir_x],0,MPI_COMM_WORLD,rr); 
  MPI_Isend(sndbuf_l,8*buf_sz,MPI_DOUBLE,leftr[hdir_x], 0,MPI_COMM_WORLD,rr+1);

  buf = 0;
  for(k=0;k<8;k++){
    for(j=0;j<y_sz;j++){
      for(i=x_end-x_gh;i<x_end;i++){
	sndbuf_r[buf++] = b_ext[k*x_sz*y_sz+j*x_sz+i];
      }  
    }
  }

  MPI_Irecv(recbuf_l,8*buf_sz,MPI_DOUBLE,leftr[hdir_x], 1,MPI_COMM_WORLD,rl);   
  MPI_Isend(sndbuf_r,8*buf_sz,MPI_DOUBLE,rightr[hdir_x],1,MPI_COMM_WORLD,rl+1);

  MPI_Waitall(2,rr,st);
  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=0;j<y_sz;j++){
      for(i=x_end;i<x_end+x_gh;i++){
	b_ext[k*x_sz*y_sz+j*x_sz+i] = recbuf_r[buf++];
      }  
    }
  }

  MPI_Waitall(2,rl,st);

  
  buf = 0;
  for(k=0;k<8;k++){
    for(j=0;j<y_sz;j++){
      for(i=0;i<x_gh;i++){
	b_ext[k*x_sz*y_sz+j*x_sz+i] = recbuf_l[buf++];
      }  
    }
  } 
}

}

#endif

