#!/bin/bash -l

#PBS -N muram_io.adios2_1stripe_2m_1644 
#PBS -A	AHAO0002
#PBS -q main
#PBS -l select=16:ncpus=64:mpiprocs=8:ngpus=4
#PBS -l walltime=00:45:00

#cd $cwd
export TMPDIR=/glade/derecho/scratch/$USER/temp
export ADIOS2_DIR=/glade/work/haiyingx/ADIOS2_derecho/install_newrc/
mkdir -p $TMPDIR

module purge
module load ncarenv
#module load nvhpc/23.5
module load nvhpc/23.5
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module load darshan-runtime
module list


#export LIBS_HOME=/glade/work/rempel/MURaM_gust/FFTW_LIBS/
#export FFTW3_HOME=$LIBS_HOME/fftw_gpu
#export HEFFTE_HOME=$LIBS_HOME/heffte_gpu
#export PATH=$FFTW3_HOME/bin:$HEFFTE_HOME/bin:$PATH
#export LD_LIBRARY_PATH=$FFTW3_HOME/lib:$HEFFTE_HOME/lib:$LD_LIBRARY_PATH
#export NVHPC_CUDA_HOME=/glade/u/apps/common/23.04/spack/opt/spack/cuda/11.7.1/
#stack
export LIBS_HOME=/glade/work/haiyingx/MURaM_derecho/
export FFTW3_HOME=$LIBS_HOME/fftw
export HEFFTE_HOME=$LIBS_HOME/heffte
export PATH=$FFTW3_HOME/bin:$HEFFTE_HOME/bin:$PATH
export LD_LIBRARY_PATH=$FFTW3_HOME/lib:$HEFFTE_HOME/lib:$LD_LIBRARY_PATH
export NVHPC_CUDA_HOME=/glade/u/apps/common/23.04/spack/opt/spack/cuda/11.7.1/

export MURaM_FFTW_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1
export CRAY_ACCEL_TARGET=nvidia80
export CUDA_VISIBLE_DEVICES=0,1,2,3
#export SstVerbose=5

#export PCAST_COMPARE=rel=7,patchall,summary,file=pcast.dat
#export PCAST_COMPARE=create,file="/glade/gust/scratch/cmille73/muram/nvhpc_2211/MURaM_main/TEST/Test_3D/golden.dat"
#export PCAST_COMPARE=compare,rel=7,patchall,summary,file="/glade/gust/scratch/cmille73/muram/nvhpc_2211/MURaM_main/TEST/Test_3D/golden.dat"


#export MPICH_MPIIO_STATS=1

#export MPICH_MPIIO_HINTS_DISPLAY=1

#export MPICH_MPIIO_TIMERS=1
export MPI_SHEPHERD=true
export DARSHAN_MOD_ENABLE="DXT_POSIX,DXT_MPIIO"
### Source the test environment
echo "job starts"
echo $cwd
echo "PATH: "
echo $PATH
echo "LD_LIBRARY_PATH: "
echo $LD_LIBRARY_PATH
lfs setstripe 3D -S 33m -c 48 
lfs setstripe 2D -S 6m -c 4 
EXECUTABLE="./mhd3d.x"
echo "executable is $EXECUTABLE"
./clean
lfs getstripe -d 3D
lfs getstripe -d 2D

ulimit -s unlimited
nvidia-smi


#export DARSHAN_LOG_DIR_PATH=/glade/derecho/scratch/haiyingx/muram_darshan_log_mpiio/
##DARSHAN_ENABLE_NONMPI=1 
#LD_PRELOAD="$NCAR_LDFLAGS_DARSHAN_RUNTIME/libdarshan.so" mpiexec -n 64 -ppn 4 --cpu-bind list:0,16,32,48 get_local_rank $EXECUTABLE & #  > writer_out
mpiexec -n 64 -ppn 4 --cpu-bind list:0,16,32,48 get_local_rank ./mhd3d_de_asd.x & #  > writer_out
#cuda-memcheck ./mhd3d.x > MURaM.out

module load conda
conda activate derecho_sperr1031
export PYTHONPATH=$ADIOS2_DIR/lib/python3.9/site-packages/adios2
export LD_LIBRARY_PATH=$PWD:/opt/cray/pe/ucx/1.14.0/ucx/lib/:$LD_LIBRARY_PATH
#DARSHAN_ENABLE_NONMPI=1 LD_PRELOAD="$NCAR_LDFLAGS_DARSHAN_RUNTIME/libdarshan.so" mpiexec -n 64 -ppn 4 --cpu-bind list:15,31,47,63 python SstReader_mpi.py > reader_out
mpiexec -n 64 -ppn 4 --cpu-bind list:15,31,47,63 python SstReader_mpi.py > reader_out
~                                                                    
