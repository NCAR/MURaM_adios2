#!/bin/bash -l

#PBS -N MURbui 
#PBS -A UCSU0085 
#PBS -q main@gusched01
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=50GB:ngpus=1
##PBS -l gpu_type=a100
#PBS -l walltime=00:10:00
#PBS -e build.err 
#PBS -o build.out 

module purge
module load ncarenv
module load nvhpc
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module list

#export NVLOCALRC=~/localrc
#export ADIOS2_DIR=/glade/work/haiyingx/ADIOS2_derecho/install_newrc
#export PATH=/glade/work/haiyingx/ADIOS2_derecho/install_newrc/bin:$PATH
#export LD_LIBRARY_PATH=/glade/work/haiyingx/ADIOS2_derecho/install_newrc/lib64:$LD_LIBRARY_PATH

export ADIOS2_DIR=/glade/derecho/scratch/haiyingx/ADIOS2_derecho/install
export PATH=$ADIOS2_DIR/bin:$PATH
export LD_LIBRARY_PATH=$ADIOS2_DIR/lib64:$LD_LIBRARY_PATH

make clean
make 

#cp src/mhd3d.x /glade/derecho/scratch/haiyingx/Run_Corona_1728x1024x1024_ASD/.
#cp backup.dat /glade/derecho/scratch/haiyingx/Run_Corona_1728x1024x1024_ASD/.
