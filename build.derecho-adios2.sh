#!/bin/bash -l

module purge
module load ncarenv
module load nvhpc
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module list

export NVLOCALRC=~/localrc
export ADIOS2_DIR=/glade/work/haiyingx/ADIOS2_derecho/install_newrc
export PATH=/glade/work/haiyingx/ADIOS2_derecho/install_newrc/bin:$PATH
export LD_LIBRARY_PATH=/glade/work/haiyingx/ADIOS2_derecho/install_newrc/lib64:$LD_LIBRARY_PATH


make clean
make 

#cp src/mhd3d.x /glade/derecho/scratch/haiyingx/Run_Corona_1728x1024x1024_ASD/.
#cp backup.dat /glade/derecho/scratch/haiyingx/Run_Corona_1728x1024x1024_ASD/.
