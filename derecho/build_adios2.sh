module purge
module load ncarenv
module load nvhpc
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module load cmake

#export NVLOCALRC=~/localrc
#export CC=mpicc CXX=mpicxx FC=mpif90
#export CMAKE_INSTALL_PREFIX=/glade/work/haiyingx/ADIOS2_derecho/install
#cmake -DCC=mpicc -DCXX=mpicxx -DFC=mpif90 -DCMAKE_INSTALL_PREFIX=/glade/work/haiyingx/ADIOS2_derecho/install -DADIOS2_USE_Python=ON  -DPython_EXECUTABLE=/glade/work/haiyingx/conda-envs/mpich_sperr_de/bin/python -DFLEX_EXECUTABLE=/usr/bin/flex -DADIOS2_USE_MPI=ON ..
cmake -DCMAKE_INSTALL_PREFIX=/glade/derecho/scratch/haiyingx/ADIOS2_derecho/install -DADIOS2_USE_Python=ON  -DPython_EXECUTABLE=/glade/work/haiyingx/conda-envs/mpich_sperr_de/bin/python -DFLEX_EXECUTABLE=/usr/bin/flex -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON ..
#cmake -DCMAKE_INSTALL_PREFIX=/glade/work/haiyingx/ADIOS2_derecho/install -DADIOS2_USE_Python=ON  -DLIBFABRIC_ROOT=/glade/work/haiyingx/libfabric-1.13.1/install -DPython_EXECUTABLE=/glade/work/haiyingx/conda-envs/mpich_sperr/bin/python -DFLEX_EXECUTABLE=/usr/bin/flex -DADIOS2_USE_MPI=ON ..
make -j 4
make install
