module purge
module load ncarenv
module load nvhpc
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module load cmake

ADIOS2_DIR=$PWD/adios2/
PYTHON_EXECUTABLE=/glade/work/haiyingx/miniconda3-py39/envs/mpich_sperr1003_de/bin/python

git clone https://github.com/ornladios/ADIOS2.git
cd ADIOS2
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=$ADIOS2_DIR -DADIOS2_USE_Python=ON  -DPython_EXECUTABLE=$PYTHON_EXECUTABLE -DFLEX_EXECUTABLE=/usr/bin/flex -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON ..

make -j 4
make install
