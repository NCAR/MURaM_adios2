module purge
module load ncarenv/23.06
module load nvhpc/23.1
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module load cmake

ADIOS2_DIR=/glade/derecho/scratch/haiyingx/adios2/
PYTHON_EXECUTABLE=/glade/work/haiyingx/conda-envs/derecho_sperr103/bin/python

git clone https://github.com/ornladios/ADIOS2.git
cd ADIOS2
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=$ADIOS2_DIR -DADIOS2_USE_Python=ON  -DPython3_EXECUTABLE=$PYTHON_EXECTUABLE -DPython_FIND_STRATEGY=LOCATION -DFLEX_EXECUTABLE=/usr/bin/flex -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON ..

make -j 4
make install
