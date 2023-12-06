#download libfabric, now it is installed by sys admin already
#wget https://github.com/ofiwg/libfabric/releases/download/v1.13.1/libfabric-1.13.1.tar.bz2 
#tar xjvf libfabric-1.13.1.tar.bz2
#cd libfabric-1.13.1
#module load cuda/11.7.1
#module load gcc
#./configure --prefix=/glade/work/haiyingx/libfabric-1.13.1/install --enable-cuda-dlopen --enable-gdrcopy-dlopen --enable-sockets=yes --enable-verbs=yes --enable-tcp=yes --enable-rxm=yes --enable-shm=yes --disable-psm --disable-psm2 --disable-psm3 --with-cuda=/glade/u/apps/common/22.12/spack/opt/spack/cuda/11.7.1


module load gcc
git clone https://github.com/ornladios/ADIOS2.git
cd ADIOS2
git checkout 08e6bc524ae019a83be03bed2ad2f9de5bc3bff1

module load conda
conda-env create -f mpich_sperr_de.yml
conda activate mpich_sperr_de
pip install /glade/work/haiyingx/numcodecs_sperr/dist/numcodecs-0.0.0-cp39-cp39-linux_x86_64.whl --force-reinstall
export CFLAGS="-noswitcherror $CFLAG"
which mpicc
env MPICC=/glade/u/apps/derecho/23.06/spack/opt/spack/ncarcompilers/1.0.0/gcc/12.2.0/gj6c/bin/mpi/mpicc pip3 install mpi4py --no-cache-dir

cd ADIOS2
mkdir build
cd build
source ../install_derecho.sh


on cpu nodes instead of login node compiled mpi4py
export CFLAGS="-noswitcherror $CFLAG"
env MPICC=/glade/u/apps/derecho/23.06/spack/opt/spack/ncarcompilers/1.0.0/nvhpc/23.1/rqst/bin/mpi/mpicc pip3 install mpi4py --no-cache-dir
