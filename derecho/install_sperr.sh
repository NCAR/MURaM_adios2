module purge
module load ncarenv
module load nvhpc
module load cuda
module load craype
module load cray-mpich
module load ncarcompilers
module load cray-libsci
module load cmake
module load conda

conda env create -f derecho_sperr1031.yml
conda activate derecho_sperr1031
export CFLAGS="-noswitcherror $CFLAG"
env mpicc=`which mpicc` pip3 install mpi4py
pip install *whl --force-reinstall

export PYTHON_PATH=`which python`
