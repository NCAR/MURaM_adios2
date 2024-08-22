module load gcc
module load conda


conda env create -f derecho_sperr1031.yml
conda activate derecho_sperr1031
export CFLAGS="-noswitcherror $CFLAG"
env mpicc=`which mpicc` pip3 install mpi4py
pip install *whl --force-reinstall --no-deps

export PYTHON_PATH=`which python`
