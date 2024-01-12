# MURaM_ADIOS2 Install Guide

## Last Update 01/11/2024

### Enviromentatl Variables

```bash
export CONDA_ENV_DIR=/global/homes/b/bz186/conda_envs/muram_adios
export ADIOS2_INSTALL_DIR=/global/homes/b/bz186/apps/gnu-11.2.0/adios2-2.10.0rc1
export HEFFTE_INSTALL_DIR=/global/homes/b/bz186/apps/nvhpc-22.7/heffte-2.4.0
```

1. Install sperr

    ```bash
    module load conda
    conda env create -f derecho/derecho_sperr1031.yml -p $CONDA_ENV_DIR
    conda activate $CONDA_ENV_DIR
    ```

    According to the [NERSC Documentation](https://docs.nersc.gov/development/languages/python/parallel-python/), install `mpi4y` on Perlmutter requires `PrgEnv-gnu` module.

    ```bash
    module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
    MPICC="cc -shared" pip install --force-reinstall --no-cache-dir --no-binary=mpi4py mpi4py

    pip install derecho/numcodecs-0.1.dev690+dirty-cp39-cp39-linux_x86_64.whl
    ```
    ignore the error message if you see "Successfully installed" in the bottom:
    ```
    "ERROR: pip's dependency resolver does not currently
    take into account all the packages that are installed. This behaviour is
    the source of the following dependency conflicts.
    zarr 2.16.1 requires numcodecs>=0.10.0, but you have numcodecs
    0.1.dev690+dirty which is incompatible.
    Successfully installed entrypoints-0.4 numcodecs-0.1.dev690+dirty
    numpy-1.26.1" 
    ```

2. Install ADIOS2

    ADIOS2 can be built with `PrgEnv-gnu` on Perlmutter

    ```bash
    module load cmake
    git clone https://github.com/ornladios/ADIOS2.git ADIOS2
    cd ADIOS2
    mkdir build && cd build
    cmake .. \
    -DCMAKE_INSTALL_PREFIX=$ADIOS2_INSTALL_DIR \
    -DADIOS2_USE_Python=ON \
    -DPython_EXECUTABLE=${ADIOS2_INSTALL_DIR}/bin/python \
    -DFLEX_EXECUTABLE=/usr/bin/flex \
    -DADIOS2_USE_Fortran=OFF \
    -DADIOS2_USE_MPI=ON

    make -j8
    make install
    ```

3. Install heffte

    heffte can be built with `PrgEnv-nvhpc` on Perlmutter. heffte needs `cray-fftw` as its backend. It can also enable `cuda/cufft` backend.

    ```bash
    module swap PrgEnv-gnu PrgEnv-nvhpc
    wget https://github.com/icl-utk-edu/heffte/archive/refs/tags/v2.4.0.tar.gz
    tar xf v2.4.0.tar.gz
    cd heffte-2.4.0/
    module load cray-fftw/3.3.10.3
    mkdir build && cd build
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=$HEFFTE_INSTALL_DIR \
    -DHeffte_ENABLE_AVX=ON \
    -DHeffte_ENABLE_FFTW=ON \
    -DFFTW_ROOT=/opt/cray/pe/fftw/3.3.10.3/x86_milan \
    -DHeffte_ENABLE_CUDA=ON \
    -DCUDA_cufft_LIBRARY=/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64/libcufft.so

    make -j8
    make install
    ```

4. Install MURaM

    MURaM is built with `PrgEnv-nvhpc` on Perlmutter. Need to modify `Make_defs` to add `-L` for cuda math library due to the special cuda math lib settings on the Perlmutter

    ```bash
    export PATH=$PATH:${ADIOS2_INSTALL_DIR}/bin
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ADIOS2_INSTALL_DIR}/lib64
    ```

    add `CUDAMATHLIBDIR` to `Make_defs` and prepend `CUDAMATHLIBDIR` to `LIBS` before `CUDALIB`:
    ```bash
    CUDAMATHLIBDIR = -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64/

    LIBS = $(MPILIBDIR) $(MPIOLIB) $(MPILIB) \
           $(HEFFTELIBDIR) $(HEFFTELIB) \
           $(FFTWLIBDIR) $(FFTWLIB) \
           $(CUDALIBDIR) $(CUDAMATHLIBDIR) $(CUDALIB) \
           $(MASSLIB) -lm \
           $(ADIOS2_LIB_DIR)
    ```