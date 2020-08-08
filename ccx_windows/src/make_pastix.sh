#!/bin/bash
if ! [[ -d build ]]; then
    mkdir build
fi
cd build


cmake   \
        -DBLAS_DIR=~/OpenBLAS_i8 \
	-DHWLOC_DIR=~/PaStiX/hwloc_i8 \
	-DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-10.2 \
	-DCMAKE_INSTALL_PREFIX=~/PaStiX/pastix_i8 \
	-DCMAKE_BUILD_TYPE=Release \
	-DPASTIX_WITH_PARSEC=ON \
	-DPARSEC_DIR=~/PaStiX/parsec_i8 \
	-DSCOTCH_DIR=~/PaStiX/scotch_i8 \
	-DPASTIX_WITH_CUDA=ON \
	-DCUDA_DIR=/usr/local/cuda-10.2 \
	-DPASTIX_ORDERING_SCOTCH=ON \
	-DCMAKE_C_COMPILER=gcc \
	-DCMAKE_CXX_COMPILER=g++ \
	-DCMAKE_Fortran_COMPILER=gfortran \
	-DCMAKE_C_FLAGS="-fopenmp" \
        ..

make -j8
make install
