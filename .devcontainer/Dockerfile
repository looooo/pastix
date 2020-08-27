FROM hpclib/gitlab-ubuntu

USER root

# Installing as root: docker images are usually set up as root.
# Since some autotools scripts might complain about this being unsafe, we set
# FORCE_UNSAFE_CONFIGURE=1 to avoid configure errors.
ENV FORCE_UNSAFE_CONFIGURE=1
ENV DEBIAN_FRONTEND noninteractive
ENV SPACK_ROOT=/usr/local

# Install common packages
RUN apt-get update && apt-get install -y \
  gnupg2 \
  libacl1-dev \
  libboost-dev \
  libboost-context-dev \
  libfftw3-dev \
  libfxt-dev \
  libgtg-dev \
  libhypre-dev \
  liblapack-dev \
  liblapacke-dev \
  libmetis-dev \
  libmumps-dev \
  libopenmpi-dev \
  libpapi-dev \
  libparmetis-dev \
  libpfm4-dev \
  libptscotch-dev \
  libscalapack-mpi-dev \
  libsuitesparse-dev \
  nvidia-cuda-dev \
  nvidia-cuda-toolkit \
  ocl-icd-opencl-dev \
  petsc-dev \
  slepc-dev \
  zlib1g-dev \
  python3-numpy \
  python3-mpi4py \
  python3-scipy && \
  apt-get autoremove -y

# Install ViTE packages
RUN apt-get update && apt-get install -y \
    qtbase5-dev qt5-qmake qttools5-dev \
    libglm-dev \
    libglew-dev

# change the default shell to be bash
SHELL ["/bin/bash", "-c"]

# Create directory where to install our specific libraries
RUN mkdir -p /builds/install

# Install QUARK
RUN cd $HOME && \
    git clone https://github.com/ecrc/quark && \
    cd quark/ && \
    sed -i -e "s#prefix=.*#prefix=/builds/install/quark#g" make.inc && \
    sed -i -e "s#CFLAGS=.*#CFLAGS= -g -O2 -DADD_ -fPIC#g" make.inc && \
    make && \
    make install && \
    cd $HOME && \
    rm quark/ -rf
ENV QUARK_DIR=/builds/install/quark

# Install PARSEC
RUN cd $HOME && \
    git clone https://bitbucket.org/mfaverge/parsec.git && \
    cd parsec && \
    git checkout mymaster && \
    git submodule update && \
    mkdir -p build-shm && \
    cd build-shm && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/builds/install/parsec-shm -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DPARSEC_GPU_WITH_CUDA=OFF && \
    make -j5 && \
    make install && \
    cd .. && \
    mkdir -p build-mpi && \
    cd build-mpi && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/builds/install/parsec-mpi -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_DIST_WITH_MPI=ON && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm parsec/ -rf
ENV PARSEC_DIR=/builds/install/parsec-shm

# Install SimGrid
RUN cd $HOME && \
    wget https://github.com/simgrid/simgrid/releases/download/v3.25/SimGrid-3.25.tar.gz && \
    tar xvzf SimGrid-3.25.tar.gz && \
    cd SimGrid-3.25 && \
    cmake . -DCMAKE_INSTALL_PREFIX=/builds/install/simgrid -DCMAKE_BUILD_TYPE=Debug -Denable_msg=ON && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm SimGrid-3.25/ SimGrid-3.25.tar.gz -rf
ENV SIMGRID_DIR=/builds/install/simgrid

# Install StarPU
RUN cd $HOME && \
    wget http://starpu.gforge.inria.fr/files/starpu-1.3.4/starpu-1.3.4.tar.gz && \
    tar xvzf starpu-1.3.4.tar.gz && \
    cd starpu-1.3.4/ && \
    ./configure --prefix=/builds/install/starpu --enable-debug && \
    make -j5 && \
    make install && \
    ./configure --prefix=/builds/install/starpu-simgrid --enable-debug --disable-cuda --disable-opencl --with-fxt=/usr/lib/x86_64-linux-gnu/ --enable-simgrid --with-simgrid-dir=/builds/install/simgrid && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm starpu-1.3.4/ starpu-1.3.4.tar.gz -rf
ENV STARPU_DIR=/builds/install/starpu
ENV STARPUSIMGRID_DIR=/builds/install/starpu-simgrid

# Install Chameleon
RUN cd $HOME && \
    git clone --recursive https://gitlab.inria.fr/solverstack/chameleon.git && \
    cd chameleon/ && \
    source ./.gitlab-ci-env.sh && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/builds/install/chameleon -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_ENABLE_DOC=ON && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm chameleon/ -rf
ENV CHAMELEON_DIR=/builds/install/chameleon

# Install HQR
RUN cd $HOME && \
    git clone --recursive https://gitlab.inria.fr/solverstack/hqr.git && \
    cd hqr/ && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/builds/install/hqr -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm hqr/ -rf
ENV HQR_DIR=/builds/install/hqr

# Install Fabulous
RUN cd $HOME && \
    git clone --recursive https://gitlab.inria.fr/solverstack/fabulous.git && \
    cd fabulous/ && \
    git checkout release/1.0 && \
    git submodule update && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/builds/install/fabulous -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DFABULOUS_BUILD_C_API=ON -DFABULOUS_BUILD_Fortran_API=ON && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm fabulous/ -rf
ENV FABULOUS_DIR=/builds/install/fabulous

# Install SPM
RUN cd $HOME && \
    git clone --recursive https://gitlab.inria.fr/solverstack/spm.git && \
    cd spm/ && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/builds/install/spm -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DSPM_INT64=OFF && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm spm/ -rf
ENV SPM_DIR=/builds/install/spm

# Install PaStiX
RUN cd $HOME && \
    git clone --recursive https://gitlab.inria.fr/solverstack/pastix.git && \
    cd pastix/ && \
    mkdir build && \
    cd build && \
    PARSEC_DIR=/builds/install/parsec-shm && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/builds/install/pastix -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DPASTIX_WITH_STARPU=ON -DPASTIX_WITH_PARSEC=ON -DPASTIX_INT64=OFF && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm pastix/ -rf
ENV PASTIX_DIR=/builds/install/pastix

# Install Maphys
RUN cd $HOME && \
    git clone --recursive https://gitlab.inria.fr/solverstack/maphys/maphys.git && \
    cd maphys/ && \
    git checkout develop && \
    git submodule update && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/builds/install/maphys -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON -DPASTIX_VER_6=ON && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm maphys/ -rf

# Install HYPRE
RUN cd $HOME && \
    wget https://github.com/hypre-space/hypre/archive/v2.19.0.tar.gz && \
    tar xf v2.19.0.tar.gz && \
    rm v2.19.0.tar.gz && \
    cd hypre-2.19.0/src && \
    ./configure --prefix=/builds/install/hypre --enable-debug --enable-shared && \
    make -j5 install

# Set gitlab environment
RUN chown -R gitlab:gitlab /builds
USER gitlab

# specific env. var.
ENV QUARK_DIR=/builds/install/quark
ENV PARSEC_DIR=/builds/install/parsec-shm
ENV STARPU_DIR=/builds/install/starpu
ENV SIMGRID_DIR=/builds/install/simgrid
ENV STARPUSIMGRID_DIR=/builds/install/starpu-simgrid
ENV CHAMELEON_DIR=/builds/install/chameleon
ENV HQR_DIR=/builds/install/hqr
ENV FABULOUS_DIR=/builds/install/fabulous
ENV SPM_DIR=/builds/install/spm
ENV SCALFMM_DIR=/builds/install/scalfmm
ENV PASTIX_DIR=/builds/install/pastix
ENV MAPHYS_DIR=/builds/install/maphys
ENV PAMPA_DIR=/builds/install/pampa
ENV HYPRE_DIR=/builds/install/hypre
WORKDIR /builds
