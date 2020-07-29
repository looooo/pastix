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
  liblapack-dev \
  liblapacke-dev \
  libmetis-dev \
  nvidia-cuda-dev \
  nvidia-cuda-toolkit \
  zlib1g-dev	      \
  python3-numpy	      \
  python3-mpi4py      \
  python3-scipy	      \
  && \
  apt-get autoremove -y

# change the default shell to be bash
SHELL ["/bin/bash", "-c"]

# Create directory where to install our specific libraries
RUN mkdir -p /builds/install

# Install OpenMPI
RUN cd $HOME && \
    wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.2.tar.bz2 && \
    tar xvjf openmpi-4.0.2.tar.bz2 && \
    cd openmpi-4.0.2/ && \
    ./configure --prefix=/builds/install/openmpi && \
    make -j 4 && \
    make install && \
    cd $HOME && \
    rm openmpi-4.0.2/ -rf
ENV OPENMPI_DIR=/builds/install/openmpi

ENV PATH=$OPENMPI_DIR/bin:$PATH
ENV LD_RUN_PATH=$OPENMPI_DIR/lib:$LD_RUN_PATH
ENV LD_LIBRARY_PATH=$OPENMPI_DIR/lib:$LD_LIBRARY_PATH
ENV INCLUDE_PATH=$OPENMPI_DIR/include:$INCLUDE_PATH
ENV PKG_CONFIG_PATH=$OPENMPI_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

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
    wget http://gforge.inria.fr/frs/download.php/file/35817/SimGrid-3.13.tar.gz && \
    tar xvzf SimGrid-3.13.tar.gz && \
    cd SimGrid-3.13 && \
    cmake . -DCMAKE_INSTALL_PREFIX=/builds/install/simgrid -DCMAKE_BUILD_TYPE=Debug && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm SimGrid-3.13/ SimGrid-3.13.tar.gz -rf
ENV SIMGRID_DIR=/builds/install/simgrid

# Install StarPU
RUN cd $HOME && \
    wget http://starpu.gforge.inria.fr/files/starpu-1.3.2/starpu-1.3.2.tar.gz && \
    tar xvzf starpu-1.3.2.tar.gz && \
    cd starpu-1.3.2/ && \
    ./configure --prefix=/builds/install/starpu --enable-debug && \
    make -j5 && \
    make install && \
    ./configure --prefix=/builds/install/starpu-simgrid --enable-debug --disable-cuda --disable-opencl --with-fxt=/usr/lib/x86_64-linux-gnu/ --enable-simgrid --with-simgrid-dir=/builds/install/simgrid && \
    make -j5 && \
    make install && \
    cd $HOME && \
    rm starpu-1.3.2/ starpu-1.3.2.tar.gz -rf
ENV STARPU_DIR=/builds/install/starpu
ENV STARPUSIMGRID_DIR=/builds/install/starpu-simgrid

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

# Install Scotch
ENV SCOTCH_DIR=/builds/install/scotch
RUN cd $HOME && \
    wget https://gforge.inria.fr/frs/download.php/file/38187/scotch_6.0.9.tar.gz && \
    tar xvzf scotch_6.0.9.tar.gz && \
    cd scotch_6.0.9/src && \
    mkdir -p $SCOTCH_DIR  && \
    cp Make.inc/Makefile.inc.i686_pc_linux3 Makefile.inc && \
    sed -i 's/__restrict//' Makefile.inc && \
    sed -i 's/CCD		= gcc/CCD		= mpicc/' Makefile.inc && \
    echo "prefix = $SCOTCH_DIR" >> Makefile.inc && \
    make scotch ptscotch && \
    make install && \
    cd $HOME && \
    rm scotch_6.0.9/ -rf

# Set gitlab environment
RUN chown -R gitlab:gitlab /builds
USER gitlab

# specific env. var.
ENV OPENMPI_DIR=/builds/install/openmpi
ENV QUARK_DIR=/builds/install/quark
ENV PARSEC_DIR=/builds/install/parsec-shm
ENV STARPU_DIR=/builds/install/starpu
ENV SIMGRID_DIR=/builds/install/simgrid
ENV STARPUSIMGRID_DIR=/builds/install/starpu-simgrid
ENV SPM_DIR=/builds/install/spm
ENV SCOTCH_DIR=/builds/install/scotch

ENV PATH=$OPENMPI_DIR/bin:$PATH
ENV LD_RUN_PATH=$OPENMPI_DIR/lib:$LD_RUN_PATH
ENV LD_LIBRARY_PATH=$OPENMPI_DIR/lib:$LD_LIBRARY_PATH
ENV INCLUDE_PATH=$OPENMPI_DIR/include:$INCLUDE_PATH
ENV PKG_CONFIG_PATH=$OPENMPI_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

WORKDIR /builds
