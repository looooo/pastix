VERSION="$1"

#OpenMPI (uncomment the following lines if OpenMPI is compiled from sources and set OPENMPI_DIR path)
#export PATH=$OPENMPI_DIR/bin:$PATH
#export LD_RUN_PATH=$OPENMPI_DIR/lib:$LD_RUN_PATH
#export LD_LIBRARY_PATH=$OPENMPI_DIR/lib:$LD_LIBRARY_PATH
#export INCLUDE_PATH=$OPENMPI_DIR/include:$INCLUDE_PATH
#export PKG_CONFIG_PATH=$OPENMPI_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

#StarPU
export LD_LIBRARY_PATH=$STARPU_DIR/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$STARPU_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
#PaRSEC
if [ "$VERSION" == "mpi" ]
then
    export PARSEC_DIR=/builds/install/parsec-mpi
else
    export PARSEC_DIR=/builds/install/parsec-shm
fi
export PATH=$PARSEC_DIR/bin:$PATH
export LD_LIBRARY_PATH=$PARSEC_DIR/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$PARSEC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
#Scotch
export LD_LIBRARY_PATH=$SCOTCH_DIR/lib:$LD_LIBRARY_PATH
#Metis
export LD_LIBRARY_PATH=$METIS_DIR/lib:$LD_LIBRARY_PATH
