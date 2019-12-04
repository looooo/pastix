VERSION="$1"

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
