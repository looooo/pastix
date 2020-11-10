## Build and Install PaStiX with CMake

** Requirements **

PaStiX needs some libraries to be build on your system :

* BLAS (MKL, OpenBlas, ...) and CBLAS (sequential version required)
* LAPACK and LAPACKE (sequential version required, with TMG enabled for testing)
* HWLOC (highly recommended)
* SCOTCH (optional)
* METIS (optional)

### Install BLAS and LAPACK

**Install BLAS and CBLAS**

BLAS (Basic Linear Algebra Subprograms) is a specification that
prescribes a set of low-level routines for performing common
linear algebra operations. These routines were originally wrote in
Fortran, so we need the CBLAS interface to compile it in C.

```sh
sudo apt-get install libopenblas-dev
```

**Install LAPACK and LAPACKE**

LAPACK and LAPACKE are librairies that provide routines for
solving systems of linear equations. They depend on the BLAS library.

```sh
sudo apt-get install liblapacke-dev
```

### Install HWLOC

The [Portable Hardware Locality (hwloc)](https://www.open-mpi.org/projects/hwloc/)
software package provides a portable abstraction of the hierarchical topology 
of modern architectures, including NUMA memory nodes, sockets, shared
caches, cores and simultaneous multithreading. It also gathers various
system attributes such as cache and memory information as well as the 
locality of I/O devices such as network interfaces, InfiniBand HCAs or GPUs.

```sh
sudo apt-get install libhwloc-dev
```

### Install an ordering library
At last, we need to install a library for the ordering part of the factorization.
PaStiX supports both [Scotch](https://www.labri.fr/perso/pelegrin/scotch/)
and [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).

#### Install Metis

You can install it simply with the following command :

```sh
sudo apt-get install libsmetis-dev
```

#### Install Scotch

You can install it simply with the following command :

```sh
sudo apt-get install libscotch-dev
```

Or, for the 6.0 version :

```sh
sudo apt-get install libscotch-6.0
```

However, you may want compile your own Scotch library. For doing so,
download the version that you want and go in your Scotch directory.

Let's define the directory where you want to install Scotch as:
```sh
export SCOTCH_DIR=/your/path/to/install
```

Note that you can also keep the default install directory, and install
in the system using root privileges. We do not recommend this latest version.

To build scotch, you have to copy the Make.inc corresponding to your architecture :
```sh
cd scotch_x.x.x/src
cp Make.inc/Makefile.inc.xxxx_xxx_xxx Makefile.in
make [-DINTSIZE32|-DINTSIZE64]
make prefix=$SCOTCH_DIR install
```
Note that either you choose INTSIZE32 or INTSIZE64, it will
define the value of -DPASTIX_INT64=[ON|OFF] for PaStiX.

Scotch is now installed on your system but may not be available in
your environment.

One way to check if it is set correctly is to verify that the following command
return the path to the freshly installed Scotch:
```sh
which gmap
```

If it does not find it, you can execute the following lines and try
again.

On Linux:
```sh
export PATH=$PATH:$SCOTCH_DIR/bin
export PKG_CONFIG_PATH=:$PKG_CONFIG_PATH:$SCOTCH_DIR/lib/pkgconfig
export LD_RUN_PATH=$LD_RUN_PATH:$SCOTCH_DIR/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SCOTCH_DIR/lib
export INCLUDE_PATH=$INCLUDE_PATH:$SCOTCH_DIR/include
```

On MacOS:
```sh
export PATH=$PATH:$SCOTCH_DIR/bin
export PKG_CONFIG_PATH=:$PKG_CONFIG_PATH:$SCOTCH_DIR/lib/pkgconfig
export INCLUDE_PATH=$INCLUDE_PATH:$SCOTCH_DIR/include
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SCOTCH_DIR/lib:
```

Note that these lines can (should) be added to your environment file
such as `.bashrc` to set it by default.

### Build and Install PaStiX

Now that you have installed all the required libraries, you should be
able to compile PaStiX:

Go in your PaStiX source directory and create a new build directory to
build the shared memory version:
```sh
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PASTIX_DIR} -DPASTIX_INT64=[ON|OFF] \
make
make install
```

Note that as before we defined a `PASTIX_DIR` variable to define the
install directory of PaStiX.

Once the compilation finished, you can setup your environment easily
by sourcing the provided file:
```sh
source $PASTIX_DIR/bin/pastix_env.sh
```

And then, you can run your favorite example:
```sh
simple -9 10:10:10
```

You setup is ready to play with Pastix. Please refer to section _How to
to use PaStiX_ to get as best results as possible.
