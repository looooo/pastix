## Installing PaStiX with GPU support

Note that GPU support is only available through runtime system
environment. StarPU and/or PaRSEC is thus a mandatory dependency to
enable GPUs support.

** Requirements **

In addition to the requirements of the sequential, you'll need:

    * An MPI implementation (MPICH, OpenMPI...)

### Download and install OpenMPI

The Open MPI Project is an open source Message Passing Interface
implementation that is developed and maintained by a consortium
of academic, research, and industry partners.

We recommend a release >= 3.0 available on the [OpenMPI
website](https://www.open-mpi.org/).
A Ubuntu package (version 3.x) is also available on Ubuntu,
and can be download with :

```sh
sudo apt-get install openmpi-bin openmpi-common
```

To verify your downloaded version, you can launch :
```sh
mpiexec --version
```

### Download and install MPICH

MPICH is a high performance and widely portable implementation
of the Message Passing Interface (MPI) standard. MPICH and its
derivatives form the most widely used implementations of MPI in the world

We recommend a release >= 3.0 available on the [MPICH
website](https://www.mpich.org/).
A Ubuntu package (version 3.x) is also available on Ubuntu,
and can be download with :

```sh
sudo apt-get install mpich
```

To verify your downloaded version, you can launch :
```sh
mpiexec --version
```

### Compile PaStiX with MPI

Now that you have installed MPI, you should be able to
compile PaStiX with MPI support:

Go in your PaStiX source directory and create a new build
directory to build the GPU version:
```sh
mkdir build_mpi
cd build_MPI
cmake .. -DCMAKE_INSTALL_PREFIX=${PASTIX_DIR} -DPASTIX_WITH_MPI=ON
make
make install
```

Note that as before we defined a `PASTIX_DIR` variable to define the
install directory of PaStiX, and you just need to turn ON/OFF the
runtime systems based on the available one(s).

Once the compilation finished, you can setup your environment easily
by sourcing the provided file:
```sh
source $PASTIX_DIR/bin/pastix_env.sh
```
And then, you can run your favorite example:
```sh
mpiexec -n 4 simple -9 10:10:10
```

### How to use pastix with MPI


