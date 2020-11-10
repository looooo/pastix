## Installing PaStiX with MPI support

Here you have the features available in hybrid shared/distributed
memory with MPI between processes, and POSIX threads within a process.

|                         | Seq   | Static | Dyn   | StarPU | PaRSEC |
|-------------------------|-------|--------|-------|--------|--------|
| POTRF (Cholesky)        | FR    | FR     | FR    | FR     | FR     |
| PXTRF (LL^t for complex)| FR    | FR     | FR    | FR     | FR     |
| HETRF (LDL^h)           | FR    | FR     | FR    | FR     | FR     |
| SYTRF (LDL^t)           | FR    | FR     | FR    | FR     | FR     |
| GETRF (LU)              | FR    | FR     | FR    | FR     | FR     |
| TRSM                    | FR    | FR     | FR    | FR     | FR     |
| DIAG                    | FR    | FR     | FR    | FR     | FR     |

**Disclaimer**

Even if PaStiX works in distributed, it is still under development.
You can face either performance issues (without runtimes) or memory
issues (with StarPU or PaRSEC)

**Requirements**

In addition to the requirements of the sequential PaStiX, you'll need:

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

### Download and install StarPU

[StarPU](https://starpu.gitlabpages.inria.fr) is a runtime system
developed by the research team STORM at Inria Bordeaux -
Sud-Ouest. It provides a sequential task flow programming model similar
to OpenMP tasks with an additional support of the GPUs and distributed
environment.

The latest release is available
[here](https://files.inria.fr/starpu/), and we recommend a release >=
1.3.7.
A detailed explanation on how to install StarPU is directly available
on [their
website](https://files.inria.fr/starpu/testing/master/doc/html/BuildingAndInstallingStarPU.html),
however we give here a basic set of commands to install the release
from the latest release.

Go into your own working directory for StarPU.
Download the latest version and extract the archive.

```sh
wget https://files.inria.fr/starpu/starpu-1.3.7/starpu-1.3.7.tar.gz
tar xvzf starpu-1.3.7.tar.gz
```

Let's define the directory where you want to StarPU everything as:
```sh
export STARPU_DIR=/your/path/to/install
```
Note that you can also keep the default install directory, and install
in the system using root privileges. We do not recommend this latest version.

Let's create a build directory outside the source directory to
configure and build StarPU with minimal flags for PaStiX:
```sh
mkdir build
cd build
../starpu-1.3.7/configure --prefix=${STARPU_DIR} --enable-mpi
make
make install
```

StarPU is now installed on your system but may not be available in
your environment.

One way to check if it is set correctly is to verify that the following command
return the path to the freshly installed StarPU:
```sh
which starpu_calibrate_bus
```

If it does not find it, you can execute the following lines and try
again.

On Linux:
```sh
export PATH=$PATH:$STARPU_DIR/bin
export PKG_CONFIG_PATH=:$PKG_CONFIG_PATH:$STARPU_DIR/lib/pkgconfig
export LD_RUN_PATH=$LD_RUN_PATH:$STARPU_DIR/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$STARPU_DIR/lib
export INCLUDE_PATH=$INCLUDE_PATH:$STARPU_DIR/include
```

On MacOS:
```sh
export PATH=$PATH:$STARPU_DIR/bin
export PKG_CONFIG_PATH=:$PKG_CONFIG_PATH:$STARPU_DIR/lib/pkgconfig
export INCLUDE_PATH=$INCLUDE_PATH:$STARPU_DIR/include
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$STARPU_DIR/lib:
```

Note that these lines can (should) be added to your environment file
such as `.bashrc` to set it by default.

### Download and install PaRSEC

[PaRSEC](http://icl.utk.edu/parsec/) is a runtime system developed by the Innovative Computing
Laboratory at the University of Tennessee. It provides multiple
programming models to interact with such as a the _Dynamic Task
Discovery (DTD)_ model (identical to the _Sequential Task Flow (STF)_
from StarPU), as well as a _Parameterized Task Graph (PTG)_ model
which is here exploited by PaStiX to implement an
heterogeneous implementation (Distributed implementation of PaStiX on
top of PaRSEC is not yet tested and validated).

The latest release available is the one from the [BitBucket
repository](https://bitbucket.org/icldistcomp/parsec). However, due to specific requirements in the context of
PaStiX, we **solely** rely on [this dedicated
fork](https://bitbucket.org/mfaverge/parsec) with the specific
`pastix-6.0.2` tag. Be careful to compile PaStiX **only** with the
specific revision to enable MPI support.

Go into your own working directory for PaRSEC.
Clone the repository, and get the specific tag.

```sh
git clone --recursive https://bitbucket.org/mfaverge/parsec
cd parsec
git checkout pastix-6.0.2
```

Let's define the directory where you want to install PaRSEC as:
```sh
export PARSEC_DIR=/your/path/to/install
```

Note that you can also keep the default install directory, and install
in the system using root privileges. We do not recommend this latest version.

Let's create a build directory outside the source directory to
configure and build PaRSEC with minimal flags for PaStiX:
```sh
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${PARSEC_DIR} -DPARSEC_DIST_WITH_MPI=ON
make
make install
```

PaRSEC is now installed on your system but may not be available in
your environment.

One way to check if it is set correctly is to verify that the following command
return the path to the freshly installed PaRSEC compiler:
```sh
which parsec_ptgpp
```

If it does not find it, you can execute the following lines and try
again.

On Linux:
```sh
export PATH=$PATH:$PARSEC_DIR/bin
export PKG_CONFIG_PATH=:$PKG_CONFIG_PATH:$PARSEC_DIR/lib/pkgconfig
export LD_RUN_PATH=$LD_RUN_PATH:$PARSEC_DIR/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PARSEC_DIR/lib
export INCLUDE_PATH=$INCLUDE_PATH:$PARSEC_DIR/include
```

On MacOS:
```sh
export PATH=$PATH:$PARSEC_DIR/bin
export PKG_CONFIG_PATH=:$PKG_CONFIG_PATH:$PARSEC_DIR/lib/pkgconfig
export INCLUDE_PATH=$INCLUDE_PATH:$PARSEC_DIR/include
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$PARSEC_DIR/lib:
```

Note that these lines can (should) be added to your environment file
such as `.bashrc` to set it by default.

### Compile PaStiX with MPI

Now that you have installed MPI, you should be able to
compile PaStiX with MPI support:

Go in your PaStiX source directory and create a new build
directory to build the GPU version:
```sh
mkdir build_mpi
cd build_MPI
cmake .. -DCMAKE_INSTALL_PREFIX=${PASTIX_DIR} -DPASTIX_WITH_MPI=ON \
         -DPASTIX_WITH_STARPU=[ON|OFF] -DPASTIX_WITH_PARSEC=[ON|OFF]
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

PaStiX works with all the schedulers whithout runtimes in distributed
memory, but if you want the best performances possible with POSIX threads,
please use the dynamic scheduler (option **-s 4**):
```sh
mpiexec -n 4 simple -s 4 -9 10:10:10
```

In the same way, if you intent to run multiple MPI instances on one node,
make sur that the number of MPI instances per node times the number of
threads per MPI instances is less or equals the number of cores of the node.

If HWLOC gives you trouble in the recognition of the topology, you can
specify the amount of threads per node with the option **-t**.
But don't forget to contact us afterwards!