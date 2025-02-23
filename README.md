# PaStiX: A sparse direct solver

PaStiX (Parallel Sparse matriX package) is a scientific library that provides a
high performance parallel solver for very large sparse linear systems based on
direct methods.  Numerical algorithms are implemented in single or double
precision (real or complex) using LLt, LDLt and LU with static pivoting (for non
symmetric matrices having a symmetric pattern).
This solver also provides some low-rank compression methods to reduce the memory footprint and/or the time-to-solution.

## Get PaStiX

To use last development state of PaStiX, please clone the master
branch. Note that PaStiX contains two `git submodule` for **spm** and **morse_cmake**.
To get sources please use these commands:

    # if git version >= 1.9
      git clone --recursive https://gitlab.inria.fr/solverstack/pastix.git
      cd pastix
    # else
      git clone https://gitlab.inria.fr/solverstack/pastix.git
      cd pastix
      git submodule init
      git submodule update

Previous releases of PaStiX are hosted on the
[https://files.inria.fr/pastix/releases/](https://files.inria.fr/pastix/releases/) for now.
Future releases will be available on this gitlab project.

## Available Features

* Features available in shared memory with POSIX threads for multicores architectures.

|                         | Seq   | Static | Dyn   | StarPU    | PaRSEC    |
|-------------------------|-------|--------|-------|-----------|-----------|
| POTRF (Cholesky)        | FR/LR | FR/LR  | FR/LR | FR/LR/GPU | FR/LR/GPU |
| PXTRF (LL^t for complex)| FR/LR | FR/LR  | FR/LR | FR/LR/GPU | FR/LR/GPU |
| HETRF (LDL^h)           | FR/LR | FR/LR  | FR/LR | FR/LR/GPU | FR/LR/GPU |
| SYTRF (LDL^t)           | FR/LR | FR/LR  | FR/LR | FR/LR/GPU | FR/LR/GPU |
| GETRF (LU)              | FR/LR | FR/LR  | FR/LR | FR/LR/GPU | FR/LR/GPU |
| TRSM                    | FR/LR | FR/LR  | FR/LR | FR/LR     | -         |
| DIAG                    | FR/LR | FR/LR  | FR/LR | FR/LR     | -         |

* Features available in hybrid shared/distributed memory with MPI between processes, and POSIX threads within a process.

|                          | Seq   | Static | Dyn   | StarPU | PaRSEC |
|--------------------------|-------|--------|-------|--------|--------|
| POTRF (Cholesky)         | FR/LR | FR/LR  | FR/LR | FR/GPU | FR/GPU |
| PXTRF (LL^t for complex) | FR/LR | FR/LR  | FR/LR | FR/GPU | FR/GPU |
| HETRF (LDL^h)            | FR/LR | FR/LR  | FR/LR | FR/GPU | FR/GPU |
| SYTRF (LDL^t)            | FR/LR | FR/LR  | FR/LR | FR/GPU | FR/GPU |
| GETRF (LU)               | FR/LR | FR/LR  | FR/LR | FR/GPU | FR/GPU |
| TRSM                     | FR/LR | FR/LR  | FR/LR | FR     | -      |
| DIAG                     | FR/LR | FR/LR  | FR/LR | FR     | -      |

* FR means Full-Rank computations without compression techniques
* LR means Low-Rank compression technique to reduce the memory footprint and/or the time-to-solution
* *WARNING* GPU kernels are not available on compressed supernodes
* *WARNING* Schur complement is not available with MPI yet.

## Documentation

The latest Doxygen documentation is available [here](http://solverstack.gitlabpages.inria.fr/pastix).

The [main steps](http://solverstack.gitlabpages.inria.fr/pastix/group__pastix__users.html) and [parameters](http://solverstack.gitlabpages.inria.fr/pastix/group__pastix__api.html) of the solver are described. Some [examples](http://solverstack.gitlabpages.inria.fr/pastix/group__pastix__examples.html) are also provided.

## Installation

### Build and install with CMake

PaStiX can be built using [CMake](https://cmake.org/). This
installation requires to have some library dependencies already
installed on the system:

* BLAS (MKL, OpenBlas, ...) and CBLAS (sequential version required)
* LAPACK and LAPACKE (sequential version required, with TMG enabled for testing)
* HWLOC (highly recommended)
* SCOTCH (optional)
* METIS (optional)
* STARPU runtime support (optional)
* PARSEC runtime support (optional)
* CUDA/CuBLAS to enable GPU functionality with runtime support (optional)
* EZTRACE to enable tracing support (optional)
* Python and Fortran compiler for wrappers and examples (optional)

For instance, on debian-like systems, dependencies can be installed with the following command:

      sudo apt-get install cmake gcc gfortran libhwloc-dev libscotch-dev libopenblas-dev liblapacke-dev python-numpy

The main options to configure the PaStiX configuration build are:

* Classic cmake options:
  * CMAKE_BUILD_TYPE: Debug, RelWithDebInfo, Release, MinSizeRel; we recommend to use the Release, or RelWithDebInfo, for performance.
  * CMAKE_INSTALL_PREFIX: Specify the prefix directory to install the library
  * BUILD_SHARED_LIBS=[OFF]: Enable the shared libraries build. This option needs to be enabled for the Python wrapper.
* Integer type:
  * PASTIX_INT64[=ON]: Enable/disable int64_t for integer arrays.
* Ordering libraries:
  * Ordering libraries must match the integer type chosen for integer arrays in PaStiX
  * PASTIX_ORDERING_SCOTCH[=ON]: Enable/Disable the support of the Scotch library to compute the ordering.
  * PASTIX_ORDERING_METIS[=OFF]: Enable/Disable the support of the Metis library to compute the ordering. Metis 5.1 is required.
* External schedulers:
  * PASTIX_WITH_PARSEC[=OFF]: Enable/disable the PaRSEC runtime support. Require to install PaRSEC tag pastix-_releasenumber_ (mymaster for master branch) from the repository <https://bitbucket.org/mfaverge/parsec> that includes a few patches on top of the original PaRSEC runtime system. PaRSEC needs to be compiled with option -DPARSEC_WITH_DEVEL_HEADERS=ON.
  * PASTIX_WITH_STARPU[=OFF]: Enable/disable the StarPU runtime support. Require to install StarPU 1.3.
* Distributed memory:
  * PASTIX_WITH_MPI=[OFF]: Enable/disable distributed memory support (See above for details). If used with the PaRSEC library, MPI should be enabled or disabled in both libraries.*
* External SpM library:
  * PASTIX_WITH_EXTERNAL_SPM=[OFF]: Enable/disable the use of an external SpM library in favor of the internal one.
* Documentation:
  * BUILD_DOCUMENTATION[=OFF] to enable the Doxygen documentation generation.

## Get involved

### Reporting an issue

We strongly recommend all users to use the issue tracker to report any problems with the software, or for any feature request. We will try our best to answer them in a short time frame.

### Contributions

<https://gitlab.inria.fr/solverstack/pastix/blob/master/CONTRIBUTING.md>

## Authors

The following people contribute or contributed to the development of PaStiX:

* Mathieu Faverge, PI
* Pierre Ramet, PI
* David Goudin
* Mathias Hastaran
* Pascal Henon
* Xavier Lacoste
* François Pellegrini
* Grégoire Pichon, Low-rank solver
* Florent Pruvost, CMake and Spack
* Theophile Terraz
* Tony Delarue, MPI implementation
* Brieuc Nicolas, mixed precision
* Alycia Lisito

If we forgot your name, please let us know that we can fix that mistake.

## Citing PaStiX

Feel free to use the following publications to reference PaStiX:

* Original paper that initiated PaStiX:
  * Pascal Hénon, Pierre Ramet, Jean Roman. Pascal Hénon, Pierre Ramet, Jean Roman. PaStiX: A High-Performance Parallel Direct Solver for Sparse Symmetric Definite Systems. Parallel Computing, Elsevier, 2002, 28 (2), pp.301--321. [INRIA HAL](https://hal.inria.fr/inria-00346017)
* Parallel incomplete factorization implemented in PaStiX:
  * Pascal Hénon, Pierre Ramet, Jean Roman. On finding approximate supernodes for an efficient ILU(k) factorization. Parallel Computing, Elsevier, 2008, 34, pp.345--362. [INRIA HAL](https://hal.inria.fr/inria-00346018)
* Reordering strategy for blocking optimization in PaStiX:
  * Grégoire Pichon, Mathieu Faverge, Pierre Ramet, Jean Roman. Reordering Strategy for Blocking Optimization in Sparse Linear Solvers. SIAM Journal on Matrix Analysis and Applications, Society for Industrial and Applied Mathematics, 2017, SIAM Journal on Matrix Analysis and Applications, 38 (1), pp.226 - 248. [INRIA HAL](https://hal.inria.fr/hal-01485507v2)
* On the use of low rank approximations in PaStiX:
  * Grégoire Pichon, Eric Darve, Mathieu Faverge, Pierre Ramet, Jean Roman. Sparse supernodal solver using block low-rank compression: Design, performance and analysis. International Journal of Computational Science and Engineering, Inderscience, 2018, 27, pp.255 - 270. [10.1016/J.JOCS.2018.06.007](http://dx.doi.org/10.1016/J.JOCS.2018.06.007) [Inria HAL](https://hal.inria.fr/hal-01824275)

## Licence

<https://gitlab.inria.fr/solverstack/pastix/blob/master/LICENSE>
