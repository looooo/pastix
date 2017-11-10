PaStiX: A sparse direct solver
==============================

PaStiX (Parallel Sparse matriX package) is a scientific library that provides a
high performance parallel solver for very large sparse linear systems based on
direct methods.  Numerical algorithms are implemented in single or double
precision (real or complex) using LLt, LDLt and LU with static pivoting (for non
symmetric matrices having a symmetric pattern).
This solver also provides some low-rank compression methods to reduce the memory footprint and/or the time-to-solution.

Get PaStiX
----------

To use last development state of PaStiX, please clone the master
branch. Note that PaStiX contains a `git submodule` **morse_cmake**.
To get sources please use these commands:

    # if git version >= 1.9
      git clone --recursive git@gitlab.inria.fr:solverstack/pastix.git
      cd pastix
    # else
      git clone git@gitlab.inria.fr:solverstack/pastix.git
      cd pastix
      git submodule init
      git submodule update

Last releases of PaStiX are hosted on the
[gforge.inria.fr](https://gforge.inria.fr/frs/?group_id=186) for now.
Future releases will be available on this gitlab project.

Available Features
------------------

|                         | Seq    | Static | Dyn    | StarPU  | PaRSEC     |
|-------------------------|--------|--------|--------|---------|------------|
| POTRF (Cholesky)        | SHM/LR | SHM/LR | -      | SHM/GPU | SHM/LR/GPU |
| PXTRF (LL^t for complex)| SHM/LR | SHM/LR | -      | SHM/GPU | SHM/LR/GPU |
| HETRF (LDL^h)           | SHM/LR | SHM/LR | -      | SHM/GPU | SHM/LR/GPU |
| SYTRF (LDL^t)           | SHM/LR | SHM/LR | -      | SHM/GPU | SHM/LR/GPU |
| GETRF (LU)              | SHM/LR | SHM/LR | -      | SHM/GPU | SHM/LR/GPU |
| TRSM                    | SHM/LR | SHM/LR | -      | -       | -          |
| DIAG                    | SHM/LR | SHM/LR | -      | -       | -          |

 * SHM means Shared Memory using POSIX theads for multicores architectures
 * LR means (block) Low-Rank compression technique to reduce the memory footprint and/or the time-to-solution
   * *WARNING* GPU kernels are not available when compression techniques are enabled
 * MPI is not available yet and will come with 6.1.0
 * StarPU support is available only when low-rank compression is disabled.
 * GPUs kernels are in the code but not exploited yet, we are targeting for a simpler scheduling that would allow everyone to get correct performance out of the box in final 6.0.0

Documentation
-------------

A temporary link to the Doxygen [documentation](http://pastix.gforge.inria.fr/doxygen/html/group__pastix__users.html). (we are waiting for availability of pages functionality in gitlab...)

Installation
------------

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
       * PASTIX_WITH_PARSEC[=OFF]: Enable/disable the PaRSEC runtime support. Require to install PaRSEC tag pastix-_releasenumber_ (mymaster for master branch) from the repository https://bitbucket.org/mfaverge/parsec that includes a few patches on top of the original PaRSEC runtime system. PaRSEC needs to be compiled with option -DPARSEC_WITH_DEVEL_HEADERS=ON.
       * PASTIX_WITH_STARPU[=OFF]: Enable/disable the StarPU runtime support. Require to install StarPU 1.2.
   * Distributed memory:
       * PASTIX_WITH_MPI=[OFF]: Distributed memory is not supported yet in PaStiX, however you might need to enable this option if your PaRSEC library has been compiled with MPI support.
   * Documentation:
       * BUILD_DOCUMENTATION[=OFF] to enable the Doxygen documentation generation


Get involved!
---------------------

### Reporting an issue

We strongly recommend all users to use the issue tracker to report any problems with the software, or for any feature request. We will try our best to answer them in a short time frame.

### Contributions

https://gitlab.inria.fr/solverstack/pastix/blob/master/CONTRIBUTING.md

### Authors

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

If we forgot your name, please let us know that we can fix that mistake.

### Citing PaStiX

Feel free to use the following publications to reference PaStiX:

* Original paper that initiated PaStiX:
  - Pascal Hénon, Pierre Ramet, Jean Roman. Pascal Hénon, Pierre Ramet, Jean Roman. PaStiX: A High-Performance Parallel Direct Solver for Sparse Symmetric Definite Systems. Parallel Computing, Elsevier, 2002, 28 (2), pp.301--321. [INRIA HAL](https://hal.inria.fr/inria-00346017)
* Parallel incomplete factorization implemented in PaStiX:
  - Pascal Hénon, Pierre Ramet, Jean Roman. On finding approximate supernodes for an efficient ILU(k) factorization. Parallel Computing, Elsevier, 2008, 34, pp.345--362. [INRIA HAL](https://hal.inria.fr/inria-00346018)
* Reordering strategy for blocking optimization in PaStiX:
  - Grégoire Pichon, Mathieu Faverge, Pierre Ramet, Jean Roman. Reordering Strategy for Blocking Optimization in Sparse Linear Solvers. SIAM Journal on Matrix Analysis and Applications, Society for Industrial and Applied Mathematics, 2017, SIAM Journal on Matrix Analysis and Applications, 38 (1), pp.226 - 248. [INRIA HAL](https://hal.inria.fr/hal-01485507v2)
* On the use of low rank approximations in PaStiX:
  - Grégoire Pichon, Eric Darve, Mathieu Faverge, Pierre Ramet, Jean Roman. Sparse Supernodal Solver Using Block Low-Rank Compression. 18th IEEE International Workshop on Parallel and Distributed Scientific and Engineering Computing (PDSEC 2017), Jun 2017, Orlando, United States. [INRIA HAL](https://hal.inria.fr/hal-01502215)

### Licence

https://gitlab.inria.fr/solverstack/pastix/blob/master/LICENCE.txt
