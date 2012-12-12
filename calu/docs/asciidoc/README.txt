PLASMA README
=============

*****************************************************************************
Univ. of Tennessee, Univ. of California Berkeley and Univ. of Colorado Denver

     __________ ____       _____    _________   _____      _____
     \______   \    |     /  _  \  /   _____/  /     \    /  _  \
      |     ___/    |    /  /_\  \ \_____  \  /  \ /  \  /  /_\  \
      |    |   |    |___/    |    \/        \/    Y    \/    |    \
      |____|   |_______ \____|__  /_______  /\____|__  /\____|__  /
                       \/       \/        \/         \/         \/

Parallel Linear Algebra Software for Multicore Architectures

http://icl.cs.utk.edu/plasma/
*****************************************************************************

Purpose of PLASMA
-----------------

The main purpose of PLASMA is to address the performance shortcomings of the
http://www.netlib.org/lapack/[LAPACK] and
http://www.netlib.org/scalapack/[ScaLAPACK]
libraries on multicore processors and multi-socket systems of multicore processors
and their inability to efficiently utilize accelerators such as Graphics Processing Units (GPUs).
PLASMA provides routines to solve dense general systems of linear equations,
symmetric positive definite systems of linear equations and linear least squares problems,
using LU, Cholesky, QR and LQ factorizations.
Real arithmetic and complex arithmetic are supported in both single precision and double precision.

PLASMA has been designed to supercede LAPACK and ScaLAPACK, principally by restructuring the software
to achieve much greater efficiency, where possible, on modern computers based on multicore processors.
PLASMA also relies on new or improved algorithms.
Currently, however, PLASMA does not serve as a complete replacement of LAPACK due to limited functionality.
Specifically, PLASMA does not support band matrices and does not solve eigenvalue and singular value problems.
Also, PLASMA does not replace ScaLAPACK as software for distributed memory computers, since it only supports
shared-memory machines.

Where to Find More Information
------------------------------

The main repository for PLASMA documentation is the distribution ./docs directory.
The directory contains important documents such as the Users' Guide and the Reference Manual.
PLASMA documentation is also available online on the PLASMA website: http://icl.cs.utk.edu/plasma/.
For installation instructions please refer to the
http://icl.cs.utk.edu/projectsfiles/plasma/html/InstallationGuide.html[Installation Guide].
In addition, the http://icl.cs.utk.edu/plasma/forum/[PLASMA User Forum] can be used
to post general questions and comments as well as to report technical problems.

Important Information about BLAS and LAPACK
-------------------------------------------

=== Optimized BLAS are Critical for Performance ===

It is absolutely critical for performance to use PLASMA in conjunction
with an optimized implementation of the Basic Linear Algebra Subroutines
(BLAS) library.
Such implementations are usually provided by the processor manufacturer
and are usually available free of charge for non-profit use,
such as academic research. Examples include:

* The VecLib from Apple,
* The AMD Core Math Library (ACML),
* The Math Kernel Library (MKL) from Intel,
* The Engineering and Scientific Software Library (ESSL) from IBM.

Open-source alternatives also exist, such as:

* http://web.tacc.utexas.edu/~kgoto/[Goto BLAS],
* http://math-atlas.sourceforge.net/[Automatically Tuned Linear Algebra Software (ATLAS)].

As the last resort, the FORTRAN implementation of BLAS
from http://www.netlib.org/blas/[Netlib] can be used (often referred to as _reference BLAS_).
However, since Netlib BLAS are completely unoptimized, PLASMA with Netlib BLAS
will deliver correct numerical results, but no performance whatsoever.

For comprehensive installation instructions please refer to the
http://icl.cs.utk.edu/projectsfiles/plasma/html/InstallationGuide.html[Installation Guide].
However, if you decide to install manually and edit the installation scripts then there is one
important issue to keep in mind. Modern optimized BLAS are not stand-alone libraries but instead
are bundled with additional software, primarily LAPACK. This makes it necessary to use proper
linking flags. Commonly these are -L (for path to library) and -l (for library name):

 -L/usr/lib -lblas

This will only pull in symbols from the BLAS library that are needed by PLASMA. The commonly
made mistake is to just specify the full path to the library:

 /usr/lib/libblas.a

This method will work for simple cases but it forces the linker to include the whole
contents of the BLAS library rather than just pull in the missing symbols. Aside from
making the binary executable files larger, this method will easily cause name clashes as the same
symbol name might be included multiple times.

=== Multithreading within BLAS Must be Disabled ===

Many Basic Linear Algebra Subroutines (BLAS) implementations exploit parallelism within BLAS
through multithreading.
PLASMA, however, utilizes BLAS for high performance implementations of single-core operations
(often referred to as _kernels_)
and exploits parallelism at the algorithmic level above the level of BLAS.
For that reason, PLASMA must not be used in conjunction with a multithreaded BLAS,
as this is likely to create more execution threads than actual cores.
The phenomenon, known as oversubscribing of cores, will completely annihilate PLASMA's performance
due to devastating impact on the operation of cache memories for dense linear algebra workloads.

PLASMA needs to be linked with a sequential BLAS library or a multithreaded BLAS library
with multithreading disabled.
To learn how to link with sequential MKL, please consult the
http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/[Intel Math Kernel Library Link Line Advisor];
To learn how to link with sequential ATLAS, please consult the
http://math-atlas.sourceforge.net/errata.html#LINK[ATLAS errata].
Alternatively, PLASMA can be linked with a multithreaded BLAS library with multithreading disabled.
Typically, disabling of multithreading can be done by setting the appropriate environment variable
from the command prompt, for instance:

 export OMP_NUM_THREADS=1
 export MKL_NUM_THREADS=1
 export GOTO_NUM_THREADS=1
 export VECLIB_MAXIMUM_THREADS=1

Currently, this disables the simultaneous use of PLASMA for some of the user's functions and vendor
BLAS for others. This cannot be easily remedied without standardization of interoperability rules
of multithreaded libraries. One consolation is that PLASMA already delivers fast parallel implementations
of all Level 3 BLAS routines and there is virtually no benefit from parallelization of Level 1 and 2 BLAS
routines on current generation of multicore platforms due to memory contention.

=== PLASMA Software Stack ===

PLASMA requires the following software packages to be installed in the system prior to PLASMA's installation:
BLAS, CBLAS, LAPACK and Netlib LAPACK C Wrapper.
CBLAS and the components of LAPACK required by PLASMA are commonly bundled with BLAS.
If this is not the case Netlib implementation of CBLAS and Netlib LAPACK can be used.
The C interface to LAPACK has not been standardized yet and the LAPACK C Wrapper from Netlib has to be used
for the time being.

BLAS is the set of Basic Linear Algebra Subprograms described in the previous section.
An unoptimized reference implementation of BLAS is available from Netlib at
http://www.netlib.org/blas/[http://www.netlib.org/blas/].
As mentioned before, it is critical for performance that optimized implementation of BLAS is used
instead of Netlib BLAS.

CBLAS is the C language interface to BLAS available with most implementations of BLAS.
A reference implementation from Netlib is also available at
http://www.netlib.org/blas/blast-forum/cblas.tgz[http://www.netlib.org/blas/blast-forum/cblas.tgz].
Since CBLAS is only a set of wrappers to the actual BLAS, the CBLAS from Netlib can be used
without any adverse effects on performance.

LAPACK is a large package of linear algebra routines for a wide range of problems.
PLASMA uses only a tiny portion of LAPACK, which is also commonly bundled with BLAS distributions.
If this is not the case, the complete LAPACK distribution from Netlib can be used, which is available
at http://www.netlib.org/lapack/[http://www.netlib.org/lapack/].
Although vendor LAPACK routines can be more optimized than those from Netlib, there should be no
adverse performance effects of using Netlib LAPACK, since PLASMA only relies on LAPACK for implementing
some of its sequential kernels.

The user can point PLASMA's installer to all the components already installed in the system.
For all the missing components Netlib equivalents will be installed.
The installer can also be forced to disregard any software already installed in the system and
use the Netlib packages instead.

=== Thou Shalt Not Mix Compilers ===

For a given processor, the user can have different compilers at his disposal.
For instance, GNU, PGI and Intel compilers are available for Intel processors.
Different compilers can have slightly different Application Binary Interfaces (ABIs)
and mixing compilers is generally a bad idea.
User's code and the PLASMA library should be compiled with the same compiler, and so should be BLAS,
CBLAS, LAPACK and LAPACK C Wrapper, if a source distribution is used.
If a binary distribution of the BLAS is used, the correct version has to be chosen (the one providing the right ABI).
For Intel processors, the
http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/[Intel Math Kernel Library Link Line Advisor]
can be used to assist with the choice.

=== Linking FORTRAN code with C Code ===

Currently PLASMA library does not contain any FORTRAN code any more.
FORTRAN is only used in PLASMA's testing suite derived from the one of LAPACK (/testing/lin/).
Because of that, neither a FORTRAN compiler nor the FORTRAN standard library has to be involved
in compiling PLASMA and linking it with C code. The following paragraph is preserved in this README
only because it is a very useful piece of information for novice users who are forced to mix FORTRAN and C.

If FORTRAN code is mixed with C code, the FORTRAN standard library has to be included.
Sometimes it can be accomplished by simply putting the standard FORTRAN library at the and of the link line,
e.g., "-lgfortran" when using GCC.
Alternatively, FORTRAN compiler can be used for linking. This will accomplish the same effect automatically.
However, the Intel IFORT compiler, when used for linking, assumes that the main program is in FORTRAN and links
_for_main.o_ into the application.
This provides the linker with two main() functions (one created by the user and one inserted by the build system),
which is cause a linker error. To prevent this from happening, the "-nofor_main" link option has to be given.

Fortran 90 Interfaces
---------------------

It is now possible to call PLASMA from modern Fortran, making use of the Fortran 2003
C interoperability features.

The benefits of using the Fortran 90 interfaces over the old-style Fortran 77 interfaces are:
        * Compile-time argument checking.
        * Native and transparent handling of pointer arguments - arrays, descriptors and handles.
        * A clean interface between Fortran and C.

In order to build the Fortran 90 interfaces add the following to your make.inc file: PLASMA_F90 = 1

To call PLASMA via the interfaces, 'Use PLASMA' in your Fortran code.

Arguments such as descriptors and handles required by the PLASMA tiled and asynchronous interfaces
are passed as type(c_ptr), which is part of the Fortran 2003 ISO C bingings module
(so you will also need to 'Use iso_c_binding').

For the LAPACK-style interfaces, arrays should be passed in as normal.

Four examples of using the Fortran 90 interfaces are given, which show how to use the module,
call auxiliary functions such as initializing PLASMA and setting options, perform tasks such
as allocating workspace and translating between layouts, and calling a computational routine:
example_sgebrd.f90            - single precision real bi-diagonal reduction using LAPACK-syle
                                interface.
example_dgetrs_tile_async.f90 - double precision real factorizaion followed by linear solve
                                using the tiled, asynchronous interface.
example_cgeqrf_tile.f90 - single precision complex QR factorization using the tiled interface.
example_zgetrf_tile.f90 - double precision complex LU factorization using the tiled interface.


The interfaces can be found in the 'control' directory:

plasma_f90.f90  plasma_cf90.F90  plasma_df90.F90  plasma_dsf90.F90
plasma_sf90.F90  plasma_zcf90.F90  plasma_zf90.F90

* Please check the subroutine wrappers (following the 'contains' statement in each module) to see
the interfaces for the routines to call from your Fortran.

A Note on Running on NUMA Systems
---------------------------------
PLASMA is a software package for shared memory systems, both Symmetric Multi-Processors (SMP) and Non-Uniform
Memory Access (NUMA) systems. PLASMA does not detect the type of the system and does not take any specific
actions in that respect. PLASMA's performance may be poor on NUMA systems if matrices are not distributed
among multiple memory nodes. The current wisdom is to use "numactl --interleave=all" when running
an application that uses PLASMA.

A Note on Running PLASMA and OpenMP
---------------------------------
PLASMA currently binds the existing threads to specific cores in order to optimize data locality during computation.
As a consequence, if you add your first OpenMP section after the call
to PLASMA_Init, all threads created by the OpenMP section will be sons
of the main thread binded on the core 0 and will thus be binded on the
same core. To avoid this problem, it is recommended to first
have an OpenMP section to create the threads even if it does not do anything
and then place the call to PLASMA_Init.
If your OpenMP section is after the call to PLASMA_Finalize, it
should not be a problem, since version 2.4.1 will unbind threads after
this call.

License Information
-------------------

PLASMA is a software package provided by University of Tennessee, University of California, Berkeley
and University of Colorado, Denver.
PLASMA's license is a BSD-style permissive free software license
(properly called modified BSD).
It allows proprietary commercial use, and for the software released under the license
to be incorporated into proprietary commercial products.
Works based on the material may be released under a proprietary license
as long as PLASMA's license requirements are maintained,
as stated in the LICENSE file, located in the main directory of the PLASMA distribution.
In contrast to copyleft licenses, like the GNU General Public License,
PLASMA's license allows for copies and derivatives of the source code
to be made available on terms more restrictive than those of the original license.

Publications
------------

A number of technical reports were written during the development of PLASMA
and published as http://www.netlib.org/lapack/lawns/downloads/[LAPACK  Working  Notes]
by the  University  of Tennessee. Almost all of these reports later appeared as journal articles.
To make a reference to PLASMA you can cite the following publications:

***************************************
_Emmanuel Agullo, Alfredo Buttari, Jack Dongarra, Mathieu Faverge, Bilel Hadri,
Azzam Haidar, Jakub Kurzak, Julien Langou, Hatem Ltaief, Piotr Luszczek, Asim YarKhan_ +
*http://icl.cs.utk.edu/projectsfiles/plasma/pdf/users_guide.pdf[PLASMA Users' Guide]* +
_Electrical Engineering and Computer Science Department_ +
_Univesity of Tennessee_
***************************************

***************************************
_Alfredo Buttari, Julien Langou, Jakub Kurzak, Jack Dongarra_ +
*A class of parallel tiled linear algebra algorithms for multicore architectures* +
_Parallel Computing 35 (2009) 38-53_ +
_http://dx.doi.org/10.1016/j.parco.2008.10.002[DOI: 10.1016/j.parco.2008.10.002]_
***************************************

***************************************
_Emmanuel Agullo, Jim Demmel, Jack Dongarra, Bilel Hadri, Jakub Kurzak, Julien Langou,
Hatem Ltaief, Piotr Luszczek, Stanimire Tomov_ +
*Numerical linear algebra on emerging architectures: The PLASMA and MAGMA projects* +
_2009 Journal of Physics: Conference Series 180 012037_ +
_http://dx.doi.org/10.1088/1742-6596/180/1/012037[DOI: 10.1088/1742-6596/180/1/012037]_
***************************************

Funding
-------

The PLASMA project is funded in part by the National Science Foundation, Department of Energy, Microsoft, and the MathWorks.
