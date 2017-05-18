PaStiX: A sparse direct solver
==============================

PaStiX (Parallel Sparse matriX package) is a scientific library that provides a
high performance parallel solver for very large sparse linear systems based on 
direct methods.  Numerical algorithms are implemented in single or double 
precision (real or complex) using LLt, LDLt and LU with static pivoting (for non 
symmetric matrices having a symmetric pattern).  This solver provides also an
adaptive blockwise iLU(k) factorization that can be used as a parallel
preconditioner using approximated supernodes to build a coarser block structure
of the incomplete factors. 

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
[gforge.inria.fr](https://gforge.inria.fr/frs/?group_id=2884) for now.
Future releases will be available on this gitlab project.

Documentation
---------------------

A temporary link to the Doxygen [documentation](http://www.labri.fr/perso/ramet/pastix/out/html/index.html). (we are waiting for availability of pages functionality in gitlab...)

Installation
---------------------

### Build and install with CMake

PaStiX can be built using [CMake](https://cmake.org/). This
installation requires to have some library dependencies already
installed on the system.

Please refer to XXXXX
to get configuration information.

### Distribution of PaStiX
To get support to install a full distribution (PaStiX +
dependencies) we encourage users to use the morse branch of
**Spack**.

Please read these documentations:

* [Spack Morse](http://morse.gforge.inria.fr/spack/spack.html).
* [Section PaStiX](http://morse.gforge.inria.fr/spack/spack.html#sec-2-2).

Get involved!
---------------------

### Mailing list

TODO

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

TODO

### Licence

https://gitlab.inria.fr/solverstack/pastix/blob/master/LICENCE.txt
