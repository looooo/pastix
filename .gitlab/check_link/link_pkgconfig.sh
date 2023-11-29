#!/bin/sh
# usage: ./link_pkgconfig.sh path_to_pastix_install [--static]
# use the --static parameter if pastix is static (.a)
# env. var. CC and FC must be defined to C and Fortran90 compilers
set -ex

export PKG_CONFIG_PATH=$1/lib/pkgconfig:$PKG_CONFIG_PATH

mkdir -p build
cd build

FLAGS=`pkg-config $2 --cflags pastix`
if [[ "$SYSTEM" == "macosx" ]]; then
  FLAGS="-Wl,-rpath,$1/lib $FLAGS"
fi
LIBS=`pkg-config $2 --libs pastix`
$CC $FLAGS ../../../example/simple.c $LIBS -o link_pastix_c
./link_pastix_c -t 2 --lap 100

FLAGS=`pkg-config $2 --cflags pastixf`
if [[ "$SYSTEM" == "macosx" ]]; then
  FLAGS="-Wl,-rpath,$1/lib $FLAGS"
fi
LIBS=`pkg-config $2 --libs pastixf`
$FC $FLAGS ../../../wrappers/fortran90/examples/fsimple.F90 $LIBS -o link_pastix_f
./link_pastix_f -t 2 --lap 100
