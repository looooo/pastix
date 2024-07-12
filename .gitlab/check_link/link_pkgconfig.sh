#!/usr/bin/env sh
###
#
#  @file link_pkgconfig.sh
#  @copyright 2023-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.4.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2024-07-05
#
# Check that linking with the project is ok when using pkg-config.
#
###
set -ex

static=""
pastix_path=""

while [ $# -gt 0 ]
do
    case $1 in
        --static | -static | -s )
            static="--static"
            ;;
        * )
            pastix_path=$1
            ;;
    esac
    shift
done

if [ -z $pastix_path ]
then
    echo """
usage: ./link_pkgconfig.sh path_to_pastix_install [--static|-static|-s]
   use the --static parameter if pastix is static (.a)
   env. var. CC and FC must be defined to C and Fortran90 compilers
"""
    exit 1
fi

export PKG_CONFIG_PATH=$pastix_path/lib/pkgconfig:$PKG_CONFIG_PATH

mkdir -p build
cd build

FLAGS=`pkg-config $static --cflags pastix`
if [[ "$SYSTEM" == "macosx" ]]; then
    FLAGS="-Wl,-rpath,$pastix_path/lib $FLAGS"
fi
LIBS=`pkg-config $static --libs pastix`
$CC $FLAGS ../../../example/simple.c $LIBS -o link_pastix_c
./link_pastix_c -t 2 --lap 100

FLAGS=`pkg-config $static --cflags pastixf`
if [[ "$SYSTEM" == "macosx" ]]; then
    FLAGS="-Wl,-rpath,$pastix_path/lib $FLAGS"
fi
LIBS=`pkg-config $static --libs pastixf`
$FC $FLAGS ../../../wrappers/fortran90/examples/fsimple.F90 $LIBS -o link_pastix_f
./link_pastix_f -t 2 --lap 100
