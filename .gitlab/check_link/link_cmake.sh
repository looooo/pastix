#!/bin/sh
# usage: ./link_cmake.sh path_to_pastix_install
set -ex
cmake -B build -DCMAKE_PREFIX_PATH=$1
cmake --build build --verbose
ctest --test-dir build --verbose
