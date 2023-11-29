#!/usr/bin/env bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

set -x

if [[ "$SYSTEM" != "windows" ]]; then
  if [[ "$SYSTEM" == "macosx" ]]; then
    if brew ls --versions scotch > /dev/null; then
      echo "Scotch is already installed with brew";
    else
      echo "Start installing Scotch with brew";
      brew install scotch;
    fi
    if brew ls --versions starpu > /dev/null; then
      echo "Starpu is already installed with brew";
    else
      echo "Start installing Starpu with brew";
      brew install --build-from-source ~/brew-repo/starpu.rb;
    fi
    # clang is used on macosx and it is not compatible with MORSE_ENABLE_COVERAGE=ON
    # to avoid the Accelerate framework and get Openblas we use BLA_PREFER_PKGCONFIG
    # we do not have parsec installed on the macosx machine
    cmake -B build -S . -DPASTIX_CI_VERSION=${VERSION} -DPASTIX_CI_BRANCH=${BRANCH} \
          -C .gitlab/ci-test-initial-cache.cmake \
          -DMORSE_ENABLE_COVERAGE=OFF -DBLA_PREFER_PKGCONFIG=ON -DPASTIX_WITH_PARSEC=OFF \
          || fatal
  else
    source ${CI_PROJECT_DIR}/.gitlab/env.sh ${VERSION}
    unset SPM_DIR # Make sure we don't use the installed spm
    if [ ! -z "$SPM_DIR" ]; then source $SPM_DIR/bin/spm_env.sh; fi
    cmake -B build -S . -DPASTIX_CI_VERSION=${VERSION} -DPASTIX_CI_BRANCH=${BRANCH} \
          -C .gitlab/ci-test-initial-cache.cmake \
          || fatal
    cp compile_commands.json compile_commands-${VERSION}.json
  fi
else
  # on windows the mpi_f08 interface is missing, see https://www.scivision.dev/windows-mpi-msys2/
  # default scotch in windows msys2 is int32
  # do not use static libraries because executables are too large and the build
  # directory can reach more than 10Go
  cmake -GNinja -B build -S . -DCMAKE_INSTALL_PREFIX=$PWD/install-${VERSION} \
        -DPASTIX_WITH_MPI=OFF -DPASTIX_INT64=OFF -DBUILD_SHARED_LIBS=ON \
        || fatal
fi
cmake --build build -j 4 || fatal
cmake --install build || fatal

# Check link to pastix
cd .gitlab/check_link/
if [[ "$SYSTEM" == "macosx" ]]; then
  export CC=clang
  if brew ls --versions pastix > /dev/null; then brew remove --force --ignore-dependencies pastix; fi
  if brew ls --versions pastix64 > /dev/null; then brew remove --force --ignore-dependencies pastix64; fi
else
  export CC=gcc
fi
export FC=gfortran
if [[ "$SYSTEM" == "linux" ]]; then
  export LIBRARY_PATH=$PWD/../../install-${VERSION}/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$PWD/../../install-${VERSION}/lib:$LD_LIBRARY_PATH
elif [[ "$SYSTEM" == "macosx" ]]; then
  export LIBRARY_PATH=$PWD/../../install-${VERSION}/lib:$LIBRARY_PATH
  export DYLD_LIBRARY_PATH=$PWD/../../install-${VERSION}/lib:$DYLD_LIBRARY_PATH
elif [[ "$SYSTEM" == "windows" ]]; then
  export PATH="/c/Windows/WinSxS/x86_microsoft-windows-m..namespace-downlevel_31bf3856ad364e35_10.0.19041.1_none_21374cb0681a6320":$PATH
  export PATH=$PWD/../../install-${VERSION}/bin:$PATH
fi
# 1) using cmake:
./link_cmake.sh $PWD/../../install-${VERSION} || fatal
# 2) using pkg-config:
./link_pkgconfig.sh $PWD/../../install-${VERSION} || fatal
rm -r build || fatal
