###
#
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.2
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2023-12-07
#
###
cmake_minimum_required(VERSION 3.5)

project(LINK_PASTIX C Fortran)

# look for PASTIX on the system
# Hint: use CMAKE_PREFIX_PATH (env. var. or cmake var.) to the installation
# directory of PASTIX if not installed in a standard path
find_package(PASTIX REQUIRED)

if (TARGET PASTIX::pastix)
  message(STATUS "PASTIX_BIN_DIR ${PASTIX_BIN_DIR}")
  message(STATUS "PASTIX_INC_DIR ${PASTIX_INC_DIR}")
  message(STATUS "PASTIX_LIB_DIR ${PASTIX_LIB_DIR}")
else()
  message(FATAL_ERROR "Target PASTIX::pastix is not found, check your PASTIXConfig.cmake.")
endif()

set(CMAKE_BUILD_RPATH ${PASTIX_LIB_DIR})

add_executable(link_pastix_c ${CMAKE_CURRENT_SOURCE_DIR}/../../example/simple.c)
target_link_libraries(link_pastix_c PRIVATE PASTIX::pastix SPM::spm)

# Hint: pastix handles shared memory parallelism, deactivate multithreading in blas/lapack
# export OPENBLAS_NUM_THREADS=1 or export MKL_NUM_THREADS=1 depending on your blas/lapack vendor
# ./link_pastix -t 2

add_executable(link_pastix_f ${CMAKE_CURRENT_SOURCE_DIR}/../../wrappers/fortran90/examples/fsimple.F90)
target_link_libraries(link_pastix_f PRIVATE PASTIX::pastixf SPM::spmf)

# Launch executables with ctest
enable_testing()
include(CTest)
add_test(link_pastix_c link_pastix_c -t 2 --lap 100)
add_test(link_pastix_f link_pastix_f -t 2 --lap 100)
