#
# Default options for the gitlab CI test configurations
#
set( BUILD_SHARED_LIBS ON CACHE BOOL "" )

set( CMAKE_INSTALL_PREFIX "$ENV{PWD}/../install-${PASTIX_CI_VERSION}" CACHE PATH "" )
set( CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "" )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "" )

set( CMAKE_BUILD_TYPE "Debug" )
set( CMAKE_C_FLAGS -O0 -g )

set( PASTIX_ORDERING_METIS ON CACHE BOOL "" )
set( PASTIX_INT64 OFF CACHE BOOL "" )

option(MORSE_ENABLE_WARNING  "Enable warning messages"        ON)
option(MORSE_ENABLE_COVERAGE "Enable flags for coverage test" ON)

set( PASTIX_WITH_STARPU ON CACHE BOOL "" )
set( PASTIX_WITH_PARSEC ON CACHE BOOL "" )

if ( "${PASTIX_CI_VERSION}" STREQUAL "doc" )
  set( BUILD_DOCUMENTATION ON CACHE BOOL "" )
  set( PASTIX_WITH_MPI     ON CACHE BOOL "" )
elseif ( "${PASTIX_CI_VERSION}" STREQUAL "mpi" )
  set( PASTIX_WITH_MPI ON  CACHE BOOL "" )
else()
  set( PASTIX_WITH_MPI OFF CACHE BOOL "" )
endif()
