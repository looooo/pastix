###
#
#  @copyright 2017-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.1
#  @author Mathieu Faverge
#  @date 2018-07-16
#
###
if ( PASTIX_WITH_MPI )
  if ( PASTIX_WITH_FORTRAN )
    message(STATUS "--- Fortran wrapper is disabled due to MPI support activated")
    set( PASTIX_WITH_FORTRAN OFF )
  endif()
endif()

if (PASTIX_WITH_FORTRAN)
  add_subdirectory( fortran90 )
endif()

if (BUILD_SHARED_LIBS)
  if ( PASTIX_WITH_MPI )
    message(STATUS "--- Python wrapper is disabled due to MPI support activated")
  else()
  add_subdirectory( python )
  endif()
else()
  message(STATUS "--- Python wrapper is disabled with static libraries")
endif()
