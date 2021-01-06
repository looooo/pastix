###
#
#  @copyright 2017-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.1.0
#  @author Mathieu Faverge
#  @date 2020-01-29
#
###

if (PASTIX_WITH_FORTRAN)
  add_subdirectory( fortran90 )
endif()

if (BUILD_SHARED_LIBS)
  add_subdirectory( python )
  add_subdirectory( julia )
else()
  message(STATUS "--- Python and Julia wrappers are disabled with static libraries")
endif()
