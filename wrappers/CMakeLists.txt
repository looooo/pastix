###
#
#  @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @date 2017-05-22
#
###

if (PASTIX_WITH_FORTRAN)
  add_subdirectory( fortran90 )
endif()
add_subdirectory( python )
