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
cmake_minimum_required (VERSION 3.1)

# Coherce CMake to install the generated .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

add_library( pastixf
  src/pastix_enums.F90
  src/spmf.f90
  src/pastixf.f90 )

if ( PASTIX_INT64 )
  set_source_files_properties(
    src/pastix_enums.F90
    PROPERTIES COMPILE_DEFINITIONS "PASTIX_INT_KIND=8")
else()
  set_source_files_properties(
    src/pastix_enums.F90
    PROPERTIES COMPILE_DEFINITIONS "PASTIX_INT_KIND=4")
endif()

target_link_libraries( pastixf pastix )
install(TARGETS pastixf
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib )


add_executable( fsimple
  examples/fsimple.f90 )

target_link_libraries( fsimple
  pastixf )

add_executable( flaplacian
  examples/flaplacian.f90 )

target_link_libraries( flaplacian
  pastixf )
