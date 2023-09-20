###
#
#  @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.0
#  @author Mathieu Faverge
#  @author Tony Delarue
#  @date 2023-02-06
#
###
cmake_minimum_required (VERSION 3.1)

if ( NOT PASTIX_WITH_FORTRAN )
  return()
endif()

if ( PASTIX_WITH_MPI AND NOT MPI_Fortran_HAVE_F08_MODULE)
  message(FATAL_ERROR "PASTIX_WITH_FORTRAN and PASTIX_WITH_MPI are both ON but MPI_Fortran_HAVE_F08_MODULE is empty.\n"
                      "Please disable MPI or Fortran.")
endif()

# Coherce CMake to install the generated .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include )

set( pastixf_sources
  src/pastix_f2c.c
  src/pastixf_bindings.f90
  src/pastixf_enums.F90
  src/pastixf_functions.f90
  src/pastixf_interfaces.f90
  src/pastixf.f90 )

add_library( pastixf
  ${pastixf_sources} )

target_include_directories( pastixf PRIVATE
  $<BUILD_INTERFACE:${PASTIX_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PASTIX_BINARY_DIR}>
  $<BUILD_INTERFACE:${PASTIX_SOURCE_DIR}/common> )
target_include_directories(pastixf INTERFACE
  $<INSTALL_INTERFACE:include>
  )

set_target_properties( pastixf PROPERTIES VERSION ${PASTIX_VERSION} )
set_target_properties( pastixf PROPERTIES SOVERSION ${PASTIX_VERSION_MAJOR} )

target_link_libraries( pastixf PRIVATE pastix SPM::spmf )

if ( PASTIX_INT64 )
  set_source_files_properties( src/pastixf_enums.F90 PROPERTIES COMPILE_DEFINITIONS "PASTIX_INT_KIND=c_int64_t" )
else()
  set_source_files_properties( src/pastixf_enums.F90 PROPERTIES COMPILE_DEFINITIONS "PASTIX_INT_KIND=c_int32_t" )
endif()

if ( PASTIX_WITH_MPI )
  target_link_libraries( pastixf PUBLIC MPI::MPI_Fortran)
endif()

# export target pastixf
# ---------------------
install(EXPORT pastixfTargets
        NAMESPACE PASTIX::
        DESTINATION lib/cmake/pastix
        )

# Installation
# ------------
install(TARGETS pastixf
        EXPORT pastixfTargets
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Add examples
# ------------
set (EXAMPLES
  fsimple.F90
  flaplacian.F90
  fstep-by-step.F90
  fmultidof.F90
  )

# List of run types
# -----------------
set( RUNTYPE shm )
if (PASTIX_WITH_MPI)
  list( APPEND RUNTYPE mpi )
endif()

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)

  add_executable( ${_name_we} examples/${_file} )
  add_binary_to_completion( ${_name_we} )
  add_documented_files( examples/${_file} )

  target_link_libraries(${_name_we} PRIVATE
    pastixf SPM::spmf )
  # Need this to use mpi and initialize it
  if ( PASTIX_WITH_MPI )
    target_compile_definitions( ${_name_we} PUBLIC "PASTIX_WITH_MPI=1" )
  endif()

  install(TARGETS ${_name_we}       RUNTIME DESTINATION examples )
  install(FILES   examples/${_file}         DESTINATION examples )

  foreach( version ${RUNTYPE} )
    unset( exe )
    if( version STREQUAL "shm" )
      set( exe "")
    endif()
    if( version STREQUAL "mpi" )
      set( exe ${MPIEXEC_EXECUTABLE} -np 4 --oversubscribe )
    endif()

    add_test(fortran_${version}_${_name_we} ${exe} ./${_name_we} )
  endforeach()

endforeach()

#
# Add the fmultilap example using OpenMP only in shared memory
#
set( _file fmultilap.F90 )

get_filename_component(_name_we ${_file} NAME_WE)
add_executable(${_name_we} examples/${_file})
add_binary_to_completion( ${_name_we} )
add_documented_files( examples/${_file} )

target_link_libraries( ${_name_we} PRIVATE
  pastixf SPM::spmf )

if ( PASTIX_WITH_MPI )
  target_compile_definitions( ${_name_we} PRIVATE "PASTIX_WITH_MPI=1" )
endif()

install(TARGETS ${_name_we}       RUNTIME DESTINATION examples )
install(FILES   examples/${_file}         DESTINATION examples )
install(FILES   examples/test_seq.in      DESTINATION examples )
install(FILES   examples/test_mt.in       DESTINATION examples )

# Add OpenMP if available
find_package(OpenMP)
if (OpenMP_Fortran_FOUND)
  target_link_libraries( ${_name_we} PRIVATE OpenMP::OpenMP_Fortran )
endif()

# Add specific test
include(FindUnixCommands)

set( _seq "${CMAKE_CURRENT_SOURCE_DIR}/examples/test_seq.in")
set( _mt  "${CMAKE_CURRENT_SOURCE_DIR}/examples/test_mt.in" )
set( _exe "${CMAKE_CURRENT_BINARY_DIR}/${_name_we}"         )

if (BASH)
  add_test( NAME fortran_shm_${_name_we}_seq COMMAND ${BASH} -c "${_exe} < ${_seq}" )
  add_test( NAME fortran_shm_${_name_we}_mt  COMMAND ${BASH} -c "${_exe} < ${_mt}"  )
endif()

#
# End of fmultilap
#
