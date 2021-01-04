###
#
#  @copyright 2017-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.3
#  @author Mathieu Faverge
#  @date 2019-12-04
#
###
cmake_minimum_required (VERSION 3.1)

if ( PASTIX_WITH_MPI AND NOT MPI_Fortran_HAVE_F08_MODULE )
  return()
endif()

# Coherce CMake to install the generated .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include )

set( pastixf_sources
  src/pastix_enums.F90
  src/pastixf.f90 )

add_library( pastixf
  ${pastixf_sources} )

set( pastixf_compile_options     "" )
set( pastixf_compile_definitions "" )
set( pastixf_include_dirs        "" )

if ( PASTIX_INT64 )
  set( pastixf_compile_definitions
    "${pastixf_compile_definitions};PASTIX_INT_KIND=c_int64_t" )
else()
  set( pastixf_compile_definitions
    "${pastixf_compile_definitions};PASTIX_INT_KIND=c_int32_t" )
endif()

if ( PASTIX_WITH_MPI )
  set( pastixf_compile_options
    "${pastixf_compile_options};${MPI_Fortran_COMPILE_OPTIONS}" )
  set( pastixf_compile_definitions
    "${pastixf_compile_definitions};${MPI_Fortran_COMPILE_DEFINITIONS};PASTIX_WITH_MPI=1" )
  set( pastixf_include_dirs
    "${pastixf_include_dirs};${MPI_Fortran_INCLUDE_DIRS}" )
endif()

set_source_files_properties(
  ${pastixf_sources}
  PROPERTIES
  COMPILE_OPTIONS     "${pastixf_compile_options}"
  COMPILE_DEFINITIONS "${pastixf_compile_definitions}"
  INCLUDE_DIRECTORIES "${pastixf_include_dirs}" )

target_link_libraries( pastixf PUBLIC
  pastix spmf )

install(TARGETS pastixf
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib )

#
# Add examples
#
set (EXAMPLES
  fsimple.f90
  flaplacian.f90
  fstep-by-step.f90
  fmultidof.f90
  )

# List of run types
set( RUNTYPE shm )
if (PASTIX_WITH_MPI)
  list( APPEND RUNTYPE mpi )
endif()

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)

  add_executable(${_name_we} examples/${_file})

  target_link_libraries(${_name_we} PRIVATE
    pastixf spmf )

  if ( PASTIX_WITH_MPI )
    set_source_files_properties(
      examples/${_file}
      PROPERTIES COMPILE_DEFINITIONS "PASTIX_WITH_MPI=1")

    target_link_libraries(${_name_we} PRIVATE
      MPI::MPI_Fortran )
  endif()

  install(TARGETS ${_name_we}       RUNTIME DESTINATION examples )
  install(FILES   examples/${_file}         DESTINATION examples )

  foreach( version ${RUNTYPE} )
    unset( exe )
    if( version STREQUAL "shm" )
      set( exe "")
    endif()
    if( version STREQUAL "mpi" )
      set( exe mpiexec -np 4 --host localhost:4 )
    endif()

    add_test(${version}_fortran_${_name_we} ${exe} ./${_name_we})
  endforeach()

endforeach()

#
# Add the fmultilap example using OpenMP only in shared memory
#
set( _file fmultilap.f90 )

get_filename_component(_name_we ${_file} NAME_WE)
add_executable(${_name_we} examples/${_file})

target_link_libraries( ${_name_we} PRIVATE
  pastixf )

if ( PASTIX_WITH_MPI )
  set_source_files_properties(
    examples/${_file}
    PROPERTIES COMPILE_DEFINITIONS "PASTIX_WITH_MPI=1")

  target_link_libraries( ${_name_we} PRIVATE
    MPI::MPI_Fortran )
endif()

install(TARGETS ${_name_we}       RUNTIME DESTINATION examples )
install(FILES   examples/${_file}         DESTINATION examples )
install(FILES   examples/test_seq.in      DESTINATION examples )
install(FILES   examples/test_mt.in       DESTINATION examples )

# Add OpenMP if available
find_package(OpenMP)
if (OpenMP_Fortran_FOUND)
  set_target_properties( ${_name_we}
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS} )
  target_link_libraries( ${_name_we} PRIVATE
    ${OpenMP_Fortran_LIBRARIES} )
endif()

# Add specific test
include(FindUnixCommands)

set( _seq "${CMAKE_CURRENT_SOURCE_DIR}/examples/test_seq.in")
set( _mt  "${CMAKE_CURRENT_SOURCE_DIR}/examples/test_mt.in" )
set( _exe "${CMAKE_CURRENT_BINARY_DIR}/${_name_we}"         )

if (BASH)
  add_test( NAME shm_fortran_${_name_we}_seq COMMAND ${BASH} -c "${_exe} < ${_seq}" )
  add_test( NAME shm_fortran_${_name_we}_mt  COMMAND ${BASH} -c "${_exe} < ${_mt}"  )
endif()

#
# End of fmultilap
#
