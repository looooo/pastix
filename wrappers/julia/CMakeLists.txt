###
#
#  @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @author Selmane Lebdaoui
#  @date 2020-07-15
#
###

# Configure enum.py
if (PASTIX_INT64)
  set(PASTIX_JULIA_INTEGER Int64)
else()
  set(PASTIX_JULIA_INTEGER Int32)
endif()

if (PASTIX_WITH_MPI)
  set(PASTIX_JULIA_MPI_ENABLED true)
else()
  set(PASTIX_JULIA_MPI_ENABLED false)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/PaStiX/src/pastix_enums.jl.in"
  "${CMAKE_CURRENT_BINARY_DIR}/PaStiX/src/pastix_enums.jl" @ONLY)

#Copy wrapper to build
if ( PASTIX_WITH_MPI )
  file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/PaStiX/Project.mpi.toml
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/PaStiX )
  file(RENAME
    ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/Project.mpi.toml
    ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/Project.toml )
else()
  file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/PaStiX/Project.shm.toml
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/PaStiX )
  file(RENAME
    ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/Project.shm.toml
    ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/Project.toml )
endif()

file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/PaStiX/src/PaStiX.jl
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/src )

# Copy examples to build
file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/simple.jl
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/step-by-step.jl
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/schur.jl
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/examples )

# Install julia package
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/Project.toml
  DESTINATION lib/julia/PaStiX )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/src/PaStiX.jl
  ${CMAKE_CURRENT_BINARY_DIR}/PaStiX/src/pastix_enums.jl
  DESTINATION lib/julia/PaStiX/src )

# Install julia examples
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/simple.jl
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/step-by-step.jl
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/schur.jl
  DESTINATION lib/julia/PaStiX/examples )

find_program( Julia_EXECUTABLE julia DOC "Julia executable" )
if(Julia_EXECUTABLE)
  execute_process(
    COMMAND "${Julia_EXECUTABLE}" --startup-file=no --version
    OUTPUT_VARIABLE Julia_VERSION_STRING
    )
  string(
    REGEX REPLACE ".*([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1"
    Julia_VERSION_STRING "${Julia_VERSION_STRING}"
    )

  message( STATUS "Julia Version: ${Julia_VERSION_STRING}")

  if ( ${Julia_VERSION_STRING} VERSION_LESS "1.4.0" )
    message( STATUS "Julia CTests: Disabled" )
    set( Julia_ENABLE_TESTS OFF )
  else()
    message( STATUS "Julia CTests: Enabled" )
    set( Julia_ENABLE_TESTS ON )
  endif()
endif()

if ( Julia_ENABLE_TESTS )

  set( Julia_DEPOT_PATH ${CMAKE_CURRENT_BINARY_DIR}/depot )
  set( Julia_MPI_BINARY "system" )
  set( Julia_ENVIRONEMENT "JULIA_DEPOT_PATH=${Julia_DEPOT_PATH}" "JULIA_MPI_BINARY=${Julia_MPI_BINARY}" )

  #
  # Make sure the PaStiX package is built before use, and if necessary
  # with the correct MPI before execution of the actual testings
  #
  if ( PASTIX_WITH_MPI )
    set( _julia_init_ ${Julia_EXECUTABLE} "--project" "-e" "using Pkg\; Pkg.activate(\"PaStiX\")\; Pkg.develop(PackageSpec(path=\"../../spm/wrappers/julia/spm\"))\; Pkg.instantiate()\; Pkg.build(\"MPI\")" )
  else()
    set( _julia_init_ ${Julia_EXECUTABLE} "--project" "-e" "using Pkg\; Pkg.activate(\"PaStiX\")\; Pkg.develop(PackageSpec(path=\"../../spm/wrappers/julia/spm\"))\; Pkg.instantiate()" )
  endif()

  add_test( shm_julia_set_mpi  ${_julia_init_} )
  set_tests_properties( shm_julia_set_mpi
    PROPERTIES ENVIRONMENT "${Julia_ENVIRONEMENT}" )
  if ( PASTIX_WITH_MPI )
    add_test( mpi_julia_set_mpi ${_julia_init_} )
    set_tests_properties( mpi_julia_set_mpi
      PROPERTIES ENVIRONMENT "${Julia_ENVIRONEMENT}" )
  endif()

  set( JULIA_TESTS
    simple step-by-step )

  foreach( example ${JULIA_TESTS} )
    set( _test_suffix_ julia_${example} )
    set( _test_file_   ${CMAKE_CURRENT_SOURCE_DIR}/examples/${example}.jl )

    set( _test_name_ shm_${_test_suffix_} )
    add_test( ${_test_name_} ${Julia_EXECUTABLE} ${_test_file_} )
    set_tests_properties( ${_test_name_}
      PROPERTIES ENVIRONMENT "${Julia_ENVIRONEMENT}" )

    if (PASTIX_WITH_MPI)
      set( _test_name_ mpi_${_test_suffix_} )
      add_test( ${_test_name_} mpiexec -np 4 --host localhost:4 ${Julia_EXECUTABLE} ${_test_file_} )
      set_tests_properties( ${_test_name_}
        PROPERTIES ENVIRONMENT "${Julia_ENVIRONEMENT}" )
    endif()
  endforeach()

  set( example schur )
  set( _test_suffix_ julia_${example} )
  set( _test_file_   ${CMAKE_CURRENT_SOURCE_DIR}/examples/${example}.jl )

  set( _test_name_ shm_${_test_suffix_} )
  add_test( ${_test_name_} ${Julia_EXECUTABLE} ${_test_file_} )
  set_tests_properties( ${_test_name_}
    PROPERTIES ENVIRONMENT "${Julia_ENVIRONEMENT}" )
endif()
