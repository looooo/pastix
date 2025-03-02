###
#
#  @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.4.0
#  @author Mathieu Faverge
#  @author Tony Delarue
#  @author Alycia Lisito
#  @author Florent Pruvost
#  @date 2024-07-05
#
###
cmake_minimum_required(VERSION 3.12)

find_package (Python3 COMPONENTS Interpreter QUIET)

if(NOT Python3_FOUND)
  return()
endif()

# Configure enum.py
if (PASTIX_INT64)
  set(PASTIX_PYTHON_INTEGER c_int64)
else()
  set(PASTIX_PYTHON_INTEGER c_int)
endif()

if (PASTIX_WITH_MPI)
  set(PASTIX_PYTHON_MPI_ENABLED 1)
else()
  set(PASTIX_PYTHON_MPI_ENABLED 0)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/examples/pypastix/enum.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/examples/pypastix/enum.py" @ONLY)

# Copy wrapper to build
file(COPY
  examples/pypastix/__init__.py
  examples/pypastix/__pastix__.py
  examples/pypastix/pastix.py
  examples/pypastix/solver.py
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/examples/pypastix )

# Copy examples to build
file(COPY
  examples/schur_obj.py
  examples/schur.py
  examples/simple_obj.py
  examples/simple.py
  examples/step-by-step.py
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/examples )

#
# Let's define the installation directory by building something similar to what
# sysconfig.get_path('platlib') returns
#
if ( NOT PYTHON_SITE_PACKAGES_DIR )
  set(PYTHON_SITE_PACKAGES_DIR lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages CACHE STRING "Python packages installation path" FORCE)
endif()

# Install python package
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/examples/pypastix/__init__.py
  ${CMAKE_CURRENT_BINARY_DIR}/examples/pypastix/__pastix__.py
  ${CMAKE_CURRENT_BINARY_DIR}/examples/pypastix/pastix.py
  ${CMAKE_CURRENT_BINARY_DIR}/examples/pypastix/enum.py
  ${CMAKE_CURRENT_BINARY_DIR}/examples/pypastix/solver.py
  DESTINATION ${PYTHON_SITE_PACKAGES_DIR}/pastix )

# Install python examples
install(FILES
  examples/simple.py
  examples/simple_obj.py
  examples/schur.py
  examples/schur_obj.py
  examples/step-by-step.py
  DESTINATION share/doc/pastix/examples/python)

# Documentation
# -------------
add_documented_files(
  examples/simple.py
  examples/simple_obj.py
  examples/schur.py
  examples/schur_obj.py
  examples/step-by-step.py
)

## CTest execution
if (Python3_Interpreter_FOUND)
  set( PYTHON_TESTS
    simple step-by-step simple_obj )

  # Add path to SPM, and to the .so in the environment to run the tests
  if( NOT SPM_FOUND )
    set( TEST_ENV_LIST
      "PYTHONPATH=$ENV{PYTHONPATH}:${CMAKE_BINARY_DIR}/spm/wrappers/python"
      "LIBRARY_PATH=$ENV{LIBRARY_PATH}:${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/spm/src"
      "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/spm/src"
      "DYLD_LIBRARY_PATH=$ENV{DYLD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/spm/src"
      )
  endif()
  foreach(example ${PYTHON_TESTS} )
    set( _test_suffix_ ${example} )
    set( _test_file_   ${CMAKE_CURRENT_BINARY_DIR}/examples/${example}.py )

    set( _test_name_ python_shm_${_test_suffix_} )
    add_test( ${_test_name_} ${Python3_EXECUTABLE} ${_test_file_} )
    set_tests_properties( ${_test_name_}
      PROPERTIES ENVIRONMENT "${TEST_ENV_LIST}" )

    if (PASTIX_WITH_MPI)
      set( _test_name_ python_mpi_${_test_suffix_} )
      add_test( ${_test_name_} ${pastix_mpiexec} ${Python3_EXECUTABLE} ${_test_file_} )
      set_tests_properties( ${_test_name_}
        PROPERTIES ENVIRONMENT "${TEST_ENV_LIST}" )
    endif()
  endforeach()

  #
  # MPI is not yeat available with the Schur complement
  #
  set( PYTHON_TESTS
    schur schur_obj )

  foreach(example ${PYTHON_TESTS} )
    set( _test_suffix_ ${example} )
    set( _test_file_   ${CMAKE_CURRENT_BINARY_DIR}/examples/${example}.py )

    set( _test_name_ python_shm_${_test_suffix_} )
    add_test( ${_test_name_} ${Python3_EXECUTABLE} ${_test_file_} )
    set_tests_properties( ${_test_name_}
      PROPERTIES ENVIRONMENT "${TEST_ENV_LIST}" )

  endforeach()
endif()
