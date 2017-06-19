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

# Configure enum.py
if (PASTIX_INT64)
  set(PASTIX_PYTHON_INTEGER c_int64)
else()
  set(PASTIX_PYTHON_INTEGER c_int)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/python/examples/pypastix/enum.py.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/python/examples/pypastix/enum.py" @ONLY)

# Install python package
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/pypastix/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/pypastix/spm.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/pypastix/pastix.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/pypastix/enum.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/pypastix/solver.py
  DESTINATION lib/python/pypastix )

# Install python examples
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/simple.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/simple_obj.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/schur.py
  ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/schur_obj.py
  DESTINATION examples
  )

## CTest execution
find_package(PythonInterp QUIET)
message(STATUS "TOTO")
if (PYTHONINTERP_FOUND)
  set( PYTHON_TESTS
    simple schur simple_obj schur_obj)

  foreach(example ${PYTHON_TESTS} )
    message(STATUS "python_${example} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/${example}.py")
    add_test(python_${example} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/examples/${example}.py)
  endforeach()
endif()

