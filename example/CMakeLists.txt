###
#
#  @copyright 2013-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @date 2013-06-24
#
###

add_subdirectory(old)
set (EXAMPLES
  analyze.c
  bench_facto.c
  compress.c
  personal.c
  reentrant.c
  schur.c
  simple.c
  step-by-step.c
  )

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix)

  install(TARGETS ${_name_we} RUNTIME DESTINATION examples )
  install(FILES   ${_file}            DESTINATION examples )
endforeach()

### Install a makefile for user compilation test
set( EXAMPLE_LIBS
  ${DEPS_LIBRARIES}
  )
list(REMOVE_DUPLICATES EXAMPLE_LIBS)
string(REPLACE ";" " " EXAMPLE_LIBS "${EXAMPLE_LIBS}")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Makefile.in"
  "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile"
  DESTINATION examples)

include(CTestLists.txt)
