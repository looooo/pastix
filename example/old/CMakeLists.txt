###
#
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.1
#  @author Mathieu Faverge
#  @author Alycia Lisito
#  @date 2023-10-25
#
###

set (EXAMPLES
  simple.c
  step-by-step.c
  )

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable("old_${_name_we}" ${_file})
  target_link_libraries("old_${_name_we}" pastix)

  install(TARGETS old_${_name_we} RUNTIME DESTINATION share/doc/pastix/examples/old )
  install(FILES   ${_file}                DESTINATION share/doc/pastix/examples/old )
endforeach()

### Install a makefile for user compilation test
set( EXAMPLE_LIBS ${DEPS_LIBRARIES} )
if ( EXAMPLE_LIBS )
  gpc_convert_libstyle_to_pkgconfig( EXAMPLE_LIBS )
  string(REPLACE ";" " " EXAMPLE_LIBS "${EXAMPLE_LIBS}")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Makefile.in"
  "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/make/Makefile"
  DESTINATION share/doc/pastix/examples/old )

