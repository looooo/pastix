# - Find SCOTCH library
# This module finds an installed  library that implements the SCOTCH
# linear-algebra interface (see http://icl.cs.utk.edu/scotch/).
# The list of libraries searched for is taken
# from the autoconf macro file, acx_blas.m4 (distributed at
# http://ac-archive.sourceforge.net/ac-archive/acx_blas.html).
#
# This module sets the following variables:
#  SCOTCH_FOUND        - set to true if a library implementing the SCOTCH
#    interface is found
#  SCOTCH_LIBRARIES    - only the libraries (w/o the '-l')
#  SCOTCH_LIBRARY_DIRS - the paths of the libraries (w/o the '-L')
#  SCOTCH_LDFLAGS      - all required linker flags
#  SCOTCH_INCLUDE_DIRS - the '-I' preprocessor flags (w/o the '-I')
#  SCOTCH_CFLAGS       - all required cflags
#
#  PTSCOTCH_FOUND        - set to true if a library implementing the SCOTCH
#    interface is found
#  PTSCOTCH_LIBRARIES    - only the libraries (w/o the '-l')
#  PTSCOTCH_LIBRARY_DIRS - the paths of the libraries (w/o the '-L')
#  PTSCOTCH_LDFLAGS      - all required linker flags
#  PTSCOTCH_INCLUDE_DIRS - the '-I' preprocessor flags (w/o the '-I')
#  PTSCOTCH_CFLAGS       - all required cflags
#
##########

find_library(SCOTCH_LIBRARY scotch DOC "Scotch library")
find_library(PTSCOTCH_LIBRARY ptscotch DOC "PT-Scotch library")
mark_as_advanced(SCOTCH_LIBRARY PTSCOTCH_LIBRARY)

set(SCOTCH_DIR "" CACHE PATH "Root directory containing SCOTCH")

find_package(PkgConfig QUIET)
if( SCOTCH_DIR )
  set(ENV{PKG_CONFIG_PATH} "${SCOTCH_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()
pkg_check_modules(PC_SCOTCH QUIET scotch)
set(SCOTCH_DEFINITIONS ${PC_SCOTCH_CFLAGS_OTHER} )

find_path(SCOTCH_INCLUDE_DIR scotch.h
	  HINTS ${SCOTCH_DIR} ${PC_SCOTCH_INCLUDEDIR} ${PC_SCOTCH_INCLUDE_DIRS}
	  PATH_SUFFIXES include
	  DOC "SCOTCH includes" )
set(SCOTCH_INCLUDE_DIRS ${SCOTCH_INCLUDE_DIR})

find_library(SCOTCH_LIBRARY scotch
	     HINTS ${SCOTCH_DIR} ${PC_SCOTCH_LIBDIR} ${PC_SCOTCH_LIBRARY_DIRS}
	     PATH_SUFFIXES lib
	     DOC "Where the SCOTCH libraries are")
set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCOTCH
    "Could NOT find SCOTCH; Options depending on SCOTCH will be disabled
  if needed, please specify the library location
    - using SCOTCH_DIR [${SCOTCH_DIR}]
    - or a combination of SCOTCH_INCLUDE_DIR [${SCOTCH_INCLUDE_DIR}] and SCOTCH_LIBRARY [${SCOTCH_LIBRARY}]"
    SCOTCH_LIBRARY SCOTCH_INCLUDE_DIR )

if(SCOTCH_FOUND)
  set(SCOTCH_SAVE_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
  list(APPEND CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIR})
  check_struct_has_member( "struct scotch_obj" parent scotch.h HAVE_SCOTCH_PARENT_MEMBER )
  check_struct_has_member( "struct scotch_cache_attr_s" size scotch.h HAVE_SCOTCH_CACHE_ATTR )
  check_c_source_compiles( "#include <scotch.h>
    int main(void) { scotch_obj_t o; o->type = SCOTCH_OBJ_PU; return 0;}" HAVE_SCOTCH_OBJ_PU)
  check_library_exists(${SCOTCH_LIBRARY} scotch_bitmap_free "" HAVE_SCOTCH_BITMAP)
  set(CMAKE_REQUIRED_INCLUDES ${SCOTCH_SAVE_CMAKE_REQUIRED_INCLUDES})
endif()


unset(SCOTCH_C_COMPILE_SUCCESS)
unset(SCOTCH_F_COMPILE_SUCCESS)

# First we try to use pkg-config to find what we're looking for
# in the directory specified by the SCOTCH_DIR or SCOTCH_PKG_DIR
if(SCOTCH_DIR)
  if(NOT SCOTCH_PKG_DIR)
    set(SCOTCH_PKG_DIR "${SCOTCH_DIR}/lib/pkgconfig")
  endif(NOT SCOTCH_PKG_DIR)
endif(SCOTCH_DIR)

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  set(ENV{PKG_CONFIG_PATH} "${SCOTCH_PKG_DIR}:$ENV{PKG_CONFIG_PATH}")
  pkg_check_modules(SCOTCH scotch)
endif(PKG_CONFIG_FOUND)

if(NOT SCOTCH_FOUND)
  #
  # No pkg-config supported on this system. Hope the user provided
  # all the required variables in order to detect SCOTCH. This include
  # either:
  # - SCOTCH_DIR: root dir to SCOTCH installation
  # - SCOTCH_PKG_DIR: directory where the scotch.pc is installed
  #
  if(NOT SCOTCH_DIR AND NOT SCOTCH_PKG_DIR)
    if(SCOTCH_FIND_REQUIRED)
      message(FATAL_ERROR "pkg-config not available. You need to provide SCOTCH_DIR or SCOTCH_PKG_DIR")
    else(SCOTCH_FIND_REQUIRED)
      message(STATUS "pkg-config not available. You need to provide SCOTCH_DIR or SCOTCH_PKG_DIR")
    endif(SCOTCH_FIND_REQUIRED)

  else(NOT SCOTCH_DIR AND NOT SCOTCH_PKG_DIR)

    if (SCOTCH_PKG_DIR)
      set(_scotch_pkg_file "${SCOTCH_PKG_DIR}/scotch.pc")
    else()
      set(_scotch_pkg_file "${SCOTCH_DIR}/lib/pkgconfig/scotch.pc")
    endif()

    if (NOT EXISTS ${_scotch_pkg_file})
      if(SCOTCH_FIND_REQUIRED)
        message(FATAL_ERROR "${_scotch_pkg_file} doesn't exist")
      else(SCOTCH_FIND_REQUIRED)
        message(STATUS "${_scotch_pkg_file} doesn't exist")
      endif(SCOTCH_FIND_REQUIRED)
    else(NOT EXISTS ${_scotch_pkg_file})
      file(STRINGS "${_scotch_pkg_file}" _cflags REGEX "Cflags:")
      file(STRINGS "${_scotch_pkg_file}" _libs   REGEX "Libs:")

      string(REGEX REPLACE "Cflags:" "" _cflags ${_cflags})
      string(REGEX REPLACE "Libs:"   "" _libs   ${_libs}  )
      string(REGEX REPLACE " +" ";" _cflags ${_cflags})
      string(REGEX REPLACE " +" ";" _libs   ${_libs}  )

      foreach(_cflag ${_cflags})
        string(REGEX REPLACE "^-I(.*)" "\\1" _incdir "${_cflag}")
        if ("${_cflag}" MATCHES "-I.*")
          list(APPEND SCOTCH_INCLUDE_DIRS ${_incdir})
        endif()
        list(APPEND SCOTCH_CFLAGS ${_cflag})
      endforeach()

      foreach(_lib ${_libs})
        string(REGEX REPLACE "^-L(.*)" "\\1" _libdir "${_lib}")
        if ("${_lib}" MATCHES "-L.*")
          list(APPEND SCOTCH_LIBRARY_DIRS ${_libdir})
        endif()

        string(REGEX REPLACE "^-l(.*)" "\\1" _onelib "${_lib}")
        if ("${_lib}" MATCHES "-l.*")
          list(APPEND SCOTCH_LIBRARIES ${_onelib})
        endif()

        list(APPEND SCOTCH_LDFLAGS ${_lib})
      endforeach()

      list(REMOVE_DUPLICATES SCOTCH_INCLUDE_DIRS)
      list(REMOVE_DUPLICATES SCOTCH_CFLAGS)
      list(REMOVE_DUPLICATES SCOTCH_LIBRARY_DIRS)
      list(REMOVE_DUPLICATES SCOTCH_LIBRARIES)
      list(REMOVE_DUPLICATES SCOTCH_LDFLAGS)

    endif(NOT EXISTS ${_scotch_pkg_file})
  endif(NOT SCOTCH_DIR AND NOT SCOTCH_PKG_DIR)
  set(SCOTCH_FOUND TRUE)
endif(NOT SCOTCH_FOUND)

if(SCOTCH_FOUND)

  #
  # There is a circular dependency in SCOTCH between the libscotch and libcoreblas.
  # Unfortunately, this cannot be handled by pkg-config (as it remove the duplicates)
  # so we have to add it by hand.
  # Those parameters are also removed by pkg-config if they are present around mkl libs
  #
  if(HAVE_LINKER_GROUP)
    list(INSERT SCOTCH_LIBRARIES 0 -Wl,--start-group)
    list(APPEND SCOTCH_LIBRARIES -Wl,--end-group)
  else()
    list(APPEND SCOTCH_LIBRARIES scotch)
  endif()

  # Validate the include file <scotch.h>
  include(CheckIncludeFile)

  set(SCOTCH_tmp_includes ${CMAKE_REQUIRED_INCLUDES})
  list(APPEND CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIRS})

  check_include_file(scotch.h SCOTCH_SCOTCH_H_FOUND)

  if ( NOT SCOTCH_SCOTCH_H_FOUND )
    if(SCOTCH_FIND_REQUIRED)
      message(FATAL_ERROR "Couln't find the scotch.h header in ${SCOTCH_INCLUDE_DIRS}")
    endif(SCOTCH_FIND_REQUIRED)
    set(SCOTCH_FOUND FALSE)
    return()
  endif()

  # Validate the library
  include(CheckCSourceCompiles)

  set(SCOTCH_tmp_libraries ${CMAKE_REQUIRED_LIBRARIES})
  set(SCOTCH_tmp_flags     ${CMAKE_REQUIRED_FLAGS})
  list(APPEND CMAKE_REQUIRED_LIBRARIES ${SCOTCH_LIBRARIES})

  # CMAKE_REQUIRED_FLAGS must be a string, not a list
  # if CMAKE_REQUIRED_FLAGS is a list (separated by ;), only the first element of the list is passed to check_c_source_compile
  # Since SCOTCH_LDFLAGS and SCOTCH_CFLAGS hold lists, we convert them by hand to a string
  foreach(arg ${SCOTCH_LDFLAGS})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${arg}")
  endforeach(arg ${SCOTCH_LDFLAGS})
  foreach(arg ${SCOTCH_CFLAGS})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${arg}")
  endforeach(arg ${SCOTCH_CFLAGS})
  foreach(arg ${SCOTCH_LIBRARY_DIRS})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${arg}")
  endforeach(arg ${SCOTCH_CFLAGS})

  check_c_source_compiles(
    "int main(int argc, char* argv[]) {
       SCOTCH_dgeqrf(); return 0;
     }"
    SCOTCH_C_COMPILE_SUCCESS
    )

  if(NOT SCOTCH_C_COMPILE_SUCCESS)
    get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
    if(NOT _LANGUAGES_ MATCHES Fortran)
      if(SCOTCH_FIND_REQUIRED)
        message(FATAL_ERROR "Find SCOTCH requires Fortran support so Fortran must be enabled.")
      else(SCOTCH_FIND_REQUIRED)
        message(STATUS "Looking for SCOTCH... - NOT found (Fortran not enabled)") #
        set(SCOTCH_FOUND FALSE)
        return()
      endif(SCOTCH_FIND_REQUIRED)
    endif(NOT _LANGUAGES_ MATCHES Fortran)
    include(CheckFortranFunctionExists)
    list(APPEND CMAKE_REQUIRED_LIBRARIES ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
    check_c_source_compiles(
      "int main(int argc, char* argv[]) {
       SCOTCH_dgeqrf(); return 0;
     }"
      SCOTCH_F_COMPILE_SUCCESS
      )

    if(NOT SCOTCH_F_COMPILE_SUCCESS)
      if(SCOTCH_FIND_REQUIRED)
        message(FATAL_ERROR "Find SCOTCH requires Fortran support so Fortran must be enabled.")
      else(SCOTCH_FIND_REQUIRED)
        message(STATUS "Looking for SCOTCH... - NOT found")
        set(SCOTCH_FOUND FALSE)
        return()
      endif(SCOTCH_FIND_REQUIRED)
    endif(NOT SCOTCH_F_COMPILE_SUCCESS)
  endif(NOT SCOTCH_C_COMPILE_SUCCESS)

  set(${CMAKE_REQUIRED_INCLUDES}  SCOTCH_tmp_includes)
  set(${CMAKE_REQUIRED_LIBRARIES} SCOTCH_tmp_libraries)
  set(${CMAKE_REQUIRED_FLAGS}     SCOTCH_tmp_flags)
  unset(SCOTCH_tmp_libraries)
  unset(SCOTCH_tmp_includes)
  unset(SCOTCH_tmp_flags)
endif(SCOTCH_FOUND)

if(NOT SCOTCH_FIND_QUIETLY)
  set(SCOTCH_status_message
    "
    SCOTCH_CFLAGS       = [${SCOTCH_CFLAGS}]
    SCOTCH_LDFLAGS      = [${SCOTCH_LDFLAGS}]
    SCOTCH_INCLUDE_DIRS = [${SCOTCH_INCLUDE_DIRS}]
    SCOTCH_LIBRARY_DIRS = [${SCOTCH_LIBRARY_DIRS}]
    SCOTCH_LIBRARIES = [${SCOTCH_LIBRARIES}]")

  if(SCOTCH_C_COMPILE_SUCCESS OR SCOTCH_F_COMPILE_SUCCESS)
    if(SCOTCH_F_COMPILE_SUCCESS)
      set(SCOTCH_REQUIRE_FORTRAN_LINKER TRUE)
      mark_as_advanced(SCOTCH_REQUIRE_FORTRAN_LINKER)
      message(STATUS "A Library with SCOTCH API found (using C compiler and Fortran linker).")
    endif(SCOTCH_F_COMPILE_SUCCESS)
    string(REGEX REPLACE ";" " " SCOTCH_LDFLAGS "${SCOTCH_LDFLAGS}")
    set(SCOTCH_FOUND TRUE)
    find_package_message(SCOTCH
      "Found SCOTCH: ${SCOTCH_status_message}"
      "[${SCOTCH_CFLAGS}][${SCOTCH_LDFLAGS}][${SCOTCH_INCLUDE_DIRS}][${SCOTCH_LIBRARY_DIRS}]")
  else(SCOTCH_C_COMPILE_SUCCESS OR SCOTCH_F_COMPILE_SUCCESS)
    if(SCOTCH_FIND_REQUIRED)
      message(FATAL_ERROR
        "A required library with SCOTCH API not found. Please specify library location.${SCOTCH_status_message}")
    else(SCOTCH_FIND_REQUIRED)
      message(STATUS
        "A library with SCOTCH API not found. Please specify library location.${SCOTCH_status_message}")
    endif(SCOTCH_FIND_REQUIRED)
  endif(SCOTCH_C_COMPILE_SUCCESS OR SCOTCH_F_COMPILE_SUCCESS)
endif(NOT SCOTCH_FIND_QUIETLY)

mark_as_advanced(SCOTCH_PKG_DIR SCOTCH_LIBRARIES SCOTCH_INCLUDE_DIRS SCOTCH_LINKER_FLAGS)
set(SCOTCH_DIR          "${SCOTCH_DIR}"          CACHE PATH   "Location of the SCOTCH library" FORCE)
set(SCOTCH_PKG_DIR      "${SCOTCH_PKG_DIR}"      CACHE PATH   "Location of the SCOTCH pkg-config decription file" FORCE)
#set(SCOTCH_LIBRARIES "${SCOTCH_LIBRARIES}" CACHE STRING "libraries to link with SCOTCH" FORCE)

unset(SCOTCH_C_COMPILE_SUCCESS)
unset(SCOTCH_F_COMPILE_SUCCESS)


# handle the QUIETLY and REQUIRED arguments and set SCOTCH_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCOTCH REQUIRED_VARS SCOTCH_EXECUTABLE
                                   VERSION_VAR SCOTCH_VERSION)

