###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find TMG include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(TMG
#               [REQUIRED]             # Fail with error if tmg is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )
# This module finds headers and tmg library.
# Results are reported in variables:
#  TMG_FOUND           - True if headers and requested libraries were found
#  TMG_INCLUDE_DIRS    - tmg include directories
#  TMG_LIBRARY_DIRS    - Link directories for tmg libraries
#  TMG_LIBRARIES       - tmg component libraries to be linked
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DTMG=path/to/tmg):
#  TMG_DIR             - Where to find the base directory of tmg
#  TMG_INCDIR          - Where to find the header files
#  TMG_LIBDIR          - Where to find the library files

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)


if (NOT TMG_FOUND)
    set(TMG_DIR "" CACHE PATH "Root directory of TMG library")
    if (NOT TMG_FIND_QUIETLY)
        message(STATUS "A cache variable, namely TMG_DIR, has been set to specify the install directory of TMG")
    endif()
endif()


# used to test a TMG function after
get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT _LANGUAGES_ MATCHES Fortran)
    include(CheckFunctionExists)
else (NOT _LANGUAGES_ MATCHES Fortran)
    include(CheckFortranFunctionExists)
endif (NOT _LANGUAGES_ MATCHES Fortran)

# TMG depends on LAPACK
# try to find it specified as COMPONENTS during the call
if (TMG_FIND_COMPONENTS)
    foreach( component ${TMG_FIND_COMPONENTS} )
        if(TMG_FIND_REQUIRED_${component})
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(TMG_${component}_FOUND TRUE)
        else()
            set(TMG_${component}_FOUND FALSE)
        endif()
    endforeach()
endif ()

# TMG depends on LAPACK anyway, try to find it
if (NOT LAPACK_FOUND)
    if(TMG_FIND_REQUIRED)
        find_package(LAPACK REQUIRED)
    else()
        find_package(LAPACK)
    endif()
endif()

# TMG depends on LAPACK
if (LAPACK_FOUND)

    # check if a tmg function exists in the LAPACK lib
    set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LIBRARIES}")
    include(CheckFunctionExists)
    include(CheckFortranFunctionExists)
    unset(TMG_WORKS CACHE)
    if (NOT _LANGUAGES_ MATCHES Fortran)
        check_function_exists(dlarnv TMG_WORKS)
    else (NOT _LANGUAGES_ MATCHES Fortran)
        check_fortran_function_exists(dlarnv TMG_WORKS)
    endif (NOT _LANGUAGES_ MATCHES Fortran)
    if (TMG_WORKS)
        unset(TMG_WORKS CACHE)
        if (NOT _LANGUAGES_ MATCHES Fortran)
            check_function_exists(dlagsy TMG_WORKS)
        else (NOT _LANGUAGES_ MATCHES Fortran)
            check_fortran_function_exists(dlagsy TMG_WORKS)
        endif (NOT _LANGUAGES_ MATCHES Fortran)
        mark_as_advanced(TMG_WORKS)
    endif()
    set(CMAKE_REQUIRED_LIBRARIES)

    if(TMG_WORKS)
        if(NOT TMG_FIND_QUIETLY)
            message(STATUS "Looking for tmg: test with lapack succeeds")
        endif()
        # test succeeds: TMG is in LAPACK
        set(TMG_LIBRARIES "${LAPACK_LIBRARIES}")
        set(TMG_LIBRARY_DIRS "${LAPACK_LIBRARY_DIRS}")
        if(LAPACK_INCLUDE_DIRS)
            set(TMG_INCLUDE_DIRS "${LAPACK_INCLUDE_DIRS}")
        endif()
    else()

        if(NOT TMG_FIND_QUIETLY)
            message(STATUS "Looking for tmg : test with lapack fails")
            message(STATUS "Looking for tmg : try to find it elsewhere")
        endif()
        # test fails: try to find TMG lib exterior to LAPACK

        # Looking for lib tmg
        # -------------------

        # Add system library paths to search lib
        # --------------------------------------
        unset(_lib_env)
        if(WIN32)
            string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
        else()
            if(APPLE)
                string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
            else()
                string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
            endif()
            list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
            list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
        endif()
        list(REMOVE_DUPLICATES _lib_env)

        # Try to find the tmg lib in the given paths
        # ----------------------------------------------

        # call cmake macro to find the lib path
        if(TMG_LIBDIR)
            set(TMG_tmg_LIBRARY "TMG_tmg_LIBRARY-NOTFOUND")
            find_library(TMG_tmg_LIBRARY
                NAMES tmg
                HINTS ${TMG_LIBDIR} )
        else()
            if(TMG_DIR)
                set(TMG_tmg_LIBRARY "TMG_tmg_LIBRARY-NOTFOUND")
                find_library(TMG_tmg_LIBRARY
                    NAMES tmg
                    HINTS ${TMG_DIR}
                    PATH_SUFFIXES lib lib32 lib64 )
            else()
                set(TMG_tmg_LIBRARY "TMG_tmg_LIBRARY-NOTFOUND")
                find_library(TMG_tmg_LIBRARY
                    NAMES tmg
                    HINTS ${_lib_env} )
            endif()
        endif()
        mark_as_advanced(TMG_tmg_LIBRARY)

        # If found, add path to cmake variable
        # ------------------------------------
        if (TMG_tmg_LIBRARY)
            get_filename_component(tmg_lib_path ${TMG_tmg_LIBRARY} PATH)
            # set cmake variables (respects naming convention)
            set(TMG_LIBRARIES    "${TMG_tmg_LIBRARY}")
            set(TMG_LIBRARY_DIRS "${tmg_lib_path}")
        else ()
            set(TMG_LIBRARIES    "TMG_LIBRARIES-NOTFOUND")
            set(TMG_LIBRARY_DIRS "TMG_LIBRARY_DIRS-NOTFOUND")
            if(NOT TMG_FIND_QUIETLY)
                message(STATUS "Looking for tmg -- lib tmg not found")
            endif()
        endif ()

        if (TMG_LIBRARY_DIRS)
            list(REMOVE_DUPLICATES TMG_LIBRARY_DIRS)
        endif ()

        # check a function to validate the find
        if (TMG_INCLUDE_DIRS)
            set(CMAKE_REQUIRED_INCLUDES  "${TMG_INCLUDE_DIRS}")
        endif()
        if (LAPACK_INCLUDE_DIRS)
            list(APPEND CMAKE_REQUIRED_INCLUDES "${LAPACK_INCLUDE_DIRS}")
        endif()
        set(CMAKE_REQUIRED_LIBRARIES "${TMG_LIBRARIES};${LAPACK_LIBRARIES}")
        if (TMG_LIBRARY_DIRS)
            set(CMAKE_REQUIRED_FLAGS "-L${TMG_LIBRARY_DIRS}")
        endif()
        if (LAPACK_LIBRARY_DIRS)
            set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${LAPACK_LIBRARY_DIRS}")
        endif()

        if(TMG_LIBRARIES)
            unset(TMG_WORKS CACHE)
            include(CheckFunctionExists)
            include(CheckFortranFunctionExists)
            if (NOT _LANGUAGES_ MATCHES Fortran)
                check_function_exists(dlarnv TMG_WORKS)
            else (NOT _LANGUAGES_ MATCHES Fortran)
                check_fortran_function_exists(dlarnv TMG_WORKS)
            endif (NOT _LANGUAGES_ MATCHES Fortran)
            if (TMG_WORKS)
                unset(TMG_WORKS CACHE)
                if (NOT _LANGUAGES_ MATCHES Fortran)
                    check_function_exists(dlagsy TMG_WORKS)
                else (NOT _LANGUAGES_ MATCHES Fortran)
                    check_fortran_function_exists(dlagsy TMG_WORKS)
                endif (NOT _LANGUAGES_ MATCHES Fortran)
                mark_as_advanced(TMG_WORKS)
            endif()

            if(TMG_WORKS)
                set(TMG_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
                set(TMG_LIBRARY_DIRS "${TMG_LIBRARY_DIRS}" "${LAPACK_LIBRARY_DIRS}")
                set(TMG_INCLUDE_DIRS ${CMAKE_REQUIRED_INCLUDES})
            else()
                if(NOT TMG_FIND_QUIETLY)
                    message(STATUS "Looking for tmg: test of dlarnv and dlagsy with tmg and lapack libraries fails")
                    message(STATUS "TMG_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
                    message(STATUS "TMG_LIBRARY_DIRS: ${CMAKE_REQUIRED_FLAGS}")
                    message(STATUS "TMG_INCLUDE_DIRS: ${CMAKE_REQUIRED_INCLUDES}")
                    message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
                    message(STATUS "Looking for tmg : set TMG_LIBRARIES to NOTFOUND")
                endif()
                set(TMG_LIBRARIES "TMG_LIBRARIES-NOTFOUND")
            endif()
            set(CMAKE_REQUIRED_INCLUDES)
            set(CMAKE_REQUIRED_FLAGS)
            set(CMAKE_REQUIRED_LIBRARIES)
        endif(TMG_LIBRARIES)

    endif()

else()

    if(NOT TMG_FIND_QUIETLY)
        message(STATUS "TMG requires LAPACK but LAPACK has not been found."
            "Please look for LAPACK first.")
    endif()

endif()


# check that TMG has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TMG DEFAULT_MSG
                                  TMG_LIBRARIES)
