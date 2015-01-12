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
# - Find SCOTCH include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(SCOTCH
#               [REQUIRED]             # Fail with error if scotch is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )  
# This module finds headers and scotch library. 
# Results are reported in variables:
#  SCOTCH_FOUND           - True if headers and requested libraries were found
#  SCOTCH_INCLUDE_DIRS    - scotch include directories
#  SCOTCH_LIBRARY_DIRS    - Link directories for scotch libraries
#  SCOTCH_LIBRARIES       - scotch component libraries to be linked
# The user can give specific paths where to find the libraries adding cmake 
# options at configure (ex: cmake path/to/project -DSCOTCH=path/to/scotch):
#  SCOTCH_DIR             - Where to find the base directory of scotch
#  SCOTCH_INCDIR          - Where to find the header files
#  SCOTCH_LIBDIR          - Where to find the library files

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


# Some macros to print status when search for headers and libs
# PrintFindStatus.cmake is in cmake_modules/morse/find directory
include(PrintFindStatus)

# SCOTCH may depend on Threads
# try to find it specified as COMPONENTS during the call
if( SCOTCH_FIND_COMPONENTS )
    foreach( component ${SCOTCH_FIND_COMPONENTS} )
        if(${SCOTCH_FIND_REQUIRED_${component}} STREQUAL 1)
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(SCOTCH_${component}_FOUND TRUE)
            if( THREADS_FOUND )
                list(APPEND EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
            endif ()
        else()
            set(SCOTCH_${component}_FOUND FALSE)
        endif()
    endforeach()
endif()


# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
if(WIN32)
    string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
else()
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
endif()
list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
list(REMOVE_DUPLICATES _inc_env)


# Try to find the scotch header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(SCOTCH_INCDIR)
    set(SCOTCH_scotch.h_DIRS "SCOTCH_scotch.h_DIRS-NOTFOUND")
    find_path(SCOTCH_scotch.h_DIRS
      NAMES scotch.h
      HINTS ${SCOTCH_INCDIR})
else()
    if(SCOTCH_DIR)
        set(SCOTCH_scotch.h_DIRS "SCOTCH_scotch.h_DIRS-NOTFOUND")
        find_path(SCOTCH_scotch.h_DIRS
          NAMES scotch.h
          HINTS ${SCOTCH_DIR}
          PATH_SUFFIXES include)        
    else()
        set(SCOTCH_scotch.h_DIRS "SCOTCH_scotch.h_DIRS-NOTFOUND")
        find_path(SCOTCH_scotch.h_DIRS
          NAMES scotch.h
          HINTS ${_inc_env})
    endif()
endif()
mark_as_advanced(SCOTCH_scotch.h_DIRS)

# Print status if not found
# -------------------------
if (NOT SCOTCH_scotch.h_DIRS AND NOT SCOTCH_FIND_QUIETLY)
    Print_Find_Header_Status(scotch scotch.h)
endif ()

# If found, add path to cmake variable
# ------------------------------------
if (SCOTCH_scotch.h_DIRS)
    set(SCOTCH_INCLUDE_DIRS "${SCOTCH_scotch.h_DIRS}")
else ()
    set(SCOTCH_INCLUDE_DIRS "SCOTCH_INCLUDE_DIRS-NOTFOUND")
    if (NOT SCOTCH_FIND_QUIETLY)
        message(STATUS "Looking for scotch -- scotch.h not found")
    endif()
endif()


# Looking for lib
# ---------------

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

# Try to find the scotch lib in the given paths
# ----------------------------------------------

set(SCOTCH_libs_to_find "scotch;scotcherrexit")

# call cmake macro to find the lib path
if(SCOTCH_LIBDIR)
    foreach(scotch_lib ${SCOTCH_libs_to_find})
        set(SCOTCH_${scotch_lib}_LIBRARY "SCOTCH_${scotch_lib}_LIBRARY-NOTFOUND")
        find_library(SCOTCH_${scotch_lib}_LIBRARY
            NAMES ${scotch_lib}
            HINTS ${SCOTCH_LIBDIR}) 
    endforeach()
else()
    if(SCOTCH_DIR)
        foreach(scotch_lib ${SCOTCH_libs_to_find})
            set(SCOTCH_${scotch_lib}_LIBRARY "SCOTCH_${scotch_lib}_LIBRARY-NOTFOUND")
            find_library(SCOTCH_${scotch_lib}_LIBRARY
                NAMES ${scotch_lib}
                HINTS ${SCOTCH_DIR}
                PATH_SUFFIXES lib lib32 lib64)
        endforeach()
    else()
        foreach(scotch_lib ${SCOTCH_libs_to_find})
            set(SCOTCH_${scotch_lib}_LIBRARY "SCOTCH_${scotch_lib}_LIBRARY-NOTFOUND")
            find_library(SCOTCH_${scotch_lib}_LIBRARY
                NAMES ${scotch_lib}
                HINTS ${_lib_env})  
        endforeach()
    endif()
endif()

# Print status if not found
# -------------------------
foreach(scotch_lib ${SCOTCH_libs_to_find})
    if (NOT SCOTCH_${scotch_lib}_LIBRARY AND NOT SCOTCH_FIND_QUIETLY)
        Print_Find_Library_Status(scotch ${scotch_lib})
    endif ()
endforeach()

set(SCOTCH_LIBRARIES "")
set(SCOTCH_LIBRARY_DIRS "")
# If found, add path to cmake variable
# ------------------------------------
foreach(scotch_lib ${SCOTCH_libs_to_find})

    if (SCOTCH_${scotch_lib}_LIBRARY)
        get_filename_component(${scotch_lib}_lib_path "${SCOTCH_${scotch_lib}_LIBRARY}" PATH)
        # set cmake variables
        list(APPEND SCOTCH_LIBRARIES "${SCOTCH_${scotch_lib}_LIBRARY}")
        list(APPEND SCOTCH_LIBRARY_DIRS "${${scotch_lib}_lib_path}")
    else ()
        list(APPEND SCOTCH_LIBRARIES "${SCOTCH_${scotch_lib}_LIBRARY}")
        if (NOT SCOTCH_FIND_QUIETLY)
            message(STATUS "Looking for scotch -- lib ${scotch_lib} not found")
        endif()
    endif ()
    
    mark_as_advanced(SCOTCH_${scotch_lib}_LIBRARY)

endforeach()

# check that SCOTCH has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCOTCH DEFAULT_MSG
                                  SCOTCH_LIBRARIES
                                  SCOTCH_INCLUDE_DIRS
                                  SCOTCH_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
