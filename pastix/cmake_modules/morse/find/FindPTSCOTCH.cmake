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
# - Find PTSCOTCH include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PTSCOTCH
#               [REQUIRED]             # Fail with error if ptscotch is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )  
# This module finds headers and ptscotch library. 
# Results are reported in variables:
#  PTSCOTCH_FOUND           - True if headers and requested libraries were found
#  PTSCOTCH_INCLUDE_DIRS    - ptscotch include directories
#  PTSCOTCH_LIBRARY_DIRS    - Link directories for ptscotch libraries
#  PTSCOTCH_LIBRARIES       - ptscotch component libraries to be linked
# The user can give specific paths where to find the libraries adding cmake 
# options at configure (ex: cmake path/to/project -DPTSCOTCH=path/to/ptscotch):
#  PTSCOTCH_DIR             - Where to find the base directory of ptscotch
#  PTSCOTCH_INCDIR          - Where to find the header files
#  PTSCOTCH_LIBDIR          - Where to find the library files

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

# PTSCOTCH may depend on MPI and Threads
# try to find it specified as COMPONENTS during the call
if( PTSCOTCH_FIND_COMPONENTS )
    foreach( component ${PTSCOTCH_FIND_COMPONENTS} )
        if(${PTSCOTCH_FIND_REQUIRED_${component}} STREQUAL 1)
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(PTSCOTCH_${component}_FOUND TRUE)
            if( THREADS_FOUND )
                list(APPEND EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
            endif ()
            if (MPI_FOUND)
                # should we have these variables available in gui modes?
                mark_as_advanced(MPI_LIBRARY)
                mark_as_advanced(MPI_EXTRA_LIBRARY)
                list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
                include_directories( ${MPI_C_INCLUDE_PATH} )                
            endif()
        else()
            set(PTSCOTCH_${component}_FOUND FALSE)
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


# Try to find the ptscotch header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(PTSCOTCH_INCDIR)
    set(PTSCOTCH_ptscotch.h_DIRS "PTSCOTCH_ptscotch.h_DIRS-NOTFOUND")
    find_path(PTSCOTCH_ptscotch.h_DIRS
      NAMES ptscotch.h
      HINTS ${PTSCOTCH_INCDIR})
else()
    if(PTSCOTCH_DIR)
        set(PTSCOTCH_ptscotch.h_DIRS "PTSCOTCH_ptscotch.h_DIRS-NOTFOUND")
        find_path(PTSCOTCH_ptscotch.h_DIRS
          NAMES ptscotch.h
          HINTS ${PTSCOTCH_DIR}
          PATH_SUFFIXES include)        
    else()
        set(PTSCOTCH_ptscotch.h_DIRS "PTSCOTCH_ptscotch.h_DIRS-NOTFOUND")
        find_path(PTSCOTCH_ptscotch.h_DIRS
          NAMES ptscotch.h
          HINTS ${_inc_env})
    endif()
endif()
mark_as_advanced(PTSCOTCH_ptscotch.h_DIRS)

# Print status if not found
# -------------------------
if (NOT PTSCOTCH_ptscotch.h_DIRS AND NOT PTSCOTCH_FIND_QUIETLY)
    Print_Find_Header_Status(ptscotch ptscotch.h)
endif ()

# If found, add path to cmake variable
# ------------------------------------
if (PTSCOTCH_ptscotch.h_DIRS)
    set(PTSCOTCH_INCLUDE_DIRS "${PTSCOTCH_ptscotch.h_DIRS}")
else ()
    set(PTSCOTCH_INCLUDE_DIRS "PTSCOTCH_INCLUDE_DIRS-NOTFOUND")
    if (NOT PTSCOTCH_FIND_QUIETLY)
        message(STATUS "Looking for ptscotch -- ptscotch.h not found")
    endif
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

# Try to find the ptscotch lib in the given paths
# ----------------------------------------------

set(PTSCOTCH_libs_to_find "ptscotch;scotch;scotcherrexit")

# call cmake macro to find the lib path
if(PTSCOTCH_LIBDIR)
    foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
        set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
        find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
            NAMES ${ptscotch_lib}
            HINTS ${PTSCOTCH_LIBDIR}) 
    endforeach()
else()
    if(PTSCOTCH_DIR)  
        foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
            set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
            find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
                NAMES ${ptscotch_lib}
                HINTS ${PTSCOTCH_DIR}
                PATH_SUFFIXES lib lib32 lib64)
        endforeach()
    else()
        foreach(scotch_lib ${PTSCOTCH_libs_to_find})
            set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
            find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
                NAMES ${ptscotch_lib}
                HINTS ${_lib_env})  
        endforeach()
    endif()
endif()

# Print status if not found
# -------------------------
foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
    if (NOT PTSCOTCH_${ptscotch_lib}_LIBRARY AND NOT PTSCOTCH_FIND_QUIETLY)
        Print_Find_Library_Status(ptscotch ${ptscotch_lib})
    endif ()
endforeach()

set(PTSCOTCH_LIBRARIES "")
set(PTSCOTCH_LIBRARY_DIRS "")
# If found, add path to cmake variable
# ------------------------------------
foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})

    if (PTSCOTCH_${ptscotch_lib}_LIBRARY)
        get_filename_component(${ptscotch_lib}_lib_path "${PTSCOTCH_${ptscotch_lib}_LIBRARY}" PATH)
        # set cmake variables
        list(APPEND PTSCOTCH_LIBRARIES "${PTSCOTCH_${ptscotch_lib}_LIBRARY}")
        list(APPEND PTSCOTCH_LIBRARY_DIRS "${${ptscotch_lib}_lib_path}")
    else ()
        list(APPEND PTSCOTCH_LIBRARIES "${PTSCOTCH_${ptscotch_lib}_LIBRARY}")
        if (NOT PTSCOTCH_FIND_QUIETLY)
            message(STATUS "Looking for ptscotch -- lib ${ptscotch_lib} not found")
        endif()
    endif ()
    
    mark_as_advanced(PTSCOTCH_${ptscotch_lib}_LIBRARY)
    
endforeach()


# check that PTSCOTCH has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH DEFAULT_MSG
                                  PTSCOTCH_LIBRARIES
                                  PTSCOTCH_INCLUDE_DIRS
                                  PTSCOTCH_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
