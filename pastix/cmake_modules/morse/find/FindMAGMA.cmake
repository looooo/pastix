# - Find MAGMA include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(MAGMA
#               [REQUIRED]             # Fail with error if magma is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )  
# This module finds headers and magma library. 
# Results are reported in variables:
#  MAGMA_FOUND           - True if headers and requested libraries were found
#  MAGMA_INCLUDE_DIRS    - magma include directories
#  MAGMA_LIBRARY_DIRS    - Link directories for magma libraries
#  MAGMA_LIBRARIES       - magma component libraries to be linked
# The user can give specific paths where to find the libraries:
#  MAGMA_DIR             - Where to find the base directory of MAGMA
#  MAGMA_INCDIR          - Where to find the header files
#  MAGMA_LIBDIR          - Where to find the library files

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


# MAGMA may depend on CUDA
# try to find it specified as COMPONENTS during the call
if( MAGMA_FIND_COMPONENTS )
    foreach( component ${MAGMA_FIND_COMPONENTS} )
        if(${MAGMA_FIND_REQUIRED_${component}} STREQUAL 1)
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(MAGMA_${component}_FOUND TRUE)
            # should we have these variables available in gui modes?
            if (CUDA_FOUND)
                mark_as_advanced(CUDA_BUILD_CUBIN)
                mark_as_advanced(CUDA_BUILD_EMULATION)
                mark_as_advanced(CUDA_SDK_ROOT_DIR)
                mark_as_advanced(CUDA_TOOLKIT_ROOT_DIR)
                mark_as_advanced(CUDA_VERBOSE_BUILD)
            endif()
        else()
            set(MAGMA_${component}_FOUND FALSE)
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
list(REMOVE_DUPLICATES _inc_env)


# Try to find the magma header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(DEFINED MAGMA_INCDIR)
    set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
    find_path(MAGMA_magma.h_DIRS
      NAMES magma.h
      HINTS ${MAGMA_INCDIR})
else()
    if(DEFINED MAGMA_DIR)
        set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
        find_path(MAGMA_magma.h_DIRS
          NAMES magma.h
          HINTS ${MAGMA_DIR}
          PATH_SUFFIXES include)
    else()
        set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
        find_path(MAGMA_magma.h_DIRS
          NAMES magma.h
          PATHS ${_inc_env})
    endif()
endif()
mark_as_advanced(MAGMA_magma.h_DIRS)

# If found, add path to cmake variable
# ------------------------------------
if (MAGMA_magma.h_DIRS)
    set(MAGMA_INCLUDE_DIRS "${MAGMA_magma.h_DIRS}")
else ()
    set(MAGMA_INCLUDE_DIRS "MAGMA_INCLUDE_DIRS-NOTFOUND")
    message(STATUS "Looking for magma -- magma.h not found")
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
    list(APPEND _lib_env "/usr/local/lib64")
    list(APPEND _lib_env "/usr/lib64")
    list(APPEND _lib_env "/usr/local/lib")
    list(APPEND _lib_env "/usr/lib")
endif()

# Try to find the magma lib in the given paths
# ----------------------------------------------

# call cmake macro to find the lib path
if(DEFINED MAGMA_LIBDIR)
    set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
    find_library(MAGMA_magma_LIBRARY
        NAMES magma
        HINTS ${MAGMA_LIBDIR})
else()
    if(DEFINED MAGMA_DIR)
        set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
        find_library(MAGMA_magma_LIBRARY
            NAMES magma
            HINTS ${MAGMA_DIR}
            PATH_SUFFIXES lib lib32 lib64)
    else()
        set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
        find_library(MAGMA_magma_LIBRARY
            NAMES magma
            PATHS ${_lib_env})
    endif()
endif()
mark_as_advanced(MAGMA_magma_LIBRARY)

# If found, add path to cmake variable
# ------------------------------------
if (MAGMA_magma_LIBRARY)
    get_filename_component(magma_lib_path "${MAGMA_magma_LIBRARY}" PATH)
    # set cmake variables
    set(MAGMA_LIBRARIES    "${MAGMA_magma_LIBRARY}")
    set(MAGMA_LIBRARY_DIRS "${magma_lib_path}")
else ()
    set(MAGMA_LIBRARIES    "MAGMA_LIBRARIES-NOTFOUND")
    set(MAGMA_LIBRARY_DIRS "MAGMA_LIBRARY_DIRS-NOTFOUND")
    message(STATUS "Looking for magma -- lib magma not found")
endif ()


# check that MAGMA has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MAGMA DEFAULT_MSG
#                                  MAGMA_FOUND
                                  MAGMA_LIBRARIES
                                  MAGMA_INCLUDE_DIRS
                                  MAGMA_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
