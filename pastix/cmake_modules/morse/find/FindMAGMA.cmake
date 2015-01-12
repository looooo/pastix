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
# The user can give specific paths where to find the libraries adding cmake 
# options at configure (ex: cmake path/to/project -DMAGMA_DIR=path/to/magma):
#  MAGMA_DIR             - Where to find the base directory of magma
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


# Some macros to print status when search for headers and libs
# PrintFindStatus.cmake is in cmake_modules/morse/find directory
include(PrintFindStatus)

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

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE)

    pkg_search_module(MAGMA magma)
    if (NOT MAGMA_FIND_QUIETLY)
        if (MAGMA_FOUND AND MAGMA_LIBRARIES)
            message(STATUS "Looking for MAGMA - found using PkgConfig")
            if(NOT MAGMA_INCLUDE_DIRS)
                message("${Magenta}MAGMA_INCLUDE_DIRS is empty using PkgConfig."
                    "Perhaps the path to magma headers is already present in your"
                    "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
            endif()
        else()
            message("${Magenta}Looking for MAGMA - not found using PkgConfig."
                "Perhaps you should add the directory containing magma.pc"
                "to the PKG_CONFIG_PATH environment variable.${ColourReset}")
        endif()
    endif()

    if (MAGMA_FIND_VERSION_EXACT STREQUAL 1)
        if( NOT (MAGMA_FIND_VERSION_MAJOR STREQUAL MAGMA_VERSION_MAJOR) OR
            NOT (MAGMA_FIND_VERSION_MINOR STREQUAL MAGMA_VERSION_MINOR) )
            if(NOT MAGMA_FIND_QUIETLY)
                message(FATAL_ERROR
                        "MAGMA version found is ${MAGMA_VERSION_STRING}"
                        "when required is ${MAGMA_FIND_VERSION}")
            endif()
        endif()
    else()
        # if the version found is older than the required then error
        if( (MAGMA_FIND_VERSION_MAJOR STRGREATER MAGMA_VERSION_MAJOR) OR
            (MAGMA_FIND_VERSION_MINOR STRGREATER MAGMA_VERSION_MINOR) )
            if(NOT MAGMA_FIND_QUIETLY)
                message(FATAL_ERROR
                        "MAGMA version found is ${MAGMA_VERSION_STRING}"
                        "when required is ${MAGMA_FIND_VERSION} or newer")
            endif()
        endif()
    endif()

endif(PKG_CONFIG_EXECUTABLE)

if(NOT MAGMA_FOUND OR NOT MAGMA_LIBRARIES)

    if (NOT MAGMA_FIND_QUIETLY)
        message(STATUS "Looking for MAGMA - PkgConfig not used")
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
    
    
    # Try to find the magma header in the given paths
    # -------------------------------------------------
    # call cmake macro to find the header path
    if(MAGMA_INCDIR)
        set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
        find_path(MAGMA_magma.h_DIRS
          NAMES magma.h
          HINTS ${MAGMA_INCDIR})
    else()
        if(MAGMA_DIR)
            set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
            find_path(MAGMA_magma.h_DIRS
              NAMES magma.h
              HINTS ${MAGMA_DIR}
              PATH_SUFFIXES include)
        else()
            set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
            find_path(MAGMA_magma.h_DIRS
              NAMES magma.h
              HINTS ${_inc_env})
        endif()
    endif()
    mark_as_advanced(MAGMA_magma.h_DIRS)
    
    # Print status if not found
    # -------------------------
    if (NOT MAGMA_magma.h_DIRS AND NOT MAGMA_FIND_QUIETLY)
        Print_Find_Header_Status(magma magma.h)
    endif ()
    
    # If found, add path to cmake variable
    # ------------------------------------
    if (MAGMA_magma.h_DIRS)
        set(MAGMA_INCLUDE_DIRS "${MAGMA_magma.h_DIRS}")
    else ()
        set(MAGMA_INCLUDE_DIRS "MAGMA_INCLUDE_DIRS-NOTFOUND")
        if(NOT MAGMA_FIND_QUIETLY)
            message(STATUS "Looking for magma -- magma.h not found")
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
    
    # Try to find the magma lib in the given paths
    # ----------------------------------------------
    
    # call cmake macro to find the lib path
    if(MAGMA_LIBDIR)
        set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
        find_library(MAGMA_magma_LIBRARY
            NAMES magma
            HINTS ${MAGMA_LIBDIR})
    else()
        if(MAGMA_DIR)
            set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
            find_library(MAGMA_magma_LIBRARY
                NAMES magma
                HINTS ${MAGMA_DIR}
                PATH_SUFFIXES lib lib32 lib64)
        else()
            set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
            find_library(MAGMA_magma_LIBRARY
                NAMES magma
                HINTS ${_lib_env})
        endif()
    endif()
    mark_as_advanced(MAGMA_magma_LIBRARY)
    
    # Print status if not found
    # -------------------------
    if (NOT MAGMA_magma_LIBRARY AND NOT MAGMA_FIND_QUIETLY)
        Print_Find_Library_Status(magma libmagma)
    endif ()
    
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
        if(NOT MAGMA_FIND_QUIETLY)
            message(STATUS "Looking for magma -- lib magma not found")
        endif()
    endif ()

endif(NOT MAGMA_FOUND OR NOT MAGMA_LIBRARIES)


# check that MAGMA has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MAGMA DEFAULT_MSG
                                  MAGMA_LIBRARIES)
#
# TODO: Add possibility to check for specific functions in the library
#
