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
# - Find LAPACKE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(LAPACKE
#               [REQUIRED]             # Fail with error if lapacke is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )
# This module finds headers and lapacke library. 
# Results are reported in variables:
#  LAPACKE_FOUND           - True if headers and requested libraries were found
#  LAPACKE_INCLUDE_DIRS    - lapacke include directories
#  LAPACKE_LIBRARY_DIRS    - Link directories for lapacke libraries
#  LAPACKE_LIBRARIES       - lapacke component libraries to be linked
# The user can give specific paths where to find the libraries adding cmake 
# options at configure (ex: cmake path/to/project -DLAPACKE_DIR=path/to/lapacke):
#  LAPACKE_DIR             - Where to find the base directory of lapacke
#  LAPACKE_INCDIR          - Where to find the header files
#  LAPACKE_LIBDIR          - Where to find the library files
# LAPACKE could be directly embedded in LAPACK library (ex: Intel MKL) so that 
# we test a lapacke function with the lapack libraries found and set LAPACKE 
# variables to LAPACK ones if test is successful. To skip this feature and 
# look for a stand alone lapacke, please add the following in your 
# CMakeLists.txt before to call find_package(LAPACKE):
# set(LAPACKE_STANDALONE TRUE)

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

# LAPACKE depends on LAPACK
# try to find it specified as COMPONENTS during the call
if (LAPACKE_FIND_COMPONENTS)
    foreach( component ${LAPACKE_FIND_COMPONENTS} )
        if(${LAPACKE_FIND_REQUIRED_${component}} STREQUAL 1)
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(LAPACKE_${component}_FOUND TRUE)
        else()
            set(LAPACKE_${component}_FOUND FALSE)
        endif()
    endforeach()
endif ()


# LAPACKE depends on LAPACK
if (LAPACK_FOUND)

    if (NOT LAPACKE_STANDALONE)
        # check if a lapacke function exists in the LAPACK lib
        include(CheckFunctionExists)
        set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LIBRARIES};${CMAKE_THREAD_LIBS_INIT};${LM}")
        unset(LAPACKE_WORKS CACHE)
        check_function_exists(LAPACKE_dgeqrf LAPACKE_WORKS)
        mark_as_advanced(LAPACKE_WORKS)
        set(CMAKE_REQUIRED_LIBRARIES)
    
        if(LAPACKE_WORKS)
            if(NOT LAPACKE_FIND_QUIETLY)
                message(STATUS "Looking for lapacke: test with lapack succeeds")
            endif()
            # test succeeds: LAPACKE is in LAPACK
            set(LAPACKE_LIBRARIES "${LAPACK_LIBRARIES}")
            set(LAPACKE_LIBRARY_DIRS "${LAPACK_LIBRARY_DIRS}")
            if(LAPACK_INCLUDE_DIRS)
                set(LAPACKE_INCLUDE_DIRS "${LAPACK_INCLUDE_DIRS}")
            endif()
        endif()
    endif (NOT LAPACKE_STANDALONE)
        
    if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
   
        if(NOT LAPACKE_WORKS AND NOT LAPACKE_FIND_QUIETLY)
            message(STATUS "Looking for lapacke : test with lapack fails")
        endif()
        # test fails: try to find LAPACKE lib exterior to LAPACK
        
        # Try to find LAPACKE lib
        #######################
        
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


        # Try to find the lapacke header in the given paths
        # -------------------------------------------------
        # call cmake macro to find the header path
        if(LAPACKE_INCDIR)
            set(LAPACKE_lapacke.h_DIRS "LAPACKE_lapacke.h_DIRS-NOTFOUND")
            find_path(LAPACKE_lapacke.h_DIRS
            NAMES lapacke.h
            HINTS ${LAPACKE_INCDIR})
        else()
            if(LAPACKE_DIR)
                set(LAPACKE_lapacke.h_DIRS "LAPACKE_lapacke.h_DIRS-NOTFOUND")
                find_path(LAPACKE_lapacke.h_DIRS
                NAMES lapacke.h
                HINTS ${LAPACKE_DIR}
                PATH_SUFFIXES include)
            else()
                set(LAPACKE_lapacke.h_DIRS "LAPACKE_lapacke.h_DIRS-NOTFOUND")
                find_path(LAPACKE_lapacke.h_DIRS
                NAMES lapacke.h
                HINTS ${_inc_env})
            endif()
        endif()
        mark_as_advanced(LAPACKE_lapacke.h_DIRS)
        
        # Print status if not found
        # -------------------------
        if (NOT LAPACKE_lapacke.h_DIRS)
            Print_Find_Header_Status(lapacke lapacke.h)
        endif ()        

        # If found, add path to cmake variable
        # ------------------------------------
        if (LAPACKE_lapacke.h_DIRS)
            set(LAPACKE_INCLUDE_DIRS "${LAPACKE_lapacke.h_DIRS}")
        else ()
            set(LAPACKE_INCLUDE_DIRS "LAPACKE_INCLUDE_DIRS-NOTFOUND")
            if(NOT LAPACKE_FIND_QUIETLY)
                message(STATUS "Looking for lapacke -- lapacke.h not found")
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

        # Try to find the lapacke lib in the given paths
        # ----------------------------------------------

        # call cmake macro to find the lib path
        if(LAPACKE_LIBDIR)
            set(LAPACKE_lapacke_LIBRARY "LAPACKE_lapacke_LIBRARY-NOTFOUND")
            find_library(LAPACKE_lapacke_LIBRARY
                NAMES lapacke
                HINTS ${LAPACKE_LIBDIR})
        else()
            if(LAPACKE_DIR)
                set(LAPACKE_lapacke_LIBRARY "LAPACKE_lapacke_LIBRARY-NOTFOUND")
                find_library(LAPACKE_lapacke_LIBRARY
                    NAMES lapacke
                    HINTS ${LAPACKE_DIR}
                    PATH_SUFFIXES lib lib32 lib64)
            else()
                set(LAPACKE_lapacke_LIBRARY "LAPACKE_lapacke_LIBRARY-NOTFOUND")
                find_library(LAPACKE_lapacke_LIBRARY
                    NAMES lapacke
                    HINTS ${_lib_env})
            endif()
        endif()
        mark_as_advanced(LAPACKE_lapacke_LIBRARY)
        
        # Print status if not found
        # -------------------------
        if (NOT LAPACKE_lapacke_LIBRARY AND NOT LAPACKE_FIND_QUIETLY)
            Print_Find_Library_Status(lapacke liblapacke)
        endif ()

        # If found, add path to cmake variable
        # ------------------------------------
        if (LAPACKE_lapacke_LIBRARY)
            get_filename_component(lapacke_lib_path "${LAPACKE_lapacke_LIBRARY}" PATH)
            # set cmake variables
            set(LAPACKE_LIBRARIES    "${LAPACKE_lapacke_LIBRARY}")
            set(LAPACKE_LIBRARY_DIRS "${lapacke_lib_path}")
        else ()
            set(LAPACKE_LIBRARIES    "LAPACKE_LIBRARIES-NOTFOUND")
            set(LAPACKE_LIBRARY_DIRS "LAPACKE_LIBRARY_DIRS-NOTFOUND")
            if (NOT LAPACKE_FIND_QUIETLY)
                message(STATUS "Looking for lapacke -- lib lapacke not found")
            endif()
        endif ()
        
    endif (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
    
else(LAPACK_FOUND)

    if (NOT LAPACKE_FIND_QUIETLY)
        message(STATUS "LAPACKE requires LAPACK but LAPACK has not been found."
            "Please look for LAPACK first.")
    endif()
        
endif(LAPACK_FOUND)


# check that LAPACKE has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE DEFAULT_MSG
                                  LAPACKE_LIBRARIES
                                  LAPACKE_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
