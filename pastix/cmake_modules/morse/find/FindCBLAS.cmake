# - Find CBLAS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(CBLAS
#               [REQUIRED]             # Fail with error if cblas is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )
# This module finds headers and cblas library. 
# Results are reported in variables:
#  CBLAS_FOUND           - True if headers and requested libraries were found
#  CBLAS_INCLUDE_DIRS    - cblas include directories
#  CBLAS_LIBRARY_DIRS    - Link directories for cblas libraries
#  CBLAS_LIBRARIES       - cblas component libraries to be linked
# The user can give specific paths where to find the libraries:
#  CBLAS_DIR             - Where to find the base directory of CBLAS
#  CBLAS_INCDIR          - Where to find the header files
#  CBLAS_LIBDIR          - Where to find the library files

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
# PrintFindStatus.cmake is in cmake_modules/morse/find directory of magmamorse
include(PrintFindStatus)

# CBLAS depends on BLAS
# try to find it specified as COMPONENTS during the call
if (CBLAS_FIND_COMPONENTS)
    foreach( component ${CBLAS_FIND_COMPONENTS} )
        if(${CBLAS_FIND_REQUIRED_${component}} STREQUAL 1)
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(CBLAS_${component}_FOUND TRUE)
        else()
            set(CBLAS_${component}_FOUND FALSE)
        endif()
    endforeach()
endif ()

# CBLAS depends on BLAS
if (BLAS_FOUND)

    # check if a cblas function exists in the BLAS lib
    include(CheckFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES "${BLAS_LIBRARIES};${CMAKE_THREAD_LIBS_INIT};${LM}")
    unset(CBLAS_WORKS CACHE)
    check_function_exists(cblas_dscal CBLAS_WORKS)
    mark_as_advanced(CBLAS_WORKS)
    set(CMAKE_REQUIRED_LIBRARIES)
    
    if(CBLAS_WORKS)
        message(STATUS "Looking for cblas: test with blas succeeds")
        # test succeeds: CBLAS is in BLAS
        set(CBLAS_LIBRARIES "${BLAS_LIBRARIES}")
        set(CBLAS_LIBRARY_DIRS "${BLAS_LIBRARY_DIRS}")
        if(BLAS_INCLUDE_DIRS)
            set(CBLAS_INCLUDE_DIRS "${BLAS_INCLUDE_DIRS}")
        endif()        
    else()
        message(STATUS "Looking for cblas : test with blas fails")
        # test fails: try to find CBLAS lib exterior to BLAS
        
        # Try to find CBLAS lib
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
        list(REMOVE_DUPLICATES _inc_env)


        # Try to find the cblas header in the given paths
        # -------------------------------------------------
        # call cmake macro to find the header path
        if(DEFINED CBLAS_INCDIR)
            set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
            find_path(CBLAS_cblas.h_DIRS
            NAMES cblas.h
            HINTS ${CBLAS_INCDIR})
        else()
            if(DEFINED CBLAS_DIR)
                set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
                find_path(CBLAS_cblas.h_DIRS
                NAMES cblas.h
                HINTS ${CBLAS_DIR}
                PATH_SUFFIXES include)
            else()
                set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
                find_path(CBLAS_cblas.h_DIRS
                NAMES cblas.h
                PATHS ${_inc_env})
            endif()
        endif()
        mark_as_advanced(CBLAS_cblas.h_DIRS)
        
        # Print status if not found
        # -------------------------
        if (NOT CBLAS_cblas.h_DIRS)
            Print_Find_Header_Status(cblas cblas.h)
        endif ()        

        # If found, add path to cmake variable
        # ------------------------------------
        if (CBLAS_cblas.h_DIRS)
            set(CBLAS_INCLUDE_DIRS "${CBLAS_cblas.h_DIRS}")
        else ()
            set(CBLAS_INCLUDE_DIRS "CBLAS_INCLUDE_DIRS-NOTFOUND")
            message(STATUS "Looking for cblas -- cblas.h not found")
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

        # Try to find the cblas lib in the given paths
        # ----------------------------------------------

        # call cmake macro to find the lib path
        if(DEFINED CBLAS_LIBDIR)
            set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
            find_library(CBLAS_cblas_LIBRARY
                NAMES cblas
                HINTS ${CBLAS_LIBDIR})
        else()
            if(DEFINED CBLAS_DIR)
                set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
                find_library(CBLAS_cblas_LIBRARY
                    NAMES cblas
                    HINTS ${CBLAS_DIR}
                    PATH_SUFFIXES lib lib32 lib64)
            else()
                set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
                find_library(CBLAS_cblas_LIBRARY
                    NAMES cblas
                    PATHS ${_lib_env})
            endif()
        endif()
        mark_as_advanced(CBLAS_cblas_LIBRARY)
        
        # Print status if not found
        # -------------------------
        if (NOT CBLAS_cblas_LIBRARY)
            Print_Find_Library_Status(cblas libcblas)
        endif ()

        # If found, add path to cmake variable
        # ------------------------------------
        if (CBLAS_cblas_LIBRARY)
            get_filename_component(cblas_lib_path "${CBLAS_cblas_LIBRARY}" PATH)
            # set cmake variables
            set(CBLAS_LIBRARIES    "${CBLAS_cblas_LIBRARY}")
            set(CBLAS_LIBRARY_DIRS "${cblas_lib_path}")
        else ()
            set(CBLAS_LIBRARIES    "CBLAS_LIBRARIES-NOTFOUND")
            set(CBLAS_LIBRARY_DIRS "CBLAS_LIBRARY_DIRS-NOTFOUND")
            message(STATUS "Looking for cblas -- lib cblas not found")
        endif ()        
    endif()
    
else()

    message(STATUS "CBLAS requires BLAS")
    
endif()


# check that CBLAS has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS DEFAULT_MSG
#                                  CBLAS_FOUND
                                  CBLAS_LIBRARIES
                                  CBLAS_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
