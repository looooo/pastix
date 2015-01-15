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
# - Find PARMETIS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PARMETIS
#               [REQUIRED]             # Fail with error if parmetis is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )  
# This module finds headers and parmetis library. 
# Results are reported in variables:
#  PARMETIS_FOUND           - True if headers and requested libraries were found
#  PARMETIS_INCLUDE_DIRS    - parmetis include directories
#  PARMETIS_LIBRARY_DIRS    - Link directories for parmetis libraries
#  PARMETIS_LIBRARIES       - parmetis component libraries to be linked
# The user can give specific paths where to find the libraries adding cmake 
# options at configure (ex: cmake path/to/project -DPARMETIS_DIR=path/to/parmetis):
#  PARMETIS_DIR             - Where to find the base directory of parmetis
#  PARMETIS_INCDIR          - Where to find the header files
#  PARMETIS_LIBDIR          - Where to find the library files

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

# PARMETIS depends on MPI
# try to find it specified as COMPONENTS during the call
if( PARMETIS_FIND_COMPONENTS )
    foreach( component ${PARMETIS_FIND_COMPONENTS} )
        if(${PARMETIS_FIND_REQUIRED_${component}} STREQUAL 1)
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(PARMETIS_${component}_FOUND TRUE)
            # should we have these variables available in gui modes?
            if (MPI_FOUND)
                mark_as_advanced(MPI_LIBRARY)
                mark_as_advanced(MPI_EXTRA_LIBRARY)
            endif()
        else()
            set(PARMETIS_${component}_FOUND FALSE)
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


# Try to find the parmetis header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(PARMETIS_INCDIR)
    set(PARMETIS_parmetis.h_DIRS "PARMETIS_parmetis.h_DIRS-NOTFOUND")
    find_path(PARMETIS_parmetis.h_DIRS
      NAMES parmetis.h
      HINTS ${PARMETIS_INCDIR})
else()
    if(PARMETIS_DIR)
        set(PARMETIS_parmetis.h_DIRS "PARMETIS_parmetis.h_DIRS-NOTFOUND")
        find_path(PARMETIS_parmetis.h_DIRS
          NAMES parmetis.h
          HINTS ${PARMETIS_DIR}
          PATH_SUFFIXES include)        
    else()
        set(PARMETIS_parmetis.h_DIRS "PARMETIS_parmetis.h_DIRS-NOTFOUND")
        find_path(PARMETIS_parmetis.h_DIRS
          NAMES parmetis.h
          HINTS ${_inc_env})
    endif()
endif()
mark_as_advanced(PARMETIS_parmetis.h_DIRS)

# Print status if not found
# -------------------------
if (NOT PARMETIS_parmetis.h_DIRS AND NOT PARMETIS_FIND_QUIETLY)
    Print_Find_Header_Status(parmetis parmetis.h)
endif ()

# If found, add path to cmake variable
# ------------------------------------
if (PARMETIS_parmetis.h_DIRS AND NOT PARMETIS_FIND_QUIETLY)
    set(PARMETIS_INCLUDE_DIRS "${PARMETIS_parmetis.h_DIRS}")
else ()
    set(PARMETIS_INCLUDE_DIRS "PARMETIS_INCLUDE_DIRS-NOTFOUND")
    if(NOT PARMETIS_FIND_QUIETLY)
        message(STATUS "Looking for parmetis -- parmetis.h not found")
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

# Try to find the parmetis lib in the given paths
# ----------------------------------------------
# call cmake macro to find the lib path
if(PARMETIS_LIBDIR)
    set(PARMETIS_parmetis_LIBRARY "PARMETIS_parmetis_LIBRARY-NOTFOUND")
    find_library(PARMETIS_parmetis_LIBRARY
        NAMES parmetis
        HINTS ${PARMETIS_LIBDIR})      
else()
    if(PARMETIS_DIR)   
        set(PARMETIS_parmetis_LIBRARY "PARMETIS_parmetis_LIBRARY-NOTFOUND")
        find_library(PARMETIS_parmetis_LIBRARY
            NAMES parmetis
            HINTS ${PARMETIS_DIR}
            PATH_SUFFIXES lib lib32 lib64)
    else()
        set(PARMETIS_parmetis_LIBRARY "PARMETIS_parmetis_LIBRARY-NOTFOUND")
        find_library(PARMETIS_parmetis_LIBRARY
            NAMES parmetis
            HINTS ${_lib_env})        
    endif()
endif()
mark_as_advanced(PARMETIS_parmetis_LIBRARY)

# Print status if not found
# -------------------------
if (NOT PARMETIS_parmetis_LIBRARY AND NOT PARMETIS_FIND_QUIETLY)
    Print_Find_Library_Status(parmetis libparmetis)
endif ()

# If found, add path to cmake variable
# ------------------------------------
if (PARMETIS_parmetis_LIBRARY)
    get_filename_component(parmetis_lib_path "${PARMETIS_parmetis_LIBRARY}" PATH)
    # set cmake variables
    set(PARMETIS_LIBRARIES    "${PARMETIS_parmetis_LIBRARY}")
    set(PARMETIS_LIBRARY_DIRS "${parmetis_lib_path}")
else ()
    set(PARMETIS_LIBRARIES    "PARMETIS_LIBRARIES-NOTFOUND")
    set(PARMETIS_LIBRARY_DIRS "PARMETIS_LIBRARY_DIRS-NOTFOUND")
    if (NOT PARMETIS_FIND_QUIETLY)
        message(STATUS "Looking for parmetis -- lib parmetis not found")
    endif
endif ()


# check that PARMETIS has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARMETIS DEFAULT_MSG
                                  PARMETIS_LIBRARIES
                                  PARMETIS_INCLUDE_DIRS
                                  PARMETIS_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
