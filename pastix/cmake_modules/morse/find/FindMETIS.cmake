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
# - Find METIS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(METIS
#               [REQUIRED]             # Fail with error if metis is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )  
# This module finds headers and metis library. 
# Results are reported in variables:
#  METIS_FOUND           - True if headers and requested libraries were found
#  METIS_INCLUDE_DIRS    - metis include directories
#  METIS_LIBRARY_DIRS    - Link directories for metis libraries
#  METIS_LIBRARIES       - metis component libraries to be linked
# The user can give specific paths where to find the libraries adding cmake 
# options at configure (ex: cmake path/to/project -DMETIS_DIR=path/to/metis):
#  METIS_DIR             - Where to find the base directory of metis
#  METIS_INCDIR          - Where to find the header files
#  METIS_LIBDIR          - Where to find the library files

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


# Try to find the metis header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(METIS_INCDIR)
    set(METIS_metis.h_DIRS "METIS_metis.h_DIRS-NOTFOUND")
    find_path(METIS_metis.h_DIRS
      NAMES metis.h
      HINTS ${METIS_INCDIR})
else()
    if(METIS_DIR)
        set(METIS_metis.h_DIRS "METIS_metis.h_DIRS-NOTFOUND")
        find_path(METIS_metis.h_DIRS
          NAMES metis.h
          HINTS ${METIS_DIR}
          PATH_SUFFIXES include)        
    else()
        set(METIS_metis.h_DIRS "METIS_metis.h_DIRS-NOTFOUND")
        find_path(METIS_metis.h_DIRS
          NAMES metis.h
          HINTS ${_inc_env})
    endif()
endif()
mark_as_advanced(METIS_metis.h_DIRS)

# Print status if not found
# -------------------------
if (NOT METIS_metis.h_DIRS AND NOT METIS_FIND_QUIETLY)
    Print_Find_Header_Status(metis metis.h)
endif ()

# If found, add path to cmake variable
# ------------------------------------
if (METIS_metis.h_DIRS)
    set(METIS_INCLUDE_DIRS "${METIS_metis.h_DIRS}")
else ()
    set(METIS_INCLUDE_DIRS "METIS_INCLUDE_DIRS-NOTFOUND")
    if(NOT METIS_FIND_QUIETLY)
        message(STATUS "Looking for metis -- metis.h not found")
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

# Try to find the metis lib in the given paths
# ----------------------------------------------
# call cmake macro to find the lib path
if(METIS_LIBDIR)
    set(METIS_metis_LIBRARY "METIS_metis_LIBRARY-NOTFOUND")
    find_library(METIS_metis_LIBRARY
        NAMES metis
        HINTS ${METIS_LIBDIR})      
else()
    if(METIS_DIR)   
        set(METIS_metis_LIBRARY "METIS_metis_LIBRARY-NOTFOUND")
        find_library(METIS_metis_LIBRARY
            NAMES metis
            HINTS ${METIS_DIR}
            PATH_SUFFIXES lib lib32 lib64)
    else()
        set(METIS_metis_LIBRARY "METIS_metis_LIBRARY-NOTFOUND")
        find_library(METIS_metis_LIBRARY
            NAMES metis
            HINTS ${_lib_env})        
    endif()
endif()
mark_as_advanced(METIS_metis_LIBRARY)

# Print status if not found
# -------------------------
if (NOT METIS_metis_LIBRARY AND NOT METIS_FIND_QUIETLY)
    Print_Find_Library_Status(metis libmetis)
endif ()

# If found, add path to cmake variable
# ------------------------------------
if (METIS_metis_LIBRARY)
    get_filename_component(metis_lib_path "${METIS_metis_LIBRARY}" PATH)
    # set cmake variables
    set(METIS_LIBRARIES    "${METIS_metis_LIBRARY}")
    set(METIS_LIBRARY_DIRS "${metis_lib_path}")
else ()
    set(METIS_LIBRARIES    "METIS_LIBRARIES-NOTFOUND")
    set(METIS_LIBRARY_DIRS "METIS_LIBRARY_DIRS-NOTFOUND")
    if(NOT METIS_FIND_QUIETLY)
        message(STATUS "Looking for metis -- lib metis not found")
    endif()
endif ()


# check that METIS has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS DEFAULT_MSG
                                  METIS_LIBRARIES
                                  METIS_INCLUDE_DIRS
                                  METIS_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
