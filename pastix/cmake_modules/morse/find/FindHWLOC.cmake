# - Find HWLOC include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(HWLOC
#               [REQUIRED]) # Fail with error if hwloc is not found
# This module finds headers and hwloc library. 
# Results are reported in variables:
#  HWLOC_FOUND           - True if headers and requested libraries were found
#  HWLOC_INCLUDE_DIRS    - hwloc include directories
#  HWLOC_LIBRARY_DIRS    - Link directories for hwloc libraries
#  HWLOC_LIBRARIES       - hwloc component libraries to be linked
# The user can give specific paths where to find the libraries:
#  HWLOC_DIR             - Where to find the base directory of hwloc
#  HWLOC_INCDIR          - Where to find the header files
#  HWLOC_LIBDIR          - Where to find the library files

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

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE)
    string(REPLACE ":" ";" PATH_PKGCONFIGPATH "$ENV{PKG_CONFIG_PATH}")
    set(HWLOC_PKG_FILE_FOUND "HWLOC_PKG_FILE_FOUND-NOTFOUND")
    find_file(HWLOC_PKG_FILE_FOUND
              NAME  hwloc.pc
              PATHS ${PATH_PKGCONFIGPATH})
    mark_as_advanced(HWLOC_PKG_FILE_FOUND)
    if(HWLOC_PKG_FILE_FOUND)
        pkg_search_module(PC_HWLOC hwloc)
    else()
        Print_Find_Pkgconfig_Status(hwloc hwloc.pc ${PATH_PKGCONFIGPATH})
        message(STATUS "Looking for HWLOC - pkgconfig not used")
    endif()

else()
    message(STATUS "Looking for HWLOC - pkgconfig not used")
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

# set paths where to look for
set(PATH_TO_LOOK_FOR "${_inc_env}")

# Try to find the hwloc header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(DEFINED HWLOC_INCDIR)
    set(HWLOC_hwloc.h_DIRS "HWLOC_hwloc.h_DIRS-NOTFOUND")
    find_path(HWLOC_hwloc.h_DIRS
      NAMES hwloc.h
      HINTS ${HWLOC_INCDIR})
else()
    if(DEFINED HWLOC_DIR)
        set(HWLOC_hwloc.h_DIRS "HWLOC_hwloc.h_DIRS-NOTFOUND")
        find_path(HWLOC_hwloc.h_DIRS
          NAMES hwloc.h
          HINTS ${HWLOC_DIR}
          PATH_SUFFIXES include)
    else()
        set(HWLOC_hwloc.h_DIRS "HWLOC_hwloc.h_DIRS-NOTFOUND")
        if(HWLOC_PKG_FILE_FOUND AND PC_HWLOC_INCLUDE_DIRS)
            find_path(HWLOC_hwloc.h_DIRS
            NAMES hwloc.h
            HINTS ${PC_HWLOC_INCLUDE_DIRS})        
        else()
            find_path(HWLOC_hwloc.h_DIRS
            NAMES hwloc.h
            PATHS ${PATH_TO_LOOK_FOR})
        endif()
    endif()
endif()
mark_as_advanced(HWLOC_hwloc.h_DIRS)

# Print status if not found
# -------------------------
if (NOT HWLOC_hwloc.h_DIRS)
    Print_Find_Header_Status(hwloc hwloc.h)
endif ()

# Add path to cmake variable
# ------------------------------------
if (HWLOC_hwloc.h_DIRS)
    set(HWLOC_INCLUDE_DIRS "${HWLOC_hwloc.h_DIRS}")
else ()
    set(HWLOC_INCLUDE_DIRS "HWLOC_INCLUDE_DIRS-NOTFOUND")
    message(STATUS "Looking for hwloc -- hwloc.h not found")
endif ()

if (HWLOC_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES HWLOC_INCLUDE_DIRS)
endif ()


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

# set paths where to look for
set(PATH_TO_LOOK_FOR "${_lib_env}")

# Try to find the hwloc lib in the given paths
# ----------------------------------------------

# call cmake macro to find the lib path
if(DEFINED HWLOC_LIBDIR)
    set(HWLOC_hwloc_LIBRARY "HWLOC_hwloc_LIBRARY-NOTFOUND")
    find_library(HWLOC_hwloc_LIBRARY
        NAMES hwloc
        HINTS ${HWLOC_LIBDIR})
else()
    if(DEFINED HWLOC_DIR)
        set(HWLOC_hwloc_LIBRARY "HWLOC_hwloc_LIBRARY-NOTFOUND")
        find_library(HWLOC_hwloc_LIBRARY
            NAMES hwloc
            HINTS ${HWLOC_DIR}
            PATH_SUFFIXES lib lib32 lib64)
    else()
        set(HWLOC_hwloc_LIBRARY "HWLOC_hwloc_LIBRARY-NOTFOUND")    
        if(HWLOC_PKG_FILE_FOUND AND PC_HWLOC_LIBRARY_DIRS)
            find_library(HWLOC_hwloc_LIBRARY
                NAMES hwloc
                HINTS ${PC_HWLOC_LIBRARY_DIRS})
        else()
            set(HWLOC_hwloc_LIBRARY "HWLOC_hwloc_LIBRARY-NOTFOUND")
            find_library(HWLOC_hwloc_LIBRARY
                NAMES hwloc
                PATHS ${PATH_TO_LOOK_FOR})        
        endif()
    endif()
endif()
mark_as_advanced(HWLOC_hwloc_LIBRARY)

# Print status if not found
# -------------------------
if (NOT HWLOC_hwloc_LIBRARY)
    Print_Find_Library_Status(hwloc libhwloc)
endif ()

# If found, add path to cmake variable
# ------------------------------------
if (HWLOC_hwloc_LIBRARY)
    get_filename_component(hwloc_lib_path ${HWLOC_hwloc_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(HWLOC_LIBRARIES    "${HWLOC_hwloc_LIBRARY}")
    set(HWLOC_LIBRARY_DIRS "${hwloc_lib_path}")
else ()
    set(HWLOC_LIBRARIES    "HWLOC_LIBRARIES-NOTFOUND")
    set(HWLOC_LIBRARY_DIRS "HWLOC_LIBRARY_DIRS-NOTFOUND")
    message(STATUS "Looking for hwloc -- lib hwloc not found")
endif ()

if (HWLOC_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES HWLOC_LIBRARY_DIRS)
endif ()


# check that HWLOC has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC DEFAULT_MSG
#                                  HWLOC_FOUND
                                  HWLOC_LIBRARIES
                                  HWLOC_INCLUDE_DIRS
                                  HWLOC_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
