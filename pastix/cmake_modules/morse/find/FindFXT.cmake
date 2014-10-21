# - Find FXT include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(FXT
#               [REQUIRED]) # Fail with error if fxt is not found
# This module finds headers and fxt library. 
# Results are reported in variables:
#  FXT_FOUND           - True if headers and requested libraries were found
#  FXT_INCLUDE_DIRS    - fxt include directories
#  FXT_LIBRARY_DIRS    - Link directories for fxt libraries
#  FXT_LIBRARIES       - fxt component libraries to be linked
# The user can give specific paths where to find the libraries:
#  FXT_DIR             - Where to find the base directory of fxt
#  FXT_INCDIR          - Where to find the header files
#  FXT_LIBDIR          - Where to find the library files

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
    set(FXT_PKG_FILE_FOUND "FXT_PKG_FILE_FOUND-NOTFOUND")
    find_file(FXT_PKG_FILE_FOUND
              NAME  fxt.pc
              PATHS ${PATH_PKGCONFIGPATH})
    mark_as_advanced(FXT_PKG_FILE_FOUND)
    if(FXT_PKG_FILE_FOUND)
        pkg_search_module(PC_FXT fxt)
    else()
        Print_Find_Pkgconfig_Status(fxt fxt.pc ${PATH_PKGCONFIGPATH})
        message(STATUS "Looking for FXT - pkgconfig not used")
    endif()

else()
    message(STATUS "Looking for FXT - pkgconfig not used")
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

# Try to find the fxt header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(DEFINED FXT_INCDIR)
    set(FXT_fxt.h_DIRS "FXT_fxt.h_DIRS-NOTFOUND")
    find_path(FXT_fxt.h_DIRS
      NAMES fxt.h
      HINTS ${FXT_INCDIR})
else()
    if(DEFINED FXT_DIR)
        set(FXT_fxt.h_DIRS "FXT_fxt.h_DIRS-NOTFOUND")
        find_path(FXT_fxt.h_DIRS
          NAMES fxt.h
          HINTS ${FXT_DIR}
          PATH_SUFFIXES include)
    else()
        set(FXT_fxt.h_DIRS "FXT_fxt.h_DIRS-NOTFOUND")
        if(FXT_PKG_FILE_FOUND AND PC_FXT_INCLUDE_DIRS)
            find_path(FXT_fxt.h_DIRS
            NAMES fxt.h
            HINTS ${PC_FXT_INCLUDE_DIRS})        
        else()
            find_path(FXT_fxt.h_DIRS
            NAMES fxt.h
            PATHS ${PATH_TO_LOOK_FOR})
        endif()
    endif()
endif()
mark_as_advanced(FXT_fxt.h_DIRS)

# Print status if not found
# -------------------------
if (NOT FXT_fxt.h_DIRS)
    Print_Find_Header_Status(fxt fxt.h)
endif ()

# Add path to cmake variable
# ------------------------------------
if (FXT_fxt.h_DIRS)
    set(FXT_INCLUDE_DIRS "${FXT_fxt.h_DIRS}")
else ()
    set(FXT_INCLUDE_DIRS "FXT_INCLUDE_DIRS-NOTFOUND")
    message(STATUS "Looking for fxt -- fxt.h not found")
endif ()

if (FXT_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES FXT_INCLUDE_DIRS)
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

# Try to find the fxt lib in the given paths
# ----------------------------------------------

# call cmake macro to find the lib path
if(DEFINED FXT_LIBDIR)
    set(FXT_fxt_LIBRARY "FXT_fxt_LIBRARY-NOTFOUND")
    find_library(FXT_fxt_LIBRARY
        NAMES fxt
        HINTS ${FXT_LIBDIR})
else()
    if(DEFINED FXT_DIR)
        set(FXT_fxt_LIBRARY "FXT_fxt_LIBRARY-NOTFOUND")
        find_library(FXT_fxt_LIBRARY
            NAMES fxt
            HINTS ${FXT_DIR}
            PATH_SUFFIXES lib lib32 lib64)
    else()
        set(FXT_fxt_LIBRARY "FXT_fxt_LIBRARY-NOTFOUND")
        if(FXT_PKG_FILE_FOUND AND PC_FXT_LIBRARY_DIRS)
            find_library(FXT_fxt_LIBRARY
                NAMES fxt
                HINTS ${PC_FXT_LIBRARY_DIRS})        
        else()        
            find_library(FXT_fxt_LIBRARY
                NAMES fxt
                PATHS ${PATH_TO_LOOK_FOR})
        endif()
    endif()
endif()
mark_as_advanced(FXT_fxt_LIBRARY)

# Print status if not found
# -------------------------
if (NOT FXT_fxt_LIBRARY)
    Print_Find_Library_Status(fxt libfxt)
endif ()

# If found, add path to cmake variable
# ------------------------------------
if (FXT_fxt_LIBRARY)
    get_filename_component(fxt_lib_path ${FXT_fxt_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(FXT_LIBRARIES    "${FXT_fxt_LIBRARY}")
    set(FXT_LIBRARY_DIRS "${fxt_lib_path}")
else ()
    set(FXT_LIBRARIES    "FXT_LIBRARIES-NOTFOUND")
    set(FXT_LIBRARY_DIRS "FXT_LIBRARY_DIRS-NOTFOUND")
    message(STATUS "Looking for fxt -- lib fxt not found")
endif ()

if (FXT_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES FXT_LIBRARY_DIRS)
endif ()


# check that FXT has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FXT DEFAULT_MSG
#                                  FXT_FOUND
                                  FXT_LIBRARIES
                                  FXT_INCLUDE_DIRS
                                  FXT_LIBRARY_DIRS)
#
# TODO: Add possibility to check for specific functions in the library
#
