# - Find STARPU include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(STARPU
#    [version] [EXACT]      # Minimum or EXACT version e.g. 1.1
#    [REQUIRED]             # Fail with error if starpu is not found
#    [COMPONENTS <libs>...] # required dependencies
#    )
# This module finds headers and starpu libraries.
# Results are reported in variables:
#  STARPU_FOUND                  - True if headers and requested libraries were found
#  STARPU_INCLUDE_DIRS           - starpu include directories
#  STARPU_LIBRARY_DIRS           - Link directories for starpu libraries
#  STARPU_SHM_LIBRARIES          - starpu component libraries to be linked (without mpi)
#  STARPU_MPI_LIBRARIES          - starpu component libraries to be linked (with mpi)
#  STARPU_${component}_FOUND     - TRUE if component has been found
#  STARPU_VERSION_STRING         - A human-readable string containing the version of the package found
#  STARPU_VERSION_MAJOR          - The major version of the package found
#  STARPU_VERSION_MINOR          - The minor version of the package found
# The user can give specific paths where to find the libraries:
#  STARPU_DIR                    - Where to find the base directory of STARPU
#  STARPU_INCDIR                 - Where to find the header files
#  STARPU_LIBDIR                 - Where to find the library files


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

# STARPU may depend on other packages (HWLOC, MPI, CUDA, ...)
# try to find them if specified as COMPONENTS during the call
set(STARPU_LOOK_FOR_MPI FALSE)
set(STARPU_LOOK_FOR_CUDA FALSE)
if( STARPU_FIND_COMPONENTS )
    foreach( component ${STARPU_FIND_COMPONENTS} )
        if(${component} STREQUAL "MPI")
            set(STARPU_LOOK_FOR_MPI TRUE)
        elseif(${component} STREQUAL "CUDA")
            set(STARPU_LOOK_FOR_CUDA TRUE)
        endif()
        if(${STARPU_FIND_REQUIRED_${component}} STREQUAL 1)
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(STARPU_${component}_FOUND TRUE)
            # should we have these variables available in gui modes?
            # lets hide them
            if (MPI_FOUND)
                mark_as_advanced(MPI_LIBRARY)
                mark_as_advanced(MPI_EXTRA_LIBRARY)
            endif()
            if (CUDA_FOUND)
                mark_as_advanced(CUDA_BUILD_CUBIN)
                mark_as_advanced(CUDA_BUILD_EMULATION)
                mark_as_advanced(CUDA_SDK_ROOT_DIR)
                mark_as_advanced(CUDA_TOOLKIT_ROOT_DIR)
                mark_as_advanced(CUDA_VERBOSE_BUILD)
            endif()
        else()
            set(STARPU_${component}_FOUND FALSE)
        endif()
    endforeach()
endif()


# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE)
    string(REPLACE ":" ";" PATH_PKGCONFIGPATH "$ENV{PKG_CONFIG_PATH}")
    set(STARPU_SHM_PKG_FILE_FOUND "STARPU_SHM_PKG_FILE_FOUND-NOTFOUND")
    find_file(STARPU_SHM_PKG_FILE_FOUND
              NAME  libstarpu.pc
              PATHS ${PATH_PKGCONFIGPATH})
    mark_as_advanced(STARPU_SHM_PKG_FILE_FOUND)
    if(STARPU_SHM_PKG_FILE_FOUND)
        pkg_search_module(PC_STARPU_SHM libstarpu)
    else()
        Print_Find_Pkgconfig_Status(starpu libstarpu.pc ${PATH_PKGCONFIGPATH})
        message(STATUS "Looking for STARPU without MPI - pkgconfig not used")
    endif()
    if(STARPU_LOOK_FOR_MPI)
        set(STARPU_MPI_PKG_FILE_FOUND "STARPU_MPI_PKG_FILE_FOUND-NOTFOUND")
        find_file(STARPU_MPI_PKG_FILE_FOUND
                  NAME  libstarpumpi.pc
                  PATHS ${PATH_PKGCONFIGPATH})
        mark_as_advanced(STARPU_MPI_PKG_FILE_FOUND)
        if(STARPU_MPI_PKG_FILE_FOUND)
            pkg_search_module(PC_STARPU_MPI libstarpumpi)
        else()
            Print_Find_Pkgconfig_Status(starpu libstarpumpi.pc 
                                        ${PATH_PKGCONFIGPATH})
            message(STATUS "Looking for STARPU with MPI- pkgconfig not used")
        endif()
    endif()

else()
    message(STATUS "Looking for STARPU - pkgconfig not used")
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
if(STARPU_SHM_PKG_FILE_FOUND)
    set(PATH_TO_LOOK_FOR "${PC_STARPU_SHM_INCLUDE_DIRS}")
endif()

# Try to find the version of StarPU in starpu_config.h file
set(STARPU_hdrs_to_find "starpu_config.h")

# call cmake macro to find the header path
if(DEFINED STARPU_INCDIR)
    foreach(starpu_hdr ${STARPU_hdrs_to_find})
        set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
        find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
                  NAMES ${starpu_hdr}
                  HINTS ${STARPU_INCDIR})
    endforeach()
else()
    if(DEFINED STARPU_DIR)
        foreach(starpu_hdr ${STARPU_hdrs_to_find})
            set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
            find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
                      NAMES ${starpu_hdr}
                      HINTS ${STARPU_DIR}
                      PATH_SUFFIXES "include/starpu/${STARPU_FIND_VERSION}")
        endforeach()
    else()
        foreach(starpu_hdr ${STARPU_hdrs_to_find})
            set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
            if(STARPU_SHM_PKG_FILE_FOUND AND PC_STARPU_SHM_INCLUDE_DIRS)
                find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
                        NAMES ${starpu_hdr}
                        HINTS ${PC_STARPU_SHM_INCLUDE_DIRS})
            else()            
                find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
                        NAMES ${starpu_hdr}
                        PATHS ${PATH_TO_LOOK_FOR})
            endif()
        endforeach()
    endif()
endif()
mark_as_advanced(STARPU_starpu_config.h_INCLUDE_DIRS)

# Print status if not found
# -------------------------
if (NOT STARPU_starpu_config.h_INCLUDE_DIRS)
    Print_Find_Header_Status(starpu starpu_config.h)
endif ()


###
#
# GET_VERSION: Get the version of the software by parsing a file
#
###
MACRO(GET_VERSION _PACKAGE _filepath)

    #message(STATUS "Looking for ${_PACKAGE} version in the file ${_filepath}")
    file(READ "${_filepath}" _file)
    string(REGEX REPLACE
        "(.*)define([ \t]*)${_PACKAGE}_MAJOR_VERSION([ \t]*)([0-9]+)(.*)"
        "\\4" ${_PACKAGE}_VERSION_MAJOR "${_file}")
    string(REGEX REPLACE
        "(.*)define([ \t]*)${_PACKAGE}_MINOR_VERSION([ \t]*)([0-9]+)(.*)"
        "\\4" ${_PACKAGE}_VERSION_MINOR "${_file}")
    set(${_PACKAGE}_VERSION_STRING
        "${${_PACKAGE}_VERSION_MAJOR}.${${_PACKAGE}_VERSION_MINOR}")
    #message(STATUS "${_PACKAGE}_VERSION_MAJOR = -${${_PACKAGE}_VERSION_MAJOR}-")
    #message(STATUS "${_PACKAGE}_VERSION_MINOR = -${${_PACKAGE}_VERSION_MINOR}-")

ENDMACRO(GET_VERSION)

# Find the version of StarPU in starpu_config.h file
# remark: the version is defined in this file since the STARPU 1.0 version
if (STARPU_starpu_config.h_INCLUDE_DIRS)
    GET_VERSION("STARPU" "${STARPU_starpu_config.h_INCLUDE_DIRS}/starpu_config.h")
    if (STARPU_FIND_VERSION)
        if (STARPU_FIND_VERSION_EXACT STREQUAL 1)
            if( NOT (STARPU_FIND_VERSION_MAJOR STREQUAL STARPU_VERSION_MAJOR) OR
                NOT (STARPU_FIND_VERSION_MINOR STREQUAL STARPU_VERSION_MINOR) )
                message(FATAL_ERROR
                        "STARPU version found is ${STARPU_VERSION_STRING} "
                        "when required is ${STARPU_FIND_VERSION}")
            endif()
        else()
            # if the version found is older than the required then error
            if( (STARPU_FIND_VERSION_MAJOR STRGREATER STARPU_VERSION_MAJOR) OR
                (STARPU_FIND_VERSION_MINOR STRGREATER STARPU_VERSION_MINOR) )
                message(FATAL_ERROR
                        "STARPU version found is ${STARPU_VERSION_STRING} "
                        "when required is ${STARPU_FIND_VERSION} or newer")
            endif()
        endif()
    endif()
endif()


# Try to find the starpu headers in the given paths
# -------------------------------------------------

# create list of headers to find
set(STARPU_hdrs_to_find "starpu.h;starpu_profiling.h")
if(STARPU_LOOK_FOR_MPI)
    list(APPEND STARPU_hdrs_to_find "starpu_mpi.h")
endif()
if(STARPU_LOOK_FOR_CUDA)
    list(APPEND STARPU_hdrs_to_find "starpu_cuda.h;starpu_scheduler.h")
endif()

# call cmake macro to find the header path
if(DEFINED STARPU_INCDIR)
    foreach(starpu_hdr ${STARPU_hdrs_to_find})
        set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
        find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
                  NAMES ${starpu_hdr}
                  HINTS ${STARPU_INCDIR})                
    endforeach()
else()
    if(DEFINED STARPU_DIR)
        set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
        foreach(starpu_hdr ${STARPU_hdrs_to_find})
            find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
                      NAMES ${starpu_hdr}
                      HINTS ${STARPU_DIR}
                      PATH_SUFFIXES "include/starpu/${STARPU_VERSION_STRING}")
        endforeach()
    else()
        foreach(starpu_hdr ${STARPU_hdrs_to_find})
            set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
            
            if(STARPU_SHM_PKG_FILE_FOUND AND PC_STARPU_SHM_INCLUDE_DIRS)
                find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
                        NAMES ${starpu_hdr}
                        HINTS ${PC_STARPU_SHM_INCLUDE_DIRS})
            else()            
                find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
                        NAMES ${starpu_hdr}
                        PATHS ${PATH_TO_LOOK_FOR})
            endif()
        endforeach()
    endif()
endif()

# Print status if not found
# -------------------------
foreach(starpu_hdr ${STARPU_hdrs_to_find})
    if (NOT STARPU_${starpu_hdr}_INCLUDE_DIRS)
        Print_Find_Header_Status(starpu ${starpu_hdr})
    endif ()
endforeach()

# If found, add path to cmake variable
# ------------------------------------
set(STARPU_INCLUDE_DIRS "")
if(STARPU_SHM_PKG_FILE_FOUND)
    set(STARPU_INCLUDE_DIRS "${PC_STARPU_SHM_INCLUDE_DIRS}")
endif()
if(STARPU_MPI_PKG_FILE_FOUND)
    if(STARPU_INCLUDE_DIRS)
        list(APPEND STARPU_INCLUDE_DIRS "${PC_STARPU_MPI_INCLUDE_DIRS}")
    else()
        set(STARPU_INCLUDE_DIRS "${PC_STARPU_MPI_INCLUDE_DIRS}")
    endif()
endif()

foreach(starpu_hdr ${STARPU_hdrs_to_find})

    if (STARPU_${starpu_hdr}_INCLUDE_DIRS)
        # set cmake variables using the pkg-config naming convention
        list(APPEND STARPU_INCLUDE_DIRS "${STARPU_${starpu_hdr}_INCLUDE_DIRS}" )
    else ()
        message(STATUS "Looking for starpu -- ${starpu_hdr} not found")
        if(starpu_hdr STREQUAL "starpu_mpi.h")
            if(NOT ${STARPU_FIND_REQUIRED_MPI} STREQUAL 1)
                message(STATUS "Looking for starpu -- ${starpu_hdr} not required")
            else()
                list(APPEND STARPU_INCLUDE_DIRS "${STARPU_${starpu_hdr}_INCLUDE_DIRS}" )
            endif()
        elseif( (starpu_hdr STREQUAL "starpu_cuda.h") OR (starpu_hdr STREQUAL "starpu_scheduler.h") )
            if(NOT ${STARPU_FIND_REQUIRED_CUDA} STREQUAL 1)
                message(STATUS "Looking for starpu -- ${starpu_hdr} not required")
            else()
                list(APPEND STARPU_INCLUDE_DIRS "${STARPU_${starpu_hdr}_INCLUDE_DIRS}" )
            endif()
        endif()
    endif ()
    mark_as_advanced(STARPU_${starpu_hdr}_INCLUDE_DIRS)

endforeach(starpu_hdr ${STARPU_hdrs_to_find})

if (STARPU_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES STARPU_INCLUDE_DIRS)
endif ()


# Looking for lib
# ---------------

set(STARPU_SHM_LIBRARIES "")
set(STARPU_MPI_LIBRARIES "")
set(STARPU_LIBRARY_DIRS "")
if(STARPU_SHM_PKG_FILE_FOUND)
    set(STARPU_SHM_LIBRARIES "${PC_STARPU_SHM_LIBRARIES}")
    set(STARPU_LIBRARY_DIRS "${PC_STARPU_SHM_LIBRARY_DIRS}")
endif()
if(STARPU_MPI_PKG_FILE_FOUND)
    set(STARPU_MPI_LIBRARIES "${PC_STARPU_MPI_LIBRARIES}")
    if(STARPU_LIBRARY_DIRS)
        list(APPEND STARPU_LIBRARY_DIRS "${PC_STARPU_MPI_LIBRARY_DIRS}")
    else()
        set(STARPU_LIBRARY_DIRS "${PC_STARPU_MPI_LIBRARY_DIRS}")
    endif()
endif()

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

# Try to find the starpu libs in the given paths
# ----------------------------------------------

# create list of libs to find
set(STARPU_libs_to_find     "starpu-${STARPU_VERSION_STRING}")
set(STARPU_SHM_libs_to_find "starpu-${STARPU_VERSION_STRING}")
if (STARPU_LOOK_FOR_MPI)
    list(APPEND STARPU_libs_to_find "starpumpi-${STARPU_VERSION_STRING}")
    set(STARPU_MPI_libs_to_find "${STARPU_libs_to_find}")
endif()

# call cmake macro to find the lib path
if(DEFINED STARPU_LIBDIR)
    foreach(starpu_lib ${STARPU_libs_to_find})
        set(STARPU_${starpu_lib}_LIBRARY "STARPU_${starpu_lib}_LIBRARY-NOTFOUND")
        find_library(STARPU_${starpu_lib}_LIBRARY
                     NAMES ${starpu_lib}
                     HINTS ${STARPU_LIBDIR})
    endforeach()
else()
    if(DEFINED STARPU_DIR)
        foreach(starpu_lib ${STARPU_libs_to_find})
            set(STARPU_${starpu_lib}_LIBRARY "STARPU_${starpu_lib}_LIBRARY-NOTFOUND")
            find_library(STARPU_${starpu_lib}_LIBRARY
                         NAMES ${starpu_lib}
                         HINTS ${STARPU_DIR}
                         PATH_SUFFIXES lib lib32 lib64)
        endforeach()
    else()
        foreach(starpu_lib ${STARPU_libs_to_find})
            set(STARPU_${starpu_lib}_LIBRARY "STARPU_${starpu_lib}_LIBRARY-NOTFOUND")
            if(STARPU_SHM_PKG_FILE_FOUND AND PC_STARPU_SHM_LIBRARY_DIRS)
                find_library(STARPU_${starpu_lib}_LIBRARY
                            NAMES ${starpu_lib}
                            HINTS ${PC_STARPU_SHM_LIBRARY_DIRS})
            else()            
                find_library(STARPU_${starpu_lib}_LIBRARY
                            NAMES ${starpu_lib}
                            PATHS ${PATH_TO_LOOK_FOR})
            endif()
        endforeach()
    endif()
endif()

# Print status if not found
# -------------------------
foreach(starpu_lib ${STARPU_libs_to_find})
    if (NOT STARPU_${starpu_lib}_LIBRARY)
        Print_Find_Library_Status(starpu ${starpu_lib})
    endif ()
endforeach()

# If found, add path to cmake variable
# ------------------------------------
foreach(starpu_lib ${STARPU_libs_to_find})

    if (STARPU_${starpu_lib}_LIBRARY)
    
        get_filename_component(${starpu_lib}_lib_path ${STARPU_${starpu_lib}_LIBRARY} PATH)
        # set cmake variables (respects naming convention)
        
        foreach(starpu_shm_lib ${STARPU_SHM_libs_to_find})
            if(starpu_shm_lib STREQUAL starpu_lib)
                if (STARPU_SHM_LIBRARIES)
                    list(APPEND STARPU_SHM_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
                else()
                    set(STARPU_SHM_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
                endif()
            endif()
        endforeach()
        if (STARPU_LOOK_FOR_MPI)
            foreach(starpu_mpi_lib ${STARPU_MPI_libs_to_find})
                if(starpu_mpi_lib STREQUAL starpu_lib)
                    if (STARPU_MPI_LIBRARIES)
                        list(APPEND STARPU_MPI_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
                    else()
                        set(STARPU_MPI_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
                    endif()
                endif()
            endforeach()
        endif ()
        if (STARPU_LIBRARY_DIRS)
            list(APPEND STARPU_LIBRARY_DIRS "${${starpu_lib}_lib_path}")
        else()
            set(STARPU_LIBRARY_DIRS "${${starpu_lib}_lib_path}")
        endif()
        
    else (STARPU_${starpu_lib}_LIBRARY)
    
        message(STATUS "Looking for starpu -- lib ${starpu_lib} not found")
        if(starpu_lib STREQUAL "starpumpi-${STARPU_VERSION_STRING}" AND
           NOT ${STARPU_FIND_REQUIRED_MPI} STREQUAL 1)
            # if MPI optional, not a problem: no NOTFOUND in list of MPI LIBRARIES
            message(STATUS "Looking for starpu -- lib ${starpu_lib} not required")
        else()
            # for any other lib, add NOTFOUND in the proper list of LIBRARIES
            foreach(starpu_shm_lib ${STARPU_SHM_libs_to_find})
                if(starpu_shm_lib STREQUAL starpu_lib)
                    if (STARPU_SHM_LIBRARIES)
                        list(APPEND STARPU_SHM_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
                    else()
                        set(STARPU_SHM_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
                    endif()
                endif()
            endforeach()
            if (STARPU_LOOK_FOR_MPI)
                foreach(starpu_mpi_lib ${STARPU_MPI_libs_to_find})
                    if(starpu_mpi_lib STREQUAL starpu_lib)
                        if (STARPU_MPI_LIBRARIES)
                            list(APPEND STARPU_MPI_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
                        else()
                            set(STARPU_MPI_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
                        endif()
                    endif()
                endforeach()
            endif ()
        endif()
        
    endif (STARPU_${starpu_lib}_LIBRARY)
    
    mark_as_advanced(STARPU_${starpu_lib}_LIBRARY)
    
endforeach(starpu_lib ${STARPU_libs_to_find})

if (STARPU_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES STARPU_SHM_LIBRARIES)
    list(REMOVE_DUPLICATES STARPU_MPI_LIBRARIES)
    list(REMOVE_DUPLICATES STARPU_LIBRARY_DIRS)
endif ()

#message(STATUS "${STARPU_LIBRARY_list}")
#message(STATUS "${STARPU_INCLUDE_DIRS}")
#message(STATUS "${STARPU_LIBRARY_DIRS}")
#message(STATUS "${STARPU_SHM_LIBRARIES}")
#message(STATUS "${STARPU_MPI_LIBRARIES}")

# check that STARPU has been found
# --------------------------------
include(FindPackageHandleStandardArgs)
message(STATUS "StarPU has been found.")
if(STARPU_LOOK_FOR_MPI)
    message(STATUS "If the mpi version of StarPU has been found then we manage"
                   "two lists of libs, one sequential and one parallel (see"
                   "STARPU_SHM_LIBRARIES and STARPU_MPI_LIBRARIES).")
endif()
message(STATUS "StarPU shared memory libraries stored in STARPU_SHM_LIBRARIES")
find_package_handle_standard_args(STARPU DEFAULT_MSG
#                                  STARPU_FOUND
                                  STARPU_SHM_LIBRARIES
                                  STARPU_INCLUDE_DIRS
                                  STARPU_LIBRARY_DIRS)
if(STARPU_LOOK_FOR_MPI)
    message(STATUS "StarPU mpi libraries stored in STARPU_MPI_LIBRARIES")
    find_package_handle_standard_args(STARPU DEFAULT_MSG
                                      STARPU_MPI_LIBRARIES)
endif()

#
# TODO: Add possibility to check for specific functions in the library
#
