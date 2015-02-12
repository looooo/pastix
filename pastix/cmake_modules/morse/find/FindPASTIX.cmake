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
# - Find PASTIX include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PASTIX
#               [REQUIRED]             # Fail with error if pastix is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )
#  Components available:
#   STARPU
#   CUDA (comes with cuBLAS): for use of GPUs
#   MPI: for use of multiple nodes of distributed memory
#   SCOTCH
#   PTSCOTCH
#   METIS
# This module finds headers and pastix library.
# Results are reported in variables:
#  PASTIX_FOUND           - True if headers and requested libraries were found
#  PASTIX_INCLUDE_DIRS    - pastix include directories
#  PASTIX_LIBRARY_DIRS    - Link directories for pastix libraries
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DPASTIX_DIR=path/to/pastix):
#  PASTIX_DIR             - Where to find the base directory of pastix
#  PASTIX_INCDIR          - Where to find the header files
#  PASTIX_LIBDIR          - Where to find the library files

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


if (NOT PASTIX_FOUND)
    set(PASTIX_DIR "" CACHE PATH "Root directory of PASTIX library")
    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "A cache variable, namely PASTIX_DIR, has been set to specify the install directory of PASTIX")
    endif()
endif()

# Try to find PASTIX dependencies if specified as COMPONENTS during the call
if( PASTIX_FIND_COMPONENTS )
    foreach( component ${PASTIX_FIND_COMPONENTS} )
        if(PASTIX_FIND_REQUIRED_${component})
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(PASTIX_${component}_FOUND TRUE)
            # should we have these variables available in gui modes?
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
            set(PASTIX_${component}_FOUND FALSE)
        endif()
    endforeach()
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE)

    pkg_search_module(PASTIX pastix)
    if (NOT PASTIX_FIND_QUIETLY)
        if (PASTIX_FOUND AND PASTIX_LIBRARIES)
            message(STATUS "Looking for PASTIX - found using PkgConfig")
            #if(NOT PASTIX_INCLUDE_DIRS)
            #    message("${Magenta}PASTIX_INCLUDE_DIRS is empty using PkgConfig."
            #        "Perhaps the path to pastix headers is already present in your"
            #        "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
            #endif()
        else()
            message("${Magenta}Looking for PASTIX - not found using PkgConfig."
                "Perhaps you should add the directory containing pastix.pc"
                "to the PKG_CONFIG_PATH environment variable.${ColourReset}")
        endif()
    endif()

    if (PASTIX_FIND_VERSION_EXACT)
        if( NOT (PASTIX_FIND_VERSION_MAJOR STREQUAL PASTIX_VERSION_MAJOR) OR
            NOT (PASTIX_FIND_VERSION_MINOR STREQUAL PASTIX_VERSION_MINOR) )
            if(NOT PASTIX_FIND_QUIETLY)
                message(FATAL_ERROR
                        "PASTIX version found is ${PASTIX_VERSION_STRING}"
                        "when required is ${PASTIX_FIND_VERSION}")
            endif()
        endif()
    else()
        # if the version found is older than the required then error
        if( (PASTIX_FIND_VERSION_MAJOR STRGREATER PASTIX_VERSION_MAJOR) OR
            (PASTIX_FIND_VERSION_MINOR STRGREATER PASTIX_VERSION_MINOR) )
            if(NOT PASTIX_FIND_QUIETLY)
                message(FATAL_ERROR
                        "PASTIX version found is ${PASTIX_VERSION_STRING}"
                        "when required is ${PASTIX_FIND_VERSION} or newer")
            endif()
        endif()
    endif()

endif(PKG_CONFIG_EXECUTABLE)

if(NOT PASTIX_FOUND OR NOT PASTIX_LIBRARIES)

    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - PkgConfig not used")
    endif()

    # Dependencies detection
    # ----------------------


    # Required dependencies
    # ---------------------

    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - Try to detect pthread")
    endif()
    if (PASTIX_FIND_REQUIRED)
        find_package(Threads REQUIRED)
    else()
        find_package(Threads)
    endif()
    set(PASTIX_EXTRA_LIBRARIES "")
    if( THREADS_FOUND )
        list(APPEND PASTIX_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
    endif ()

    # Add math library to the list of extra
    # it normally exists on all common systems provided with a C compiler
    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - Try to detect libm")
    endif()
    set(PASTIX_M_LIBRARIES "")
    if(UNIX OR WIN32)
        find_library(
            PASTIX_M_m_LIBRARY
            NAMES m
            )
        mark_as_advanced(PASTIX_M_m_LIBRARY)
        if (PASTIX_M_m_LIBRARY)
            list(APPEND PASTIX_M_LIBRARIES "${PASTIX_M_m_LIBRARY}")
            list(APPEND PASTIX_EXTRA_LIBRARIES "${PASTIX_M_m_LIBRARY}")
        else()
            if (PASTIX_FIND_REQUIRED)
                message(FATAL_ERROR "Could NOT find libm on your system."
                    "Are you sure to a have a C compiler installed?")
            endif()
        endif()
    endif()

    # Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
    # on Unix systems except Apple ones because it does not exist on it
    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - Try to detect librt")
    endif()
    set(PASTIX_RT_LIBRARIES "")
    if(UNIX AND NOT APPLE)
        find_library(
            PASTIX_RT_rt_LIBRARY
            NAMES rt
            )
        mark_as_advanced(PASTIX_RT_rt_LIBRARY)
        if (PASTIX_RT_rt_LIBRARY)
            list(APPEND PASTIX_RT_LIBRARIES "${PASTIX_RT_rt_LIBRARY}")
            list(APPEND PASTIX_EXTRA_LIBRARIES "${PASTIX_RT_rt_LIBRARY}")
        else()
            if (PASTIX_FIND_REQUIRED)
                message(FATAL_ERROR "Could NOT find librt on your system")
            endif()
        endif()
    endif()

    # PASTIX depends on BLAS
    #-----------------------
    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - Try to detect BLAS")
    endif()
    if (PASTIX_FIND_REQUIRED)
        find_package(BLASEXT REQUIRED)
    else()
        find_package(BLASEXT)
    endif()

    # PASTIX depends on HWLOC
    #------------------------
    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - Try to detect HWLOC")
    endif()
    if (PASTIX_FIND_REQUIRED)
        find_package(HWLOC REQUIRED)
    else()
        find_package(HWLOC)
    endif()

    # Optional dependencies
    # ---------------------

    # PASTIX may depend on CUDA/CUBLAS
    #---------------------------------
    if (NOT CUDA_FOUND)
        if (NOT PASTIX_FIND_QUIETLY)
            message(STATUS "Looking for PASTIX - Try to detect CUDA/cuBLAS")
        endif()
        if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_CUDA)
            find_package(CUDA REQUIRED)
        else()
            find_package(CUDA)
        endif()
    endif()

    # PASTIX may depend on MPI
    #-------------------------
    if (NOT MPI_FOUND)
        if (NOT PASTIX_FIND_QUIETLY)
            message(STATUS "Looking for PASTIX - Try to detect MPI")
        endif()
        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER mpicc)
        endif()
        if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_MPI)
            find_package(MPI REQUIRED)
        else()
            find_package(MPI)
        endif()
    endif (NOT MPI_FOUND)

    if( NOT STARPU_FOUND )

        if (NOT PASTIX_FIND_QUIETLY)
            message(STATUS "Looking for PASTIX - Try to detect StarPU")
        endif()

        set(PASTIX_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")

        # create list of components in order to make a single call to find_package(starpu...)
        # we explicitly need a StarPU version built with hwloc
        set(STARPU_COMPONENT_LIST "HWLOC")

        # StarPU may depend on MPI
        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER mpicc)
        endif()
        if (PASTIX_FIND_REQUIRED_MPI)
            if(PASTIX_FIND_REQUIRED_MPI)
                list(APPEND STARPU_COMPONENT_LIST "MPI")
            endif()
        endif()
        if (PASTIX_FIND_REQUIRED_CUDA)
            if(PASTIX_FIND_REQUIRED_CUDA)
                list(APPEND STARPU_COMPONENT_LIST "CUDA")
            endif()
        endif()
        # set the list of optional dependencies we may discover
        set(STARPU_OPTIONAL_COMPONENT_LIST "MPI" "CUDA" "MAGMA")
        if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_STARPU)
            find_package(STARPU ${PASTIX_STARPU_VERSION} REQUIRED
                         COMPONENTS ${STARPU_COMPONENT_LIST}
                         OPTIONAL_COMPONENTS ${STARPU_OPTIONAL_COMPONENT_LIST})
        else()
            find_package(STARPU ${PASTIX_STARPU_VERSION}
                         COMPONENTS ${STARPU_COMPONENT_LIST}
                         OPTIONAL_COMPONENTS ${STARPU_OPTIONAL_COMPONENT_LIST})
        endif()

    endif()

    # PASTIX may depends on SCOTCH
    #-----------------------------
    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - Try to detect SCOTCH")
    endif()
    if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_SCOTCH)
        find_package(SCOTCH REQUIRED)
    else()
        find_package(SCOTCH)
    endif()

    # PASTIX may depends on PTSCOTCH
    #-------------------------------
    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - Try to detect PTSCOTCH")
    endif()
    if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_PTSCOTCH)
        find_package(PTSCOTCH REQUIRED)
    else()
        find_package(PTSCOTCH)
    endif()

    # PASTIX may depends on METIS
    #----------------------------
    if (NOT PASTIX_FIND_QUIETLY)
        message(STATUS "Looking for PASTIX - Try to detect METIS")
    endif()
    if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_METIS)
        find_package(METIS REQUIRED)
    else()
        find_package(METIS)
    endif()

    # Error if pastix required and no partitioning lib found
    if (PASTIX_FIND_REQUIRED AND NOT SCOTCH_FOUND AND NOT PTSCOTCH_FOUND AND NOT METIS_FOUND)
        message(FATAL_ERROR "Could NOT find any partitioning library on your system"
        " (install scotch, ptscotch or metis)")
    endif()


    # Looking for PaStiX
    # ------------------

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


    # Try to find the pastix header in the given paths
    # ---------------------------------------------------
    # call cmake macro to find the header path
    if(PASTIX_INCDIR)
        set(PASTIX_pastix.h_DIRS "PASTIX_pastix.h_DIRS-NOTFOUND")
        find_path(PASTIX_pastix.h_DIRS
          NAMES pastix.h
          HINTS ${PASTIX_INCDIR})
    else()
        if(PASTIX_DIR)
            set(PASTIX_pastix.h_DIRS "PASTIX_pastix.h_DIRS-NOTFOUND")
            find_path(PASTIX_pastix.h_DIRS
              NAMES pastix.h
              HINTS ${PASTIX_DIR}
              PATH_SUFFIXES "include" "include/pastix")
        else()
            set(PASTIX_pastix.h_DIRS "PASTIX_pastix.h_DIRS-NOTFOUND")
            find_path(PASTIX_pastix.h_DIRS
              NAMES pastix.h
              HINTS ${_inc_env})
        endif()
    endif()
    mark_as_advanced(PASTIX_pastix.h_DIRS)

    # If found, add path to cmake variable
    # ------------------------------------
    if (PASTIX_pastix.h_DIRS)
        set(PASTIX_INCLUDE_DIRS "${PASTIX_pastix.h_DIRS}")
    else ()
        set(PASTIX_INCLUDE_DIRS "PASTIX_INCLUDE_DIRS-NOTFOUND")
        if(NOT PASTIX_FIND_QUIETLY)
            message(STATUS "Looking for pastix -- pastix.h not found")
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

    # Try to find the pastix lib in the given paths
    # ------------------------------------------------

    # create list of libs to find
    set(PASTIX_libs_to_find "pastix_murge;pastix")

    # call cmake macro to find the lib path
    if(PASTIX_LIBDIR)
        foreach(pastix_lib ${PASTIX_libs_to_find})
            set(PASTIX_${pastix_lib}_LIBRARY "PASTIX_${pastix_lib}_LIBRARY-NOTFOUND")
            find_library(PASTIX_${pastix_lib}_LIBRARY
                         NAMES ${pastix_lib}
                         HINTS ${PASTIX_LIBDIR})
        endforeach()
    else()
        if(PASTIX_DIR)
            foreach(pastix_lib ${PASTIX_libs_to_find})
                set(PASTIX_${pastix_lib}_LIBRARY "PASTIX_${pastix_lib}_LIBRARY-NOTFOUND")
                find_library(PASTIX_${pastix_lib}_LIBRARY
                             NAMES ${pastix_lib}
                             HINTS ${PASTIX_DIR}
                             PATH_SUFFIXES lib lib32 lib64)
            endforeach()
        else()
            foreach(pastix_lib ${PASTIX_libs_to_find})
                set(PASTIX_${pastix_lib}_LIBRARY "PASTIX_${pastix_lib}_LIBRARY-NOTFOUND")
                find_library(PASTIX_${pastix_lib}_LIBRARY
                             NAMES ${pastix_lib}
                             HINTS ${_lib_env})
            endforeach()
        endif()
    endif()

    # If found, add path to cmake variable
    # ------------------------------------
    foreach(pastix_lib ${PASTIX_libs_to_find})

        get_filename_component(${pastix_lib}_lib_path ${PASTIX_${pastix_lib}_LIBRARY} PATH)
        # set cmake variables (respects naming convention)
        if (PASTIX_LIBRARIES)
            list(APPEND PASTIX_LIBRARIES "${PASTIX_${pastix_lib}_LIBRARY}")
        else()
            set(PASTIX_LIBRARIES "${PASTIX_${pastix_lib}_LIBRARY}")
        endif()
        if (PASTIX_LIBRARY_DIRS)
            list(APPEND PASTIX_LIBRARY_DIRS "${${pastix_lib}_lib_path}")
        else()
            set(PASTIX_LIBRARY_DIRS "${${pastix_lib}_lib_path}")
        endif()
        mark_as_advanced(PASTIX_${pastix_lib}_LIBRARY)

    endforeach(pastix_lib ${PASTIX_libs_to_find})

    if(PASTIX_LIBRARIES)
        # check a function to validate the find
        if (PASTIX_INCLUDE_DIRS)
            set(CMAKE_REQUIRED_INCLUDES "${PASTIX_INCLUDE_DIRS}")
        endif()
        set(CMAKE_REQUIRED_FLAGS)
        foreach(libdir ${PASTIX_LIBRARY_DIRS})
            if (libdir)
                set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
            endif()
        endforeach()
        set(CMAKE_REQUIRED_LIBRARIES "${PASTIX_LIBRARIES}")
        if (STARPU_FOUND)
            if (STARPU_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${STARPU_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES}")
        endif()
        if (CUDA_FOUND)
            if (CUDA_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${CUDA_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${CUDA_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${CUDA_CUBLAS_LIBRARIES};${CUDA_LIBRARIES}")
        endif()
        if (MPI_FOUND)
            if (MPI_C_INCLUDE_PATH)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH}")
            endif()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${MPI_C_LIBRARIES}")
        endif()
        if (HWLOC_FOUND)
            if (HWLOC_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${HWLOC_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${HWLOC_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${HWLOC_LIBRARIES}")
        endif()
        if (BLAS_FOUND)
            if (BLAS_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${BLAS_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${BLAS_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${BLAS_LIBRARIES}")
        endif()
        if (SCOTCH_FOUND)
            if (SCOTCH_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${SCOTCH_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${SCOTCH_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${SCOTCH_LIBRARIES}")
        endif()
        if (PTSCOTCH_FOUND)
            if (PTSCOTCH_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${PTSCOTCH_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${PTSCOTCH_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${PTSCOTCH_LIBRARIES}")
        endif()
        if (METIS_FOUND)
            if (METIS_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${METIS_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${METIS_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${METIS_LIBRARIES}")
        endif()
        list(APPEND CMAKE_REQUIRED_LIBRARIES ${PASTIX_EXTRA_LIBRARIES})

        unset(PASTIX_WORKS CACHE)
        include(CheckFunctionExists)
        check_function_exists(pastix PASTIX_WORKS)
        mark_as_advanced(PASTIX_WORKS)

        if(PASTIX_WORKS)
            string(REPLACE " -L" ";" CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
            set(PASTIX_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
            set(PASTIX_LIBRARY_DIRS "${CMAKE_REQUIRED_FLAGS}")
            set(PASTIX_INCLUDE_DIRS "${CMAKE_REQUIRED_INCLUDES}")
        else()
            if(NOT PASTIX_FIND_QUIETLY)
                message(STATUS "Looking for PASTIX : test of pastix() fails")
                message(STATUS "PASTIX_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
                message(STATUS "PASTIX_LIBRARY_DIRS: ${CMAKE_REQUIRED_FLAGS}")
                message(STATUS "PASTIX_INCLUDE_DIRS: ${CMAKE_REQUIRED_INCLUDES}")
                message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
                message(STATUS "Looking for pastix : set PASTIX_LIBRARIES to NOTFOUND")
            endif()
            set(PASTIX_LIBRARIES "PASTIX_LIBRARIES-NOTFOUND")
        endif()
        set(CMAKE_REQUIRED_INCLUDES)
        set(CMAKE_REQUIRED_FLAGS)
        set(CMAKE_REQUIRED_LIBRARIES)
    endif(PASTIX_LIBRARIES)

endif(NOT PASTIX_FOUND OR NOT PASTIX_LIBRARIES)


# check that PASTIX has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PASTIX DEFAULT_MSG
                                  PASTIX_LIBRARIES)
