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
# - Find CHAMELEON include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(CHAMELEON
#               [REQUIRED]             # Fail with error if chameleon is not found
#               [COMPONENTS <libs>...] # required dependencies
#              )
#  Components available:
#   STARPU (default) or QUARK (not both can be activated): choose one of these runtimes
#   CUDA (comes with cuBLAS): for use of GPUs
#   MAGMA: for linear algebra on GPUs
#   MPI: for use of multiple nodes of distributed memory
# This module finds headers and chameleon library.
# Results are reported in variables:
#  CHAMELEON_FOUND           - True if headers and requested libraries were found
#  CHAMELEON_INCLUDE_DIRS    - chameleon include directories
#  CHAMELEON_LIBRARY_DIRS    - Link directories for chameleon libraries
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DCHAMELEON_DIR=path/to/chameleon):
#  CHAMELEON_DIR             - Where to find the base directory of chameleon
#  CHAMELEON_INCDIR          - Where to find the header files
#  CHAMELEON_LIBDIR          - Where to find the library files

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


if (NOT CHAMELEON_FOUND)
    set(CHAMELEON_DIR "" CACHE PATH "Root directory of CHAMELEON library")
    if (NOT CHAMELEON_FIND_QUIETLY)
        message(STATUS "A cache variable, namely CHAMELEON_DIR, has been set to specify the install directory of CHAMELEON")
    endif()
endif()

# Try to find CHAMELEON dependencies if specified as COMPONENTS during the call
if( CHAMELEON_FIND_COMPONENTS )
    foreach( component ${CHAMELEON_FIND_COMPONENTS} )
        if(CHAMELEON_FIND_REQUIRED_${component})
            find_package(${component} REQUIRED)
        else()
            find_package(${component})
        endif()
        if(${component}_FOUND)
            set(CHAMELEON_${component}_FOUND TRUE)
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
            set(CHAMELEON_${component}_FOUND FALSE)
        endif()
    endforeach()
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE)

    pkg_search_module(CHAMELEON chameleon)
    if (NOT CHAMELEON_FIND_QUIETLY)
        if (CHAMELEON_FOUND AND CHAMELEON_LIBRARIES)
            message(STATUS "Looking for CHAMELEON - found using PkgConfig")
            #if(NOT CHAMELEON_INCLUDE_DIRS)
            #    message("${Magenta}CHAMELEON_INCLUDE_DIRS is empty using PkgConfig."
            #        "Perhaps the path to chameleon headers is already present in your"
            #        "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
            #endif()
        else()
            message("${Magenta}Looking for CHAMELEON - not found using PkgConfig."
                "Perhaps you should add the directory containing chameleon.pc"
                "to the PKG_CONFIG_PATH environment variable.${ColourReset}")
        endif()
    endif()

    if (CHAMELEON_FIND_VERSION_EXACT)
        if( NOT (CHAMELEON_FIND_VERSION_MAJOR STREQUAL CHAMELEON_VERSION_MAJOR) OR
            NOT (CHAMELEON_FIND_VERSION_MINOR STREQUAL CHAMELEON_VERSION_MINOR) )
            if(NOT CHAMELEON_FIND_QUIETLY)
                message(FATAL_ERROR
                        "CHAMELEON version found is ${CHAMELEON_VERSION_STRING}"
                        "when required is ${CHAMELEON_FIND_VERSION}")
            endif()
        endif()
    else()
        # if the version found is older than the required then error
        if( (CHAMELEON_FIND_VERSION_MAJOR STRGREATER CHAMELEON_VERSION_MAJOR) OR
            (CHAMELEON_FIND_VERSION_MINOR STRGREATER CHAMELEON_VERSION_MINOR) )
            if(NOT CHAMELEON_FIND_QUIETLY)
                message(FATAL_ERROR
                        "CHAMELEON version found is ${CHAMELEON_VERSION_STRING}"
                        "when required is ${CHAMELEON_FIND_VERSION} or newer")
            endif()
        endif()
    endif()

endif(PKG_CONFIG_EXECUTABLE)

if(NOT CHAMELEON_FOUND OR NOT CHAMELEON_LIBRARIES)

    if (NOT CHAMELEON_FIND_QUIETLY)
        message(STATUS "Looking for CHAMELEON - PkgConfig not used")
    endif()

    # Dependencies detection
    # ----------------------

    if (NOT CHAMELEON_FIND_QUIETLY)
        message(STATUS "Looking for CHAMELEON - Try to detect pthread")
    endif()
    if (CHAMELEON_FIND_REQUIRED)
        find_package(Threads REQUIRED)
    else()
        find_package(Threads)
    endif()
    set(CHAMELEON_EXTRA_LIBRARIES "")
    if( THREADS_FOUND )
        list(APPEND CHAMELEON_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
    endif ()

    # Add math library to the list of extra
    # it normally exists on all common systems provided with a C compiler
    if (NOT CHAMELEON_FIND_QUIETLY)
        message(STATUS "Looking for CHAMELEON - Try to detect libm")
    endif()
    set(CHAMELEON_M_LIBRARIES "")
    if(UNIX OR WIN32)
        find_library(
            CHAMELEON_M_m_LIBRARY
            NAMES m
            )
        mark_as_advanced(CHAMELEON_M_m_LIBRARY)
        if (CHAMELEON_M_m_LIBRARY)
            list(APPEND CHAMELEON_M_LIBRARIES "${CHAMELEON_M_m_LIBRARY}")
            list(APPEND CHAMELEON_EXTRA_LIBRARIES "${CHAMELEON_M_m_LIBRARY}")
        else()
            if (CHAMELEON_FIND_REQUIRED)
                message(FATAL_ERROR "Could NOT find libm on your system."
                    "Are you sure to a have a C compiler installed?")
            endif()
        endif()
    endif()

    # Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
    # on Unix systems except Apple ones because it does not exist on it
    if (NOT CHAMELEON_FIND_QUIETLY)
        message(STATUS "Looking for CHAMELEON - Try to detect librt")
    endif()
    set(CHAMELEON_RT_LIBRARIES "")
    if(UNIX AND NOT APPLE)
        find_library(
            CHAMELEON_RT_rt_LIBRARY
            NAMES rt
            )
        mark_as_advanced(CHAMELEON_RT_rt_LIBRARY)
        if (CHAMELEON_RT_rt_LIBRARY)
            list(APPEND CHAMELEON_RT_LIBRARIES "${CHAMELEON_RT_rt_LIBRARY}")
            list(APPEND CHAMELEON_EXTRA_LIBRARIES "${CHAMELEON_RT_rt_LIBRARY}")
        else()
            if (CHAMELEON_FIND_REQUIRED)
                message(FATAL_ERROR "Could NOT find librt on your system")
            endif()
        endif()
    endif()

    # CHAMELEON depends on CBLAS
    #---------------------------
    if (NOT CHAMELEON_FIND_QUIETLY)
        message(STATUS "Looking for CHAMELEON - Try to detect CBLAS (depends on BLAS)")
    endif()
    if (CHAMELEON_FIND_REQUIRED)
        find_package(CBLAS REQUIRED COMPONENTS BLASEXT)
    else()
        find_package(CBLAS COMPONENTS BLASEXT)
    endif()

    # CHAMELEON depends on LAPACKE
    #-----------------------------

    # standalone version of lapacke seems useless for now
    # let the comment in case we meet some problems of non existing lapacke
    # functions in lapack library such as mkl, acml, ...
    #set(LAPACKE_STANDALONE TRUE)
    if (NOT CHAMELEON_FIND_QUIETLY)
        message(STATUS "Looking for CHAMELEON - Try to detect LAPACKE (depends on LAPACK)")
    endif()
    if (CHAMELEON_FIND_REQUIRED)
        find_package(LAPACKE REQUIRED COMPONENTS LAPACKEXT)
    else()
        find_package(LAPACKE COMPONENTS LAPACKEXT)
    endif()

    # CHAMELEON depends on TMG
    #-------------------------
    if (NOT CHAMELEON_FIND_QUIETLY)
        message(STATUS "Looking for CHAMELEON - Try to detect TMG (depends on LAPACK)")
    endif()
    if (CHAMELEON_FIND_REQUIRED)
        find_package(TMG REQUIRED)
    else()
        find_package(TMG)
    endif()

    # CHAMELEON may depend on CUDA/CUBLAS
    #------------------------------------
    if (NOT CUDA_FOUND)
        find_package(CUDA)
    endif()

    # CHAMELEON may depend on MAGMA gpu kernels
    # call our cmake module to test (in cmake_modules)
    # change this call position if not appropriated
    #-------------------------------------------------
    if ( CUDA_FOUND )
        set(CHAMELEON_MAGMA_VERSION "1.4" CACHE STRING "oldest MAGMA version desired")
        find_package(MAGMA ${CHAMELEON_MAGMA_VERSION} COMPONENTS CBLAS LAPACK CUDA)
    endif()

    # CHAMELEON depends on MPI
    #-------------------------
    if (NOT MPI_FOUND)

        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER mpicc)
        endif()
        find_package(MPI)

    endif (NOT MPI_FOUND)

    if( NOT STARPU_FOUND AND NOT QUARK_FOUND )

        set(CHAMELEON_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")

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
        if (CHAMELEON_FIND_REQUIRED_MPI)
            if(CHAMELEON_FIND_REQUIRED_MPI)
                list(APPEND STARPU_COMPONENT_LIST "MPI")
            endif()
        endif()
        if (CHAMELEON_FIND_REQUIRED_CUDA)
            if(CHAMELEON_FIND_REQUIRED_CUDA)
                list(APPEND STARPU_COMPONENT_LIST "CUDA")
            endif()
        endif()
        # set the list of optional dependencies we may discover
        set(STARPU_OPTIONAL_COMPONENT_LIST "MPI" "CUDA" "MAGMA")
        find_package(STARPU ${CHAMELEON_STARPU_VERSION}
                     COMPONENTS ${STARPU_COMPONENT_LIST}
                     OPTIONAL_COMPONENTS ${STARPU_OPTIONAL_COMPONENT_LIST})

    endif()

    if( NOT STARPU_FOUND AND NOT QUARK_FOUND )

        # try to find quark runtime
        find_package(QUARK COMPONENTS HWLOC)

        if (NOT QUARK_FOUND)
            message(FATAL_ERROR "Nor StarPU (default runtime) neither QUARK"
            "runtimes have been found while at least one of them should be installed")
        endif()

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


    # Try to find the chameleon header in the given paths
    # ---------------------------------------------------
    # call cmake macro to find the header path
    if(CHAMELEON_INCDIR)
        set(CHAMELEON_morse.h_DIRS "CHAMELEON_morse.h_DIRS-NOTFOUND")
        find_path(CHAMELEON_morse.h_DIRS
          NAMES morse.h
          HINTS ${CHAMELEON_INCDIR})
    else()
        if(CHAMELEON_DIR)
            set(CHAMELEON_morse.h_DIRS "CHAMELEON_morse.h_DIRS-NOTFOUND")
            find_path(CHAMELEON_morse.h_DIRS
              NAMES morse.h
              HINTS ${CHAMELEON_DIR}
              PATH_SUFFIXES "include" "include/chameleon")
        else()
            set(CHAMELEON_morse.h_DIRS "CHAMELEON_morse.h_DIRS-NOTFOUND")
            find_path(CHAMELEON_morse.h_DIRS
              NAMES morse.h
              HINTS ${_inc_env})
        endif()
    endif()
    mark_as_advanced(CHAMELEON_morse.h_DIRS)

    # If found, add path to cmake variable
    # ------------------------------------
    if (CHAMELEON_morse.h_DIRS)
        set(CHAMELEON_INCLUDE_DIRS "${CHAMELEON_morse.h_DIRS}")
    else ()
        set(CHAMELEON_INCLUDE_DIRS "CHAMELEON_INCLUDE_DIRS-NOTFOUND")
        if(NOT CHAMELEON_FIND_QUIETLY)
            message(STATUS "Looking for chameleon -- morse.h not found")
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

    # Try to find the chameleon lib in the given paths
    # ------------------------------------------------

    # create list of libs to find
    set(CHAMELEON_libs_to_find "chameleon")
    if (STARPU_FOUND)
        list(APPEND CHAMELEON_libs_to_find "chameleon_starpu")
    elseif (QUARK_FOUND)
        list(APPEND CHAMELEON_libs_to_find "chameleon_quark")
    endif()
    list(APPEND CHAMELEON_libs_to_find "coreblas")

    # call cmake macro to find the lib path
    if(CHAMELEON_LIBDIR)
        foreach(chameleon_lib ${CHAMELEON_libs_to_find})
            set(CHAMELEON_${chameleon_lib}_LIBRARY "CHAMELEON_${chameleon_lib}_LIBRARY-NOTFOUND")
            find_library(CHAMELEON_${chameleon_lib}_LIBRARY
                         NAMES ${chameleon_lib}
                         HINTS ${CHAMELEON_LIBDIR})
        endforeach()
    else()
        if(CHAMELEON_DIR)
            foreach(chameleon_lib ${CHAMELEON_libs_to_find})
                set(CHAMELEON_${chameleon_lib}_LIBRARY "CHAMELEON_${chameleon_lib}_LIBRARY-NOTFOUND")
                find_library(CHAMELEON_${chameleon_lib}_LIBRARY
                             NAMES ${chameleon_lib}
                             HINTS ${CHAMELEON_DIR}
                             PATH_SUFFIXES lib lib32 lib64)
            endforeach()
        else()
            foreach(chameleon_lib ${CHAMELEON_libs_to_find})
                set(CHAMELEON_${chameleon_lib}_LIBRARY "CHAMELEON_${chameleon_lib}_LIBRARY-NOTFOUND")
                find_library(CHAMELEON_${chameleon_lib}_LIBRARY
                             NAMES ${chameleon_lib}
                             HINTS ${_lib_env})
            endforeach()
        endif()
    endif()

    # If found, add path to cmake variable
    # ------------------------------------
    foreach(chameleon_lib ${CHAMELEON_libs_to_find})

        get_filename_component(${chameleon_lib}_lib_path ${CHAMELEON_${chameleon_lib}_LIBRARY} PATH)
        # set cmake variables (respects naming convention)
        if (CHAMELEON_LIBRARIES)
            list(APPEND CHAMELEON_LIBRARIES "${CHAMELEON_${chameleon_lib}_LIBRARY}")
        else()
            set(CHAMELEON_LIBRARIES "${CHAMELEON_${chameleon_lib}_LIBRARY}")
        endif()
        if (CHAMELEON_LIBRARY_DIRS)
            list(APPEND CHAMELEON_LIBRARY_DIRS "${${chameleon_lib}_lib_path}")
        else()
            set(CHAMELEON_LIBRARY_DIRS "${${chameleon_lib}_lib_path}")
        endif()
        mark_as_advanced(CHAMELEON_${chameleon_lib}_LIBRARY)

    endforeach(chameleon_lib ${CHAMELEON_libs_to_find})

    if(CHAMELEON_LIBRARIES)
        # check a function to validate the find
        if (CHAMELEON_INCLUDE_DIRS)
            set(CMAKE_REQUIRED_INCLUDES "${CHAMELEON_INCLUDE_DIRS}")
        endif()
        set(CMAKE_REQUIRED_FLAGS)
        foreach(libdir ${CHAMELEON_LIBRARY_DIRS})
            if (libdir)
                set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
            endif()
        endforeach()
        set(CMAKE_REQUIRED_LIBRARIES "${CHAMELEON_LIBRARIES}")
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
        if (QUARK_FOUND)
            if (QUARK_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "$QUARK_INCLUDE_DIRS}")
            endif()
            if(QUARK_LIBRARY_DIRS)
                set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${QUARK_LIBRARY_DIRS}")
            endif()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${QUARK_LIBRARIES}")
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
        if (MAGMA_FOUND)
            if (EXISTS MAGMA_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${MAGMA_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${MAGMA_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${${MAGMA_LIBRARIES}}")
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
        if (TMG_FOUND)
            if (TMG_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${TMG_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${TMG_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${TMG_LIBRARIES}")
        endif()
        if (LAPACKE_FOUND)
            if (LAPACKE_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${LAPACKE_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${LAPACKE_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${LAPACKE_LIBRARIES}")
        endif()
        if (CBLAS_FOUND)
            if (CBLAS_INCLUDE_DIRS)
                list(APPEND CMAKE_REQUIRED_INCLUDES "${CBLAS_INCLUDE_DIRS}")
            endif()
            foreach(libdir ${CBLAS_LIBRARY_DIRS})
                if (libdir)
                    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -L${libdir}")
                endif()
            endforeach()
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${CBLAS_LIBRARIES}")
        endif()
        list(APPEND CMAKE_REQUIRED_LIBRARIES ${CHAMELEON_EXTRA_LIBRARIES})

        unset(CHAMELEON_WORKS CACHE)
        include(CheckFunctionExists)
        check_function_exists(MORSE_Init CHAMELEON_WORKS)
        mark_as_advanced(CHAMELEON_WORKS)

        if(CHAMELEON_WORKS)
            string(REPLACE " -L" ";" CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
            set(CHAMELEON_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
            set(CHAMELEON_LIBRARY_DIRS "${CMAKE_REQUIRED_FLAGS}")
            set(CHAMELEON_INCLUDE_DIRS "${CMAKE_REQUIRED_INCLUDES}")
        else()
            if(NOT CHAMELEON_FIND_QUIETLY)
                message(STATUS "Looking for chameleon : test of MORSE_Init fails")
                message(STATUS "CHAMELEON_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
                message(STATUS "CHAMELEON_LIBRARY_DIRS: ${CMAKE_REQUIRED_FLAGS}")
                message(STATUS "CHAMELEON_INCLUDE_DIRS: ${CMAKE_REQUIRED_INCLUDES}")
                message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
                message(STATUS "Looking for chameleon : set CHAMELEON_LIBRARIES to NOTFOUND")
            endif()
            set(CHAMELEON_LIBRARIES "CHAMELEON_LIBRARIES-NOTFOUND")
        endif()
        set(CMAKE_REQUIRED_INCLUDES)
        set(CMAKE_REQUIRED_FLAGS)
        set(CMAKE_REQUIRED_LIBRARIES)
    endif(CHAMELEON_LIBRARIES)

endif(NOT CHAMELEON_FOUND OR NOT CHAMELEON_LIBRARIES)


# check that CHAMELEON has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CHAMELEON DEFAULT_MSG
                                  CHAMELEON_LIBRARIES)
