###
#
# @file CMakeLists.txt
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @author Asim YarKhan
# @date 2010-11-15
#
###

project( PLASMA )
cmake_minimum_required( VERSION 2.6 )
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

set( PLASMA_VERSION_MAJOR "2" )
set( PLASMA_VERSION_MINOR "4" )
set( PLASMA_VERSION_PATCH "6" )
set( PLASMA_VERSION "${PLASMA_VERSION_MAJOR}.${PLASMA_VERSION_MINOR}.${PLASMA_VERSION_PATCH}" )

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Disallow in-source build
string( COMPARE EQUAL "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}" INSOURCE )
if( INSOURCE )
  message( FATAL_ERROR "\n"
  "PLASMA requires an out-of-source-directory "
  "build. Please create a separate build directory and "
  "use CMake to configure and build there. \n"
  "e.g. mkdir build_dir; cd build_dir; cmake ..; make" )
endif()

# Release build is the default
if( NOT DEFINED CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
endif()

# Enable language compilers used in this project
enable_language( Fortran C )

# Set default install prefix to root of the build directory
#if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   set( CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/INSTALL" )
#endif ()

## Are we building 32 or 64 bit libraries
if (CMAKE_SIZEOF_VOID_P MATCHES 4)
  set( BUILD32BIT 1)
else ()
  set( BUILD64BIT 1)
endif ()

# --------------------------------------------------
# WINDOWS specific setup
if (WIN32 AND NOT UNIX)

  ## Setup Windows output locations
  set( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib/" CACHE PATH "Output directory for all libraries" )
  ## set( EXECUTABLE_OUTPUT_PATH "${LIBRARY_OUTPUT_PATH}/bin" CACHE PATH "Output directory for all executables")

  # Setup BLAS libaries, specifically MKL blas
  if (BUILD64BIT)     # 64 bit MKL BLAS
    set( CMAKE_FORCE_WIN64 )
    set( BLAS_LIBS mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib )
  else ()     # 32 bit MKL BLAS
    set( BLAS_LIBS mkl_intel_c.lib mkl_sequential.lib mkl_core.lib )
  endif ()

  set( BLAS_LIBRARIES "${BLAS_LIBS}" CACHE PATH "BLAS libraries" )
  set( CBLAS_LIBRARIES $ENV{LIBCBLAS} CACHE PATH "CBLAS BLAS C wrapper libraries" )
  set( LAPACKE_LIBRARIES $ENV{LIBLAPACKE} CACHE FILEPATH "LAPACKE cwrapper libraries" )
  set( LAPACKE_INCLUDES $ENV{INCLAPACKE} CACHE PATH "LAPACKE cwrapper include directory" )
  set( LAPACK_LIBRARIES $ENV{LIBLAPACK} CACHE PATH "LAPACK libraries" )
  set( LAPACK_INCLUDES $ENV{INCLAPACK} CACHE PATH "LAPACK include directory" )

  set( LINK_LIBS plasma coreblas quark ${LAPACKE_LIBRARIES} ${CBLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} libcmt.lib ${LINK_LIBS} )
  set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE /INCREMENTAL:NO" )

  set( INCLUDES ${INCLUDES} ${CLAPACK_INCLUDES} ${LAPACK_INCLUDES} )

  # Change /MD flags to /MT flags to get multi-threaded static libraries
  # Change any /O1 optimization flags to /O2
  foreach(flag_var
    CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
    CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_Fortran_FLAGS CMAKE_Fortran_FLAGS_DEBUG CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_MINSIZEREL CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
    if(${flag_var} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
    endif(${flag_var} MATCHES "/MD")
    if(${flag_var} MATCHES "/O1")
      string(REGEX REPLACE "/O1" "/O2" ${flag_var} "${${flag_var}}")
    endif(${flag_var} MATCHES "/O1")
  endforeach(flag_var)

  # Don't include MSVCRT or other dynamic runtime libraries by default
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:MSVCRT.lib /NODEFAULTLIB:libifcoremd.lib /NODEFAULTLIB:libmmd.lib")

  # Set C compiler specific flags
  if ( "${CMAKE_C_COMPILER}" MATCHES "icl" )
    # Intel C; Force C99 mode
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qstd=c99 /Zi /D_CRT_SECURE_NO_WARNINGS /Qdiag-disable:51 /MT" )
  elseif ( "${CMAKE_C_COMPILER}" MATCHES "cl" )
    # MS C compiler; Force C++ mode, error handler
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /TP /EHsc /MT" )
  else ()
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MT" )
  endif ()

  # Set Windows Fortran compiler specific flags
  if ( "${CMAKE_Fortran_COMPILER}" MATCHES "ifort" )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /fltconsistency /Zi /names:uppercase /MT" )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /DUPCASE" )
    set( Fortran_LINKER_C_MAIN_FLAGS " " )
  else ()
    set( Fortran_LINKER_C_MAIN_FLAGS " " )
  endif ()

# --------------------------------------------------
# Non Windows builds (e.g. Unix, Linux)

else() # Not Windows build

  set( BLAS_LIBRARIES "$ENV{LIBBLAS}" CACHE PATH "BLAS libraries" )
  set( CBLAS_LIBRARIES $ENV{LIBCBLAS} CACHE PATH "CBLAS BLAS C wrapper libraries" )
  set( LAPACKE_LIBRARIES $ENV{LIBLAPACKE} CACHE PATH "LAPACKE C wrapper libraries" )
  set( LAPACKE_INCLUDES $ENV{INCLAPACKE} CACHE PATH "LAPACK C wrapper include directory" )
  set( LAPACK_LIBRARIES $ENV{LIBLAPACK} CACHE PATH "LAPACK libraries" )
  set( LAPACK_INCLUDES $ENV{INCLAPACK} CACHE PATH "LAPACK include directory" )
  set( PTHREADS_LIBRARY $ENV{LIBPTHREADS} CACHE PATH "PTHREADS libraries" )

  set( INCLUDES ${INCLUDES} ${CLAPACK_INCLUDES} ${LAPACK_INCLUDES} )

  find_library( PTHREADS NAMES pthread pthreads DOC "PTHREADS library" REQUIRED)

#  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUPCASE" )
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DADD_" )

  set( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib/" CACHE PATH "Output directory for all libraries" )
  #set( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/" CACHE PATH "Output directory for all executables" )
  set( LINK_LIBS plasma coreblas quark ${LAPACKE_LIBRARIES} ${CBLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${PTHREADS_LIBRARY} m )

  if ( "${CMAKE_Fortran_COMPILER}" MATCHES "ifort" )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fltconsistency -fp_port " )
    set( Fortran_LINKER_C_MAIN_FLAGS " -nofor_main " )
  elseif ( "${CMAKE_Fortran_COMPILER}" MATCHES "ftn" )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  " )
    set( Fortran_LINKER_C_MAIN_FLAGS " -Mnomain " )
  else()
    set( CMAKE_Fortran_FLAGS " " )
    set( Fortran_LINKER_C_MAIN_FLAGS " " )
  endif ()

endif()

# --------------------------------------------------
# Multi architecture macros

# Macro to setup compile flags for binaries with C main files
macro( PLASMA_SETUP_C_MAIN target )
  target_link_libraries ( ${target} ${LINK_LIBS} )
  set_target_properties( ${target} PROPERTIES LINKER_LANGUAGE Fortran)
  set_target_properties( ${target} PROPERTIES LINK_FLAGS ${Fortran_LINKER_C_MAIN_FLAGS} )
endmacro( PLASMA_SETUP_C_MAIN )


# Macro to setup compile flags for binaries with Fortran main files
macro( PLASMA_SETUP_Fortran_MAIN target )
  target_link_libraries ( ${target} ${LINK_LIBS} )
  set_target_properties( ${target} PROPERTIES LINKER_LANGUAGE Fortran)
endmacro( PLASMA_SETUP_Fortran_MAIN )

# --------------------------------------------------
# Subdirectories that need to be processed

include_directories( include control quark ${INCLUDES} )
subdirs( include quark compute core_blas examples testing timing )

# Install some additional source files in the distribution
INSTALL( FILES CMakeLists.txt LICENSE docs/asciidoc/README.txt docs/asciidoc/ReleaseNotes.txt docs/asciidoc/CMakeBuildNotes.txt DESTINATION . )
# The following would be included in a Unix release, but are not required for Windows
#INSTALL( FILES Makefile Makefile.gen Makefile.internal Makefile.tau make.inc.example DESTINATION . )

# Install all the files in the docs directory in the distribution
INSTALL( DIRECTORY docs/pdf DESTINATION ./docs PATTERN ".svn" EXCLUDE PATTERN "*~" EXCLUDE PATTERN "#*" EXCLUDE)
# The following may be included in a full documentation release
# INSTALL( DIRECTORY docs/html DESTINATION ./docs PATTERN ".svn" EXCLUDE PATTERN "*~" EXCLUDE PATTERN "#*" EXCLUDE PATTERN "*htmlbrowsing*" EXCLUDE)

# --------------------------------------------------
# CPACK Packaging

SET(CPACK_PACKAGE_NAME "PLASMA")
SET(CPACK_PACKAGE_VENDOR "Innovative Computing Lab - University of Tennessee")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "PLASMA ${PLASMA_VERSION} - Parallel Linear Algebra Software for Multicore Architectures")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/docs/asciidoc/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION_MAJOR ${PLASMA_VERSION_MAJOR} )
SET(CPACK_PACKAGE_VERSION_MINOR ${PLASMA_VERSION_MINOR} )
SET(CPACK_PACKAGE_VERSION_PATCH ${PLASMA_VERSION_PATCH} )
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "PLASMA")

# Windows specific packaging
IF (WIN32 AND NOT UNIX)
  SET( CPACK_GENERATOR "NSIS")
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  ##SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
  ##SET(CPACK_NSIS_INSTALLED_ICON_NAME "PLASMA Linear Algebra Library")
  ##SET(CPACK_NSIS_MODIFY_PATH ON)
  SET(CPACK_NSIS_DISPLAY_NAME "PLASMA-${PLASMA_VERSION}")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\icl.eecs.utk.edu\\\\plasma")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\icl.eecs.utk.edu\\\\plasma")
  SET(CPACK_NSIS_CONTACT "plasma@eecs.utk.edu")
  set(CPACK_PACKAGE_RELOCATABLE "true")
  ## FIND_PATH(FREDIST_PATH libifcoremt.lib )
  IF(BUILD64BIT)
    SET(CPACK_SYSTEM_NAME win64)
    SET(FREDIST_PATH "c:/Program Files (x86)/Intel/Compiler/11.1/067/lib/intel64" )
    SET(LOCAL_PATH "d:/plasmawin/intel64" )
  ENDIF()
  IF(BUILD32BIT)
    SET(CPACK_SYSTEM_NAME win32)
    SET(FREDIST_PATH "c:/Program Files (x86)/Intel/Compiler/11.1/067/lib/ia32" )
    SET(LOCAL_PATH "d:/plasmawin/ia32" )
  ENDIF()
  INSTALL( PROGRAMS
    ${FREDIST_PATH}/ifconsol.lib
    ${FREDIST_PATH}/libifcoremt.lib
    ${FREDIST_PATH}/libifport.lib
    ${FREDIST_PATH}/libifportmd.lib
    ${FREDIST_PATH}/libirc.lib
    ${FREDIST_PATH}/libmmd.lib
    ${FREDIST_PATH}/libmmt.lib
    ${FREDIST_PATH}/svml_disp.lib
    ${FREDIST_PATH}/svml_dispmd.lib
    ${FREDIST_PATH}/svml_dispmt.lib
    ${FREDIST_PATH}/libdecimal.lib
    ${LOCAL_PATH}/lib/blas_cwrapper.lib
    ${LOCAL_PATH}/lib/lapacke.lib
    ${LOCAL_PATH}/lib/BLAS.lib
    ${LOCAL_PATH}/lib/LAPACK.lib
    ${LOCAL_PATH}/lib/MATGEN.lib
    DESTINATION lib )
  INCLUDE(InstallRequiredSystemLibraries)
  INSTALL( FILES
    # ${LOCAL_PATH}/include/cblas.h
    ${LOCAL_PATH}/include/lapacke.h
    ${LOCAL_PATH}/include/lapacke_config.h
    ${LOCAL_PATH}/include/lapacke_utils.h
    DESTINATION include )

ELSE()
  # Not windows
  SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_SOURCE_GENERATOR TGZ)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "plasma-${PLASMA_VERSION}" )
  SET(CPACK_SOURCE_IGNORE_FILES ~$ ./run* .svn ${CPACK_SOURCE_IGNORE_FILES} )
  MESSAGE(${CPACK_SOURCE_IGNORE_FILES} )
ENDIF()
INCLUDE(CPack)


message("------------------------------------------------------------------------")
message("Configuration:")
message("")
message(" Source code location: ${CMAKE_SOURCE_DIR}")
message(" Build type: ${CMAKE_BUILD_TYPE}")
message(" C Compiler: ${CMAKE_C_COMPILER}")
message(" C Compiler flags: ${CMAKE_C_FLAGS}")
message(" Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
message(" Fortran flags: ${CMAKE_Fortran_FLAGS}")
message(" Install path: ${CMAKE_INSTALL_PREFIX}")
message(" INCLUDES: ${INCLUDES}")
message(" LINK_LIBS: ${LINK_LIBS}")
message(" CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
message(" CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}" )
message("")
message("------------------------------------------------------------------------")
