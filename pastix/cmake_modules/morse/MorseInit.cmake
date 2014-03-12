###
#
# -- Inria
# -- (C) Copyright 2013
#
# This software is a computer program whose purpose is to process
# Matrices Over Runtime Systems @ Exascale (MORSE). More information
# can be found on the following website: http://www.inria.fr/en/teams/morse.
#
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
#
###
#
#  @file MorseInit.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 2.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###

list(APPEND CMAKE_MODULE_PATH ${MORSE_CMAKE_MODULE_PATH}/find)

include (CMakeDetermineSystem)
include (CheckCCompilerFlag)
include (CheckFunctionExists)
include (CheckSymbolExists)
include (CheckIncludeFiles)

#
# check the capabilities of the system we are building for
#

# Check we are using the reentrant version of xlc/xlf if used
STRING(REGEX MATCH ".*xlc$" _match_xlc ${CMAKE_C_COMPILER})
IF (_match_xlc)
     MESSAGE(ERROR "Please use the thread-safe version of the xlc compiler (xlc_r)")
ENDIF (_match_xlc)
STRING(REGEX MATCH ".*xlf$" _match_xlf ${CMAKE_Fortran_COMPILER})
IF (_match_xlf)
     MESSAGE(ERROR "Please use the thread-safe version of the xlf compiler (xlf_r)")
ENDIF (_match_xlf)

# Fortran tricks (BLAS library require to link with fortran linker, while main in in C)
STRING(REGEX MATCH "Intel" _match_intel ${CMAKE_Fortran_COMPILER_ID})
IF (_match_intel)
     MESSAGE(STATUS "Add -nofor_main to the Fortran linker.")
     SET(LOCAL_FORTRAN_LINK_FLAGS "${LOCAL_FORTRAN_LINK_FLAGS} -nofor_main")
ENDIF (_match_intel)

STRING(REGEX MATCH "PGI$" _match_pgi ${CMAKE_Fortran_COMPILER_ID})
IF (_match_pgi)
    MESSAGE(STATUS "Add -Mnomain to the Fortran linker.")
    SET(LOCAL_FORTRAN_LINK_FLAGS "${LOCAL_FORTRAN_LINK_FLAGS} -Mnomain -Bstatic")
ENDIF (_match_pgi)

STRING(REGEX MATCH "XL" _match_xlc ${CMAKE_C_COMPILER_ID})
IF (_match_xlc AND BUILD_64bits)
     MESSAGE(STATUS "Add -q64 to the C compiler/linker.")
     SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -q64")
ENDIF (_match_xlc AND BUILD_64bits)

STRING(REGEX MATCH "XL$" _match_xlf ${CMAKE_Fortran_COMPILER_ID})
IF (_match_xlf)
  SET(arch_flags "-q32")
  IF(BUILD_64bits)
     SET(arch_flags "-q64")
  ENDIF(BUILD_64bits)
  MESSAGE(STATUS "Add ${arch_flags} and -nofor_main to the Fortran linker.")
  SET(LOCAL_FORTRAN_LINK_FLAGS "${LOCAL_FORTRAN_LINK_FLAGS} ${arch_flags} -nofor_main")
ENDIF (_match_xlf)

#
# check for the CPU we build for
#
MESSAGE(STATUS "Building for target ${CMAKE_SYSTEM_PROCESSOR}")
STRING(REGEX MATCH "(i.86-*)|(athlon-*)|(pentium-*)" _mach_x86 ${CMAKE_SYSTEM_PROCESSOR})
IF (_mach_x86)
    MESSAGE(STATUS "Found target for X86")
    SET(ARCH_X86 1)
ENDIF (_mach_x86)

STRING(REGEX MATCH "(x86_64-*)|(X86_64-*)|(AMD64-*)|(amd64-*)" _mach_x86_64 ${CMAKE_SYSTEM_PROCESSOR})
IF (_mach_x86_64)
    MESSAGE(STATUS "Found target X86_64")
    SET(ARCH_X86_64 1)
ENDIF (_mach_x86_64)

STRING(REGEX MATCH "(ppc-*)|(powerpc-*)" _mach_ppc ${CMAKE_SYSTEM_PROCESSOR})
IF (_mach_ppc)
    MESSAGE(STATUS "Found target for PPC")
    SET(ARCH_PPC 1)
ENDIF (_mach_ppc)

#
# Fix the building system for 32 or 64 bits.
#
# On MAC OS X there is a easy solution, by setting the
# CMAKE_OSX_ARCHITECTURES to a subset of the following values:
# ppc;ppc64;i386;x86_64.
# On Linux this is a little bit tricky. We have to check that the
# compiler supports the -m32/-m64 flags as well as the linker.
# Once this issue resolved the CMAKE_C_FLAGS and CMAKE_C_LDFLAGS
# have to be updated accordingly.
#
# TODO: Same trick for the Fortran compiler...
#       no idea how to correctly detect if the required/optional
#          libraries are in the correct format.
#
set(SAVE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
if (BUILD_64bits)
  if( _match_xlc)
    set( ARCH_BUILD "-q64" )
  else (_match_xlc)
    set( ARCH_BUILD "-m64" )
  endif(_match_xlc)
else (BUILD_64bits)
  if( _match_xlc)
    set( ARCH_BUILD "-q32" )
  else (_match_xlc)
    set( ARCH_BUILD "-m32" )
  endif(_match_xlc)
endif (BUILD_64bits)

set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${ARCH_BUILD}")
check_c_compiler_flag(${ARCH_BUILD} C_M32or64)

if (C_M32or64)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ARCH_BUILD}")
  set(CMAKE_C_LDFLAGS "${CMAKE_C_LDFLAGS} ${ARCH_BUILD}")
  set(LOCAL_FORTRAN_LINK_FLAGS "${LOCAL_FORTRAN_LINK_FLAGS} ${ARCH_BUILD}")
else (C_M32or64)
  set(CMAKE_REQUIRED_FLAGS "${SAVE_CMAKE_REQUIRED_FLAGS}")
endif (C_M32or64)
unset( SAVE_CMAKE_REQUIRED_FLAGS )

#
# Check compiler flags and capabilities
#
IF( NOT _match_xlc )
  CHECK_C_COMPILER_FLAG( "-std=c99" HAVE_STD_C99)
  IF( HAVE_STD_C99 )
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99" )
  ENDIF( HAVE_STD_C99 )
ELSE( NOT _match_xlc )
  CHECK_C_COMPILER_FLAG( "-qlanglvl=extc99" HAVE_STD_C99)
  IF( HAVE_STD_C99 )
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qlanglvl=extc99" )
  ENDIF( HAVE_STD_C99 )
ENDIF( NOT _match_xlc )

# Set warnings for debug builds
CHECK_C_COMPILER_FLAG( "-Wall" HAVE_WALL )
IF( HAVE_WALL )
    SET( C_WFLAGS "${C_WFLAGS} -Wall" )
ENDIF( HAVE_WALL )
CHECK_C_COMPILER_FLAG( "-Wextra" HAVE_WEXTRA )
IF( HAVE_WEXTRA )
    SET( C_WFLAGS "${C_WFLAGS} -Wextra" )
ENDIF( HAVE_WEXTRA )
# flags for the overly verbose icc
CHECK_C_COMPILER_FLAG( "-wd424" HAVE_WD )
IF( HAVE_WD )
    # 424: checks for duplicate ";"
    # 981: every volatile triggers a "unspecified evaluation order", obnoxious
    #      but might be useful for some debugging sessions.
    # 1419: warning about extern functions being declared in .c
    #       files
    # 1572: cuda compares floats with 0.0f.
    SET( C_WFLAGS "${C_WFLAGS} -wd424 -wd981 -wd1419 -wd1572" )
ENDIF( HAVE_WD )
CHECK_C_COMPILER_FLAG( "-g3" HAVE_G3 )
IF( HAVE_G3 )
    SET( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g3 ${C_WFLAGS}" )
ELSE()
    SET( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g ${C_WFLAGS}" )
ENDIF( HAVE_G3)
SET( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${C_WFLAGS}" )
SET( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -DNDEBUG" )

if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set( MAC_OS_X 1 CACHE INTERNAL "Compile on MAC OS X")
endif(CMAKE_SYSTEM_NAME MATCHES "Darwin")

#
# Threads and Atomics
#

# Atomics
#include (cmake_modules/CheckAtomicIntrinsic.cmake)

# Threads
find_package(Threads)
if(Threads_FOUND)
  set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};${CMAKE_THREAD_LIBS_INIT}")
  check_function_exists(pthread_create HAVE_PTHREAD)
  if(HAVE_PTHREAD)
    list(APPEND EXTRA_LIBS ${CMAKE_THREAD_LIBS_INIT})
  endif(HAVE_PTHREAD)
endif(Threads_FOUND)

check_function_exists(sched_setaffinity HAVE_SCHED_SETAFFINITY)
if( NOT HAVE_SCHED_SETAFFINITY )
  check_library_exists(rt sched_setaffinity "" HAVE_SCHED_SETAFFINITY)
endif( NOT HAVE_SCHED_SETAFFINITY )

#
# timeval, timespec, realtime clocks, etc
#
include(CheckStructHasMember)
check_struct_has_member("struct timespec" tv_nsec time.h HAVE_TIMESPEC_TV_NSEC)
if( NOT HAVE_TIMESPEC_TV_NSEC )
  add_definitions(-D_GNU_SOURCE)
  check_struct_has_member("struct timespec" tv_nsec time.h HAVE_TIMESPEC_TV_NSEC)
endif( NOT HAVE_TIMESPEC_TV_NSEC )
check_library_exists(rt clock_gettime "" HAVE_CLOCK_GETTIME)
if( HAVE_CLOCK_GETTIME )
  list(APPEND EXTRA_LIBS rt)
endif( HAVE_CLOCK_GETTIME )

#
# stdlib, stdio, string, getopt, etc
#
check_include_files(stdarg.h HAVE_STDARG_H)
# va_copy is special as it is not required to be a function.
if (HAVE_STDARG_H)
  check_c_source_compiles("
      #include <stdarg.h>
      int main(void) {
          va_list a, b;
          va_copy(a, b);
          return 0;
      }"
      HAVE_VA_COPY
      )

  if (NOT HAVE_VA_COPY)
    check_c_source_compiles("
        #include <stdarg.h>
        int main(void) {
            va_list a, b;
            __va_copy(a, b);
            return 0;
        }"
        HAVE_UNDERSCORE_VA_COPY
        )
  endif (NOT HAVE_VA_COPY)
endif (HAVE_STDARG_H)
check_function_exists(asprintf HAVE_ASPRINTF)
check_function_exists(vasprintf HAVE_VASPRINTF)
check_include_files(getopt.h HAVE_GETOPT_H)
check_include_files(unistd.h HAVE_UNISTD_H)
check_function_exists(getopt_long HAVE_GETOPT_LONG)
check_include_files(errno.h HAVE_ERRNO_H)
check_include_files(stddef.h HAVE_STDDEF_H)
check_function_exists(getrusage HAVE_GETRUSAGE)
check_include_files(limits.h HAVE_LIMITS_H)
check_include_files(string.h HAVE_STRING_H)
check_include_files(complex.h HAVE_COMPLEX_H)

##
## @end file MorseInit.cmake
##
