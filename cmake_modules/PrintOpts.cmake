###
#  @file PrintOpts.cmake
#
#  @copyright 2013-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.3
#  @author Tony Delarue
#  @date 2019-11-21
#
###

set(dep_message "\nConfiguration of Pastix:\n"
        "       PASTIX_VERSION ......: ${PASTIX_VERSION_MAJOR}.${PASTIX_VERSION_MINOR}.${PASTIX_VERSION_MICRO}\n"
        "       BUILDNAME ...........: ${BUILDNAME}\n"
        "       SITE ................: ${SITE}\n"
        "\n"
        "       Compiler: C .........: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})\n"
        "       Compiler: Fortran ...: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})\n")
if(PASTIX_USE_MPI)
  set(dep_message "${dep_message}"
  "       Compiler: MPI .......: ${MPI_C_COMPILER}\n"
  "       compiler flags ......: ${MPI_C_COMPILE_FLAGS}\n")
endif()
set(dep_message "${dep_message}"
"       Linker: .............: ${CMAKE_LINKER}\n"
"\n"
"       Build type ..........: ${CMAKE_BUILD_TYPE}\n"
"       Build shared ........: ${BUILD_SHARED_LIBS}\n"
"       CFlags ..............: ${CMAKE_C_FLAGS}\n"
"       LDFlags .............: ${CMAKE_SHARED_LINKER_FLAGS}\n"
"       EXE LDFlags .........: ${CMAKE_EXE_LINKER_FLAGS}\n"
"\n"
"       Implementation paradigm\n"
"       CUDA ................: ${PASTIX_WITH_CUDA}\n"
"       MPI .................: ${PASTIX_WITH_MPI}\n"
"       SCOTCH ..............: ${PASTIX_ORDERING_SCOTCH}\n"
"       PTSCOTCH ............: ${PASTIX_ORDERING_PTSCOTCH}\n"
"       METIS ...............: ${PASTIX_ORDERING_METIS}\n"

"\n"
"       Runtime specific\n"
"       PARSEC ..............: ${PASTIX_WITH_PARSEC}\n"
"       STARPU ..............: ${PASTIX_WITH_STARPU}\n"
"\n"
"       Kernels specific\n"
"       CBLAS ...............: ${CBLAS_FOUND}\n"
"       LAPACKE .............: ${LAPACKE_FOUND}\n"
"       HWLOC ...............: ${HWLOC_FOUND}\n"
"\n"
"       Blend Trace .........: ${PASTIX_BLEND_GENTRACE}\n"
"       Trace ...............: ${PASTIX_WITH_EZTRACE}\n"
"\n"
"       Binaries to build\n"
"       documentation .......: ${BUILD_DOCUMENTATION}\n"
"       testing .............: ${BUILD_TESTING}\n"
"       Precisions ..........: ${PASTIX_PRECISIONS}\n"
"\n")
# "       PASTIX dependencies :\n"
# foreach (_dep ${PASTIX_LIBRARIES_DEP})
#     set(dep_message "${dep_message}"
#     "                                 ${_dep}\n")
# endforeach ()
set(dep_message "${dep_message}"
"\n"
"       INSTALL_PREFIX ......: ${CMAKE_INSTALL_PREFIX}\n\n")

string(REPLACE ";" " " dep_message_wsc "${dep_message}")
message(${dep_message})
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config.log "${dep_message_wsc}")
# message(STATUS "Configuration is done - A summary of the current configuration"
# "\n   has been written in ${CMAKE_CURRENT_BINARY_DIR}/config.log")
# installation
# ------------
#INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.log DESTINATION share/pastix)
