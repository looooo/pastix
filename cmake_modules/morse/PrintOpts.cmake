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
#  @file PrintOpts.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver. 
#
#  @version 0.9.0
#  @author Florent Pruvost
#  @date 10-11-2014
#
###
string(TIMESTAMP the_date "%Y-%m-%d %H:%M:%S")
set(dep_message "\nConfiguration of Chameleon (${the_date}) :\n"
        "       BUILDNAME ...........: ${BUILDNAME}\n"
        "       SITE ................: ${SITE}\n"
        "\n"
        "       Compiler: C .........: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})\n"
        "       Compiler: Fortran ...: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})\n")
if(CHAMELEON_USE_MPI)
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
"       LDFlags .............: ${CMAKE_C_LINK_FLAGS}\n"
"\n"
"       Implementation paradigm\n"
"       CUDA ................: ${CHAMELEON_USE_CUDA}\n"
"       MPI .................: ${CHAMELEON_USE_MPI}\n"
"\n"
"       Runtime specific\n"
"       QUARK ...............: ${CHAMELEON_SCHED_QUARK}\n"
"       StarPU ..............: ${CHAMELEON_SCHED_STARPU}\n"
"       FxT .................: ${CHAMELEON_USE_FXT}\n"
"\n"
"       Kernels specific\n"
"       BLAS ................: ${BLA_VENDOR}\n"
"       MAGMA ...............: ${CHAMELEON_USE_MAGMA}\n"
"\n"
"       Simulation mode .....: ${CHAMELEON_SIMULATION}\n"
"\n"
"       Binaries to build\n"
"       documentation ........: ${CHAMELEON_ENABLE_DOCS}\n"
"       example ..............: ${CHAMELEON_ENABLE_EXAMPLE}\n"
"       testing ..............: ${CHAMELEON_ENABLE_TESTING}\n"
"       timing ...............: ${CHAMELEON_ENABLE_TIMING}\n"
"\n"
"       CHAMELEON dependencies :\n")
foreach (_dep ${CHAMELEON_DEP})
    set(dep_message "${dep_message}"
    "                                 ${_dep}\n")
endforeach ()
set(dep_message "${dep_message}"
"\n"
"       Definitions: ${CHAMELEON_DEFINITIONS_LIST}\n")
set(dep_message "${dep_message}"
"\n"
"       INSTALL_PREFIX ......: ${CMAKE_INSTALL_PREFIX}\n\n")

string(REPLACE ";" " " dep_message_wsc "${dep_message}")
#message(${dep_message})
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config.log "${dep_message_wsc}")
message(STATUS "Configuration is done - A summary of the current configuration"
"has been written in ${CMAKE_CURRENT_BINARY_DIR}/config.log")
# installation
# ------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.log DESTINATION share/chameleon)
