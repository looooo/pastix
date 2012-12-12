###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @author Asim YarKhan
# @date 2010-11-15
#
###

INCLUDE( ${CMAKE_SOURCE_DIR}/GenerateZDCS.cmake OPTIONAL )

set( timing_HEADERS
        auxiliary.h
        sauxiliary.h
        dauxiliary.h
        cauxiliary.h
        zauxiliary.h
        timing.h
        flops.h
)

# timing.c is included in other files via preprocessor #include commands
# timing.c

set( timing_NAMES
        gemm
        posv
        potrf
        gesv
        getrf
        gesv_incpiv
        getrf_incpiv
        gels
        geqrf
        gemm_tile
        posv_tile
        potrf_tile
        potri_tile
        gesv_tile
        getrf_tile
        gesv_incpiv_tile
        getrf_incpiv_tile
        gels_tile
        geqrf_tile
#        gecfi
#       heev_tile
#       hegv_tile
        gesvd_tile
#        lapack2tile
        getrf_reclap
        getrf_rectil
)

set( timing_AUXILIARY_NAMES
     gecfi
     lapack2tile
)

foreach( basename ${timing_NAMES} )
    ADD_EXECUTABLE( time_z${basename} time_z${basename}.c zauxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_z${basename} )
    ADD_EXECUTABLE( time_c${basename} time_c${basename}.c cauxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_c${basename} )
    ADD_EXECUTABLE( time_d${basename} time_d${basename}.c dauxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_d${basename} )
    ADD_EXECUTABLE( time_s${basename} time_s${basename}.c sauxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_s${basename} )
    SET ( timing_EXES ${timing_EXES} time_z${basename} time_c${basename} time_d${basename} time_s${basename} )
    SET ( timing_SRCS ${timing_SRCS} time_z${basename}.c time_c${basename}.c time_d${basename}.c time_s${basename}.c )
endforeach()

foreach( basename ${timing_AUXILIARY_NAMES} )
    ADD_EXECUTABLE( time_z${basename} time_z${basename}.c auxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_z${basename} )
    ADD_EXECUTABLE( time_c${basename} time_c${basename}.c auxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_c${basename} )
    ADD_EXECUTABLE( time_d${basename} time_d${basename}.c auxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_d${basename} )
    ADD_EXECUTABLE( time_s${basename} time_s${basename}.c auxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_s${basename} )
    SET ( timing_EXES ${timing_EXES} time_z${basename} time_c${basename} time_d${basename} time_s${basename} )
    SET ( timing_SRCS ${timing_SRCS} time_z${basename}.c time_c${basename}.c time_d${basename}.c time_s${basename}.c )
endforeach()

set( symm_NAMES
     ev_tile
     gv_tile
)

foreach( basename ${symm_NAMES} )
    ADD_EXECUTABLE( time_zhe${basename} time_zhe${basename}.c auxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_zhe${basename} )
    ADD_EXECUTABLE( time_che${basename} time_che${basename}.c auxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_che${basename} )
    ADD_EXECUTABLE( time_dsy${basename} time_dsy${basename}.c auxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_dsy${basename} )
    ADD_EXECUTABLE( time_ssy${basename} time_ssy${basename}.c auxiliary.c ${timing_HEADERS} )
    PLASMA_SETUP_C_MAIN( time_ssy${basename} )
    SET ( timing_EXES ${timing_EXES} time_zhe${basename} time_che${basename} time_dsy${basename} time_ssy${basename} )
    SET ( timing_SRCS ${timing_SRCS} time_zhe${basename}.c time_che${basename}.c time_dsy${basename}.c time_ssy${basename}.c )
endforeach()


SET ( timing_SRCS ${timing_SRCS} zauxiliary.c cauxiliary.c dauxiliary.c sauxiliary.c auxiliary.c timing.c auxiliary.c ${timing_HEADERS} )

ADD_CUSTOM_TARGET( timing ALL DEPENDS ${timing_EXES} )

INSTALL( FILES ${timing_SRCS} ${timing_HEADERS} CMakeLists.txt DESTINATION timing )
#INSTALL( TARGETS ${timing_EXES} DESTINATION timing )
INSTALL( TARGETS time_dgesv DESTINATION timing )
