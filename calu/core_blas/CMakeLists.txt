###
#
# @file CMakeLists.txt
#
# PLASMA is a software package provided by Univ. of Tennessee,
# Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.5.0
# @author Asim YarKhan
# @date 2010-11-15
#
###

INCLUDE( ${CMAKE_SOURCE_DIR}/GenerateZDCS.cmake OPTIONAL )

set( NAMES
        geqrt tsqrt tsmqr
        gelqt tslqt tsmlq
                      ttqrt ttmqr
                      ttlqt ttmlq
        getrf_incpiv tstrf gessm ssssm
        potrf trtri lauum
        gemm  symm  syrk trsm
        trmm  geadd  syr2k
        lacpy lange lansy
        ### core_dzasum unmqr unmlq
        ### hemm  herk  lanhe her2k
        ## lag2c
        swpab getrip shift
        plrnt plgsy ### plghe
        laswp getrf
        getrf_rectil getrf_reclap
        laset laset2
        ### herfb hegst
        larfx_tbrd trdalg brdalg
        hbelr hbrce hblrx
        gbelr gbrce gblrx
        ## tsmqr_hetra1
        tsmqr_corner
        ### tsmlq_hetra1
        tsmlq_corner
        parfb pamm pemv
)

foreach ( basename ${NAMES} )
  set( ZSRC ${ZSRC} core_z${basename}.c )
  set( CSRC ${CSRC} core_c${basename}.c )
  set( DSRC ${DSRC} core_d${basename}.c )
  set( SSRC ${SSRC} core_s${basename}.c )
endforeach()
set( ZSRC ${ZSRC} core_dzasum.c core_zunmqr.c core_zunmlq.c core_zherfb.c core_zhegst.c core_ztsmqr_hetra1.c core_ztsmlq_hetra1.c )
set( CSRC ${CSRC} core_scasum.c core_cunmqr.c core_cunmlq.c core_cherfb.c core_chegst.c core_ctsmqr_hetra1.c core_ctsmlq_hetra1.c )
set( DSRC ${DSRC} core_dasum.c core_dormqr.c core_dormlq.c core_dsyrfb.c core_dsygst.c core_dtsmqr_sytra1.c core_dtsmlq_sytra1.c )
set( SSRC ${SSRC} core_sasum.c core_sormqr.c core_sormlq.c core_ssyrfb.c core_ssygst.c core_stsmqr_sytra1.c core_stsmlq_sytra1.c )
set( ZSRC ${ZSRC} core_zhemm.c core_zherk.c core_zher2k.c core_zlanhe.c core_zplghe.c )
set( CSRC ${CSRC} core_chemm.c core_cherk.c core_cher2k.c core_clanhe.c core_cplghe.c )
set( ZSRC ${ZSRC} core_zlag2c.c )
set( DSRC ${DSRC} core_dlag2s.c )

# Introducing these dependencies causes the include files to be generated (if missing)
set( coreblas_SRC ${ZSRC} ${CSRC} ${DSRC} ${SSRC} core_alloc.c global.c )

ADD_LIBRARY( coreblas ${coreblas_SRC} )

# Custom target causes various header files to be generated
add_dependencies( coreblas include )

INSTALL( TARGETS coreblas ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin )

INSTALL( FILES ${coreblas_SRC} CMakeLists.txt DESTINATION core_blas )
