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

add_subdirectory(lin)

INCLUDE( ${CMAKE_SOURCE_DIR}/GenerateZDCS.cmake OPTIONAL )

set( testing_NAMES
        main
        gesv
        getri
        gesvd
        gels
        posv
        potri
        gemm
#        heev
#        hegst
#        hegv
        syrk
        syr2k
        symm
        trsm
        trmm
        lange
##      herk
##      her2k
##      hemm
        gecfi
        getmi
##      cposv
##      cgesv
##      cungesv
        pemv
)
set( testing_MNAMES
     posv
     gesv
     ungesv
)

foreach( name ${testing_NAMES} )
  set( testing_ZSRC ${testing_ZSRC} testing_z${name}.c )
  set( testing_CSRC ${testing_CSRC} testing_c${name}.c )
  set( testing_DSRC ${testing_DSRC} testing_d${name}.c )
  set( testing_SSRC ${testing_SSRC} testing_s${name}.c )
endforeach()

foreach( mname ${testing_MNAMES} )
  set( testing_ZSRC ${testing_ZSRC} testing_zc${mname}.c )
  set( testing_DSRC ${testing_DSRC} testing_ds${mname}.c )
endforeach()

set( testing_ZSRC ${testing_ZSRC} testing_zherk.c testing_zher2k.c testing_zhemm.c testing_zmain.h testing_zheev.c testing_zhegst.c testing_zhegv.c )
set( testing_CSRC ${testing_CSRC} testing_cherk.c testing_cher2k.c testing_chemm.c testing_cmain.h testing_cheev.c testing_chegst.c testing_chegv.c )
set( testing_DSRC ${testing_DSRC} testing_dmain.h testing_dsyev.c testing_dsygst.c testing_dsygv.c )
set( testing_SSRC ${testing_SSRC} testing_smain.h testing_ssyev.c testing_ssygst.c testing_ssygv.c )

set( testing_SRCS ${testing_ZSRC} ${testing_CSRC} ${testing_DSRC} ${testing_SSRC} )

add_executable( ztesting ${testing_ZSRC} )
add_executable( ctesting ${testing_CSRC} )
add_executable( dtesting ${testing_DSRC} )
add_executable( stesting ${testing_SSRC} )

set( EXES ztesting ctesting dtesting stesting )

foreach( exe ${EXES} )
  PLASMA_SETUP_C_MAIN( ${exe} )
endforeach()

INSTALL( TARGETS ${EXES} DESTINATION testing )
INSTALL( PROGRAMS plasma_testing.py DESTINATION testing )

INSTALL( FILES ${testing_SRCS} CMakeLists.txt DESTINATION testing )
