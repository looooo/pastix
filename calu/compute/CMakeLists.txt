###
#
# @file
#
#  CMake Information
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @author Asim YarKhan
# @date 2010-11-15
#
###

INCLUDE( ${CMAKE_SOURCE_DIR}/GenerateZDCS.cmake OPTIONAL )

set( plasma_compute_ZSRC_BASE
        gels.c  gelqf.c gelqs.c geqrf.c geqrs.c
        gesv_incpiv.c getrf_incpiv.c getrs_incpiv.c
        gesv.c  getrf.c getrs.c laswp.c laswpc.c
        posv.c  potrf.c potrs.c
        trsm.c  gemm.c  symm.c syrk.c trsmpl.c trsmrv.c
        potri.c trtri.c lauum.c trmm.c syr2k.c
        getri.c
        lange.c lansy.c laset.c lacpy.c
        gecfi2.h
        gecfi2.c gecfi.c getmi.c
        plrnt.c plgsy.c
        gesvd.c gebrd.c
)
set( plasma_compute_PZSRC_BASE
        gelqf.c geqrf.c potrf.c trtri.c lauum.c
        gelqfrh.c geqrfrh.c
        trsm.c  gemm.c  symm.c  syrk.c  syr2k.c trsmpl.c trsmrv.c
        lange.c lansy.c
        lacpy.c geadd.c getrf_incpiv.c trmm.c
        laset.c laset2.c
        pack.c shift.c getmi2.c
        plrnt.c plgsy.c
        getrf_reclap.c getrf_rectil.c laswp.c laswpc.c
        barrier.c
        gerbb.c
        gbrdb.c
)

# Start with some manually setup sources
set( plasma_compute_ZSRCS ${plasma_compute_ZSRCS}
        zunmqr.c zunmlq.c zungqr.c zunglq.c
        pzunmqr.c pzunmlq.c pzungqr.c pzunglq.c
        pzunmqrrh.c pzungqrrh.c
        pzunmlqrh.c pzunglqrh.c
        zhemm.c  zherk.c zher2k.c zlanhe.c
        pzlanhe.c pzhemm.c pzherk.c pzher2k.c
        zplghe.c pzplghe.c
        pzlag2c.c
        zcposv.c zcgesv.c
        zcungesv.c
        zhegv.c zheev.c zhegst.c zhetrd.c
        pzhetrd_he2bhe.c pzhegst.c
        pzhbrdt.c
)
set( plasma_compute_CSRCS ${plasma_compute_CSRCS}
        cunmqr.c cunmlq.c cungqr.c cunglq.c
        pcunmqr.c pcunmlq.c pcungqr.c pcunglq.c
        pcunmqrrh.c pcungqrrh.c
        pcunmlqrh.c pcunglqrh.c
        chemm.c  cherk.c  cher2k.c clanhe.c
        pclanhe.c pchemm.c pcherk.c pcher2k.c
        cplghe.c pcplghe.c
        chegv.c cheev.c chegst.c chetrd.c
        pchetrd_he2bhe.c pchegst.c
        pchbrdt.c
        )
set( plasma_compute_DSRCS ${plasma_compute_DSRCS}
        dormqr.c dormlq.c dorgqr.c dorglq.c
        pdormqr.c pdormlq.c pdorgqr.c pdorglq.c
        pdormqrrh.c pdorgqrrh.c
        pdormlqrh.c pdorglqrh.c
        pdlag2s.c
        dsungesv.c
        dsposv.c dsgesv.c
        dsygv.c dsyev.c dsygst.c dsytrd.c
        pdsysytrd_sy2bsy.c pdsygst.c
        pdsbrdt.c
)
set( plasma_compute_SSRCS ${plasma_compute_SSRCS}
        sormqr.c sormlq.c sorgqr.c sorglq.c
        psormqr.c psormlq.c psorgqr.c psorglq.c
        psormqrrh.c psorgqrrh.c
        psormlqrh.c psorglqrh.c
        pssysytrd_sy2bsy.c pssygst.c
        ssygv.c ssyev.c ssygst.c ssytrd.c
        pssbrdt.c
)


foreach ( base ${plasma_compute_ZSRC_BASE} )
  set( plasma_compute_ZSRCS ${plasma_compute_ZSRCS} z${base} )
  set( plasma_compute_CSRCS ${plasma_compute_CSRCS} c${base} )
  set( plasma_compute_DSRCS ${plasma_compute_DSRCS} d${base} )
  set( plasma_compute_SSRCS ${plasma_compute_SSRCS} s${base} )
endforeach()

foreach ( base ${plasma_compute_PZSRC_BASE} )
  set( plasma_compute_ZSRCS ${plasma_compute_ZSRCS} pz${base} )
  set( plasma_compute_CSRCS ${plasma_compute_CSRCS} pc${base} )
  set( plasma_compute_DSRCS ${plasma_compute_DSRCS} pd${base} )
  set( plasma_compute_SSRCS ${plasma_compute_SSRCS} ps${base} )
endforeach()

set( plasma_compute_SRCS ${plasma_compute_ZSRCS} ${plasma_compute_CSRCS} ${plasma_compute_DSRCS} ${plasma_compute_SSRCS} )


##################################################

# Setup files in control directory

set( plasma_control_HEADERS
        allocate.h
        async.h
        auxiliary.h
        common.h
        context.h
        control.h
        descriptor.h
        global.h
        tile.h
        workspace.h primes.h gkkleader.h
#####        workspace_z.h workspace_d.h workspace_c.h workspace_s.h
        compute_z.h   compute_d.h   compute_c.h   compute_s.h
)

set( plasma_control_ALLSRC
        allocate.c async.c    auxiliary.c
        context.c  control.c  descriptor.c
        tile.c     plasmaos.c workspace.c
        primes.c   gkkleader.c
        plasma_f77.c plasma_mf77.c
        plasma_zf77.c plasma_df77.c plasma_cf77.c plasma_sf77.c
        workspace_z.c workspace_d.c workspace_c.c workspace_s.c
        pztile.c      pdtile.c      pctile.c      pstile.c
        ztile.c       dtile.c       ctile.c       stile.c
)

foreach( base ${plasma_control_HEADERS} )
  set( plasma_control_SRCS ${plasma_control_SRCS} ${CMAKE_SOURCE_DIR}/control/${base} )
endforeach()
foreach( base ${plasma_control_ALLSRC} )
  set( plasma_control_SRCS ${plasma_control_SRCS} ${CMAKE_SOURCE_DIR}/control/${base} )
endforeach()
if(WIN32)
  set( plasma_control_SRCS ${plasma_control_SRCS} ${CMAKE_SOURCE_DIR}/control/plasmawinthread.c ${CMAKE_SOURCE_DIR}/control/plasmawinthread.h )
endif( WIN32 )

set( plasma_SRCS ${plasma_control_SRCS} ${plasma_compute_SRCS} )

add_library( plasma ${plasma_SRCS} )

add_dependencies( plasma include control_include )

INSTALL( TARGETS plasma DESTINATION ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin )
INSTALL( FILES ${plasma_compute_SRCS} CMakeLists.txt DESTINATION compute )
INSTALL( FILES ${plasma_control_SRCS} DESTINATION control )
