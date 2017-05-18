###
#
#  @copyright 2013-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @date 2013-06-24
#
###
add_subdirectory(old)
set (EXAMPLES
  analyze.c
  simple.c
  reentrant.c
  bench_facto.c
  step-by-step.c
  schur.c
  personal.c
  )

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix)

  install(TARGETS ${_name_we} RUNTIME DESTINATION examples )
  install(FILES   ${_file}            DESTINATION examples )
endforeach()

### CTest execution
set( PASTIX_DRIVERS
  analyze simple step-by-step schur personal) # reentrant fails with Scotch parser !
set( PASTIX_TESTS
  simple )
if (PASTIX_WITH_PARSEC)
  set( PASTIX_SCHEDS
    0 1 2) # sequential, static, parsec
else()
  set( PASTIX_SCHEDS
    0 1)   # sequential, static
endif (PASTIX_WITH_PARSEC)

### All drivers with Laplacian and default parameters
foreach(example ${PASTIX_DRIVERS} )
  add_test(drv_${example} ./${example} -lap 10:10:10)
endforeach()

### Factotype
foreach(example ${PASTIX_TESTS} )
  # RSA
  add_test(rsa_${example} ./${example} -rsa ${CMAKE_SOURCE_DIR}/test/matrix/small.rsa   -f 2) # LLt fails with negative pivot !
  # Matrix Market
  add_test(mm_${example}  ./${example} -mm ${CMAKE_SOURCE_DIR}/test/matrix/young4c.mtx  -f 2) # LLt fails with negative pivot !
  # Harwell Boeing
  add_test(hb_${example}  ./${example} -hb ${CMAKE_SOURCE_DIR}/test/matrix/orsirr.rua   -f 2)
  # Matrix Market - Hermitian
  add_test(mm2_${example} ./${example} -mm ${CMAKE_SOURCE_DIR}/test/matrix/mhd1280b.mtx -f 2)
endforeach()

### Refinement
foreach(example ${PASTIX_TESTS} )
  add_test(cg_${example}       ./${example} -rsa ${CMAKE_SOURCE_DIR}/test/matrix/small.rsa -f 2 -i iparm_refinement pastixrefinecg) # LLt fails with negative pivot !
  add_test(gmres_${example}    ./${example} -hb ${CMAKE_SOURCE_DIR}/test/matrix/orsirr.rua -f 2 -i iparm_refimenent pastixrefinegmres)
  add_test(bicgstab_${example} ./${example} -hb ${CMAKE_SOURCE_DIR}/test/matrix/orsirr.rua -f 2 -i iparm_refinement pastixrefinebicgstab)
endforeach()

### Distribution
if (PASTIX_WITH_PARSEC)
  foreach(example ${PASTIX_TESTS} )
    # 2D
    add_test(2d_${example}   ./${example} -lap 30:30:30 -s 2 -i iparm_distribution_level 0 )
    # 1D/2D
    add_test(1d2d_${example} ./${example} -lap 30:30:30 -s 2 -i iparm_distribution_level 16)
  endforeach()
endif (PASTIX_WITH_PARSEC)

### Scheduling
foreach(example ${PASTIX_TESTS} )
  foreach(scheduler ${PASTIX_SCHEDS} )
    add_test(sched${scheduler}_${example}   ./${example} -lap 30:30:30 -s ${scheduler})
  endforeach()
endforeach()

# Low Rank
foreach(example ${PASTIX_TESTS} )
  foreach(scheduler ${PASTIX_SCHEDS} )
    # SVD Begin
    add_test(svdbegin_sched${scheduler}_${example}  ./${example} -lap 30:30:30 -s ${scheduler} -i iparm_compress_when pastixcompresswhenbegin -i iparm_compress_method pastixcompressmethodsvd)
    # SVD End
    add_test(svdend_sched${scheduler}_${example}    ./${example} -lap 30:30:30 -s ${scheduler} -i iparm_compress_when pastixcompresswhenend   -i iparm_compress_method pastixcompressmethodsvd)
    # RRQR Begin
    add_test(rrqrbegin_sched${scheduler}_${example} ./${example} -lap 30:30:30 -s ${scheduler} -i iparm_compress_when pastixcompresswhenbegin -i iparm_compress_method pastixcompressmethodrrqr)
    # RRQR End
    add_test(rrqrend_sched${scheduler}_${example}   ./${example} -lap 30:30:30 -s ${scheduler} -i iparm_compress_when pastixcompresswhenend   -i iparm_compress_method pastixcompressmethodrrqr)
    endforeach()
endforeach()
