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
include(RulesPrecisions)

## reset variables
set(generated_files "")
set(generated_files2 "")

## Generate the dsparse wrappers for all required precisions
set(LIB_SOURCES
  z_spm_tests.c
  z_bcsc_tests.c
)

set(TESTS_SOURCES
  z_ge2lr_tests.c
  z_rradd_tests.c
)

precisions_rules_py(generated_libfiles
   "${LIB_SOURCES}"
   PRECISIONS "p;s;d;c;z")

precisions_rules_py(generated_exfiles
   "${TESTS_SOURCES}"
   PRECISIONS "p;s;d;c;z")


add_library(spm_test ${generated_libfiles} )

target_link_libraries( spm_test
  pastix_kernels
  pastix_bcsc
  pastix_spm
  ${LAPACKE_LIBRARIES}
  ${LAPACK_SEQ_LIBRARIES} )


## Generate all test executables
set (TESTS
  spm_convert_tests.c
  spm_norm_tests.c
  spm_matvec_tests.c
  spm_dof_expand_tests.c
  spm_dof_norm_tests.c
  spm_dof_matvec_tests.c
  bcsc_norm_tests.c
  bcsc_matvec_tests.c
  )

foreach (_file ${TESTS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix spm_test)
endforeach()

foreach (_file ${generated_exfiles})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we}
    pastix_kernels
    pastix
    ${LAPACKE_LIBRARIES}
    ${TMG_LIBRARIES}
    ${LAPACK_SEQ_LIBRARIES}
    ${BLAS_SEQ_LIBRARIES}
    ${EXTRA_LIBS}
    m)
endforeach()


## CTest execution
set( SPM_TESTS
  spm_convert_tests spm_norm_tests spm_matvec_tests )
set( SPM_DOF_TESTS
  spm_dof_expand_tests spm_dof_norm_tests spm_dof_matvec_tests)
set( BCSC_TESTS
  bcsc_norm_tests bcsc_matvec_tests )

set( GE2LR_TESTS
  s_ge2lr_tests d_ge2lr_tests c_ge2lr_tests z_ge2lr_tests )
set( RRADD_TESTS
  s_rradd_tests d_rradd_tests c_rradd_tests z_rradd_tests )

# GE2LR
foreach(example ${GE2LR_TESTS} )
  add_test(test_${example} ./${example})
endforeach()

# RRADD
foreach(example ${RRADD_TESTS} )
  add_test(test_${example} ./${example})
endforeach()

# Laplacian
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} ${BCSC_TESTS} )
  add_test(test_lap_${example} ./${example} --lap 10:10:10)
endforeach()

# RSA
foreach(example ${SPM_TESTS} ${BCSC_TESTS} )
  add_test(test_rsa_${example} ./${example} --rsa ${CMAKE_SOURCE_DIR}/test/matrix/small.rsa)
endforeach()

# Matrix Market
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} ${BCSC_TESTS} )
  add_test(test_mm_${example} ./${example} --mm ${CMAKE_SOURCE_DIR}/test/matrix/young4c.mtx)
endforeach()

# Harwell Boeing
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} ${BCSC_TESTS} )
  add_test(test_hb_${example} ./${example} --hb ${CMAKE_SOURCE_DIR}/test/matrix/orsirr.rua)
endforeach()

# Matrix Market - Hermitian
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} ${BCSC_TESTS} )
  add_test(test_mm2_${example} ./${example} --mm ${CMAKE_SOURCE_DIR}/test/matrix/mhd1280b.mtx)
endforeach()
