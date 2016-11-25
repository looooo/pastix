include(RulesPrecisions)

## reset variables
set(generated_files "")

## Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_spm_tests.c
  z_bcsc_tests.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "p;s;d;c;z")

add_library(spm_test ${generated_files} )

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
  target_link_libraries(${_name_we} spm_test pastix)
endforeach()


## CTest execution
set( SPM_TESTS
  spm_convert_tests spm_norm_tests spm_matvec_tests )
set( SPM_DOF_TESTS
  spm_dof_expand_tests spm_dof_norm_tests spm_dof_matvec_tests)
set( BCSC_TESTS
  bcsc_norm_tests bcsc_matvec_tests )

# Laplacian
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} ${BCSC_TESTS} )
  add_test(lap_${example} ./${example} -lap 10:10:10)
endforeach()

# RSA
foreach(example ${SPM_TESTS} ${BCSC_TESTS} )
  add_test(rsa_${example} ./${example} -rsa ${CMAKE_SOURCE_DIR}/test/matrix/small.rsa)
endforeach()

# Matrix Market
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} ${BCSC_TESTS} )
  add_test(mm_${example} ./${example} -mm ${CMAKE_SOURCE_DIR}/test/matrix/young4c.mtx)
endforeach()

# Harwell Boeing
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} ${BCSC_TESTS} )
  add_test(hb_${example} ./${example} -hb ${CMAKE_SOURCE_DIR}/test/matrix/orsirr.rua)
endforeach()

# Matrix Market - Hermitian
foreach(example ${SPM_TESTS} ${SPM_DOF_TESTS} ${BCSC_TESTS} )
  add_test(mm2_${example} ./${example} -mm ${CMAKE_SOURCE_DIR}/test/matrix/mhd1280b.mtx)
endforeach()
