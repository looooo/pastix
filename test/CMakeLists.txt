include(RulesPrecisions)

## reset variables
set(generated_files "")

## Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_spm_matvec_test.c
  z_bcsc_tests.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "p;s;d;c;z")

add_library(spm_test ${generated_files} )

## Generte all test executable
add_executable(spm_convert_tests spm_convert_tests.c)
target_link_libraries(spm_convert_tests pastix)

set (TESTS
  spm_matvec_tests.c
  spm_norm_tests.c
  bcsc_norm_tests.c
  bcsc_matvec_tests.c
  )

foreach (_file ${TESTS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} spm_test pastix)
endforeach()


## CTest execution
# Laplacian
foreach(example spm_convert_tests spm_matvec_tests spm_norm_tests bcsc_norm_tests bcsc_matvec_tests)
  add_test(lap_${example} ./${example} -lap 1000)
endforeach()

# RSA
foreach(example spm_convert_tests spm_matvec_tests spm_norm_tests bcsc_norm_tests bcsc_matvec_tests)
  add_test(rsa_${example} ./${example} -rsa ${CMAKE_SOURCE_DIR}/test/matrix/small.rsa)
endforeach()

# Matrix Market
foreach(example spm_convert_tests spm_matvec_tests spm_norm_tests bcsc_norm_tests bcsc_matvec_tests)
  add_test(mm_${example} ./${example} -mm ${CMAKE_SOURCE_DIR}/test/matrix/young4c.mtx)
endforeach()

# Harwell Boeing
foreach(example spm_convert_tests spm_matvec_tests spm_norm_tests bcsc_norm_tests bcsc_matvec_tests)
  add_test(hb_${example} ./${example} -hb ${CMAKE_SOURCE_DIR}/test/matrix/orsirr.rua)
endforeach()
