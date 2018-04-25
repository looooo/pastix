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

## Generate the test library files for all required precisions
set(LIB_SOURCES
  z_bcsc_tests.c
)

## reset variables
set(generated_libfiles "")

precisions_rules_py(generated_libfiles
   "${LIB_SOURCES}"
   PRECISIONS "p;s;d;c;z")

add_library(bcsc_test
   ${generated_libfiles}
   )

target_link_libraries( bcsc_test
  pastix_kernels
  pastix_bcsc
  ${SPM_LIBRARY}
  ${LAPACKE_LIBRARIES_DEP}
  )

## Generate all test executables
set (TESTS
  bcsc_norm_tests.c
  bcsc_matvec_tests.c
  )

foreach (_file ${TESTS})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix bcsc_test)
endforeach()

## CTest execution
set( BCSC_TESTS
  bcsc_norm_tests bcsc_matvec_tests )

# Laplacian
foreach(example ${BCSC_TESTS} )
  foreach(arithm ${PASTIX_PRECISIONS} )
    add_test(test_lap_${arithm}_${example} ./${example} --lap ${arithm}:10:10:10)
  endforeach()
endforeach()

# RSA
if (PASTIX_WITH_FORTRAN)
  foreach(example ${BCSC_TESTS} )
    add_test(test_rsa_${example} ./${example} --rsa ${CMAKE_SOURCE_DIR}/test/matrix/small.rsa)
  endforeach()
endif()

# Matrix Market
foreach(example ${BCSC_TESTS} )
  add_test(test_mm_${example} ./${example} --mm ${CMAKE_SOURCE_DIR}/test/matrix/young4c.mtx)
endforeach()

# Harwell Boeing
foreach(example ${BCSC_TESTS} )
  add_test(test_hb_${example} ./${example} --hb ${CMAKE_SOURCE_DIR}/test/matrix/orsirr.rua)
endforeach()

# Matrix Market - Hermitian
foreach(example ${BCSC_TESTS} )
  add_test(test_mm2_${example} ./${example} --mm ${CMAKE_SOURCE_DIR}/test/matrix/mhd1280b.mtx)
endforeach()

###################################
#   Low-rank testings
###

## Tests requires TMG lib to generate random matrices.
find_package(TMG)

if(TMG_FOUND)
  ## Generate the test library files for all required precisions
  set(LIB_SOURCES
    z_lowrank_tests.c
    )

  ## reset variables
  set(generated_libfiles "")

  precisions_rules_py(generated_libfiles
    "${LIB_SOURCES}"
    PRECISIONS "p;s;d;c;z")

  add_library(lowrank_test
    ${generated_libfiles}
    )

  target_link_libraries( lowrank_test
    pastix_kernels
    ${LAPACKE_LIBRARIES_DEP}
    ${TMG_LIBRARIES_DEP}
    ${EXTRA_LIBS}
    m)

  set(TESTS_SOURCES
    z_ge2lr_tests.c
    z_rradd_tests.c
    z_lrmm_tests.c
    )

  ## reset variables
  set(generated_exfiles "")

  precisions_rules_py(generated_exfiles
    "${TESTS_SOURCES}"
    PRECISIONS "p;s;d;c;z")

  foreach (_file ${generated_exfiles})
    get_filename_component(_name_we ${_file} NAME_WE)
    add_executable(${_name_we} ${_file})
    target_link_libraries(${_name_we}
      pastix_kernels
      pastix
      lowrank_test
      ${LAPACKE_LIBRARIES_DEP}
      ${TMG_LIBRARIES_DEP}
      ${EXTRA_LIBS}
      m)

    add_test(test_${_name_we} ./${_name_we})
  endforeach()

endif()

