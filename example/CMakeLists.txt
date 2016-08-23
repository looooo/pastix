include(RulesPrecisions)
include(array2d)

# reset variables
set(generated_files "")
set(generated_headers "")

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_spm_matvec_test.c
  z_bcsc_tests.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "p;s;d;c;z")

add_library(spm_test ${generated_files} )

include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}/")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/")
link_directories(
  ${LAPACKE_LIBRARY_DIRS}
  ${LAPACK_LIBRARY_DIRS}
  ${CBLAS_LIBRARY_DIRS}
  ${BLAS_LIBRARY_DIRS}
)

message(STATUS "LAPACKE_LIBRARIES:    ${LAPACKE_LIBRARIES}")
message(STATUS "LAPACK_SEQ_LIBRARIES: ${LAPACK_SEQ_LIBRARIES}")
message(STATUS "CBLAS_LIBRARIES:      ${CBLAS_LIBRARIES}")
message(STATUS "BLAS_SEQ_LIBRARIES:   ${BLAS_SEQ_LIBRARIES}")

set(PASTIX_EXAMPLES_LIBRARIES
  pastix
  ${LAPACKE_LIBRARIES}
  ${LAPACK_SEQ_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  ${EXTRA_LIBS}
  m
)

if(PASTIX_WITH_STARPU)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    pastix_starpu)
endif(PASTIX_WITH_STARPU)

message(STATUS "Parsec libs: ${PARSEC_LIBRARIES}")
if(PASTIX_WITH_PARSEC)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    pastix_parsec
    ${PARSEC_LIBRARIES}
    dague_distribution
    dague_distribution_matrix
    ${PARSEC_EXTRA_LIBRARIES}
    )
endif(PASTIX_WITH_PARSEC)

if (PASTIX_ORDERING_SCOTCH)
  link_directories(${SCOTCH_LIBRARY_DIRS})

  if (PASTIX_ORDERING_PTSCOTCH)
    set(PASTIX_EXAMPLES_LIBRARIES
      ${PASTIX_EXAMPLES_LIBRARIES}
      ${PTSCOTCH_LIBRARIES}
      )
    # TODO IF 6.0.0 APPEND SCOTCH_LIBRARIES
  else()
    set(PASTIX_EXAMPLES_LIBRARIES
      ${PASTIX_EXAMPLES_LIBRARIES}
      ${SCOTCH_LIBRARIES}
      )
  endif()
endif()

if (PASTIX_ORDERING_METIS)
  link_directories(${METIS_DIR}/lib)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    ${METIS_LIBRARIES}
    )
endif()

if (PASTIX_BLEND_GENTRACE)
  link_directories(${GTG_DIR}/lib)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    ${GTG_LIBRARIES}
    )
endif()

if(PASTIX_WITH_STARPU)
  link_directories(${STARPU_LIBRARY_DIRS})
  if (PASTIX_WITH_MPI)
    set(PASTIX_EXAMPLES_LIBRARIES
      ${PASTIX_EXAMPLES_LIBRARIES}
      ${STARPU_MPI_LIBRARIES}
      )
  else(PASTIX_WITH_MPI)
    set(PASTIX_EXAMPLES_LIBRARIES
      ${PASTIX_EXAMPLES_LIBRARIES}
      ${STARPU_SHM_LIBRARIES}
      )
  endif(PASTIX_WITH_MPI)
endif(PASTIX_WITH_STARPU)

if (PASTIX_WITH_MPI)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    ${MPI_LIBRARIES}
    )
endif()

set(HEADERS z_utils.h)
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")
add_custom_target(examples_headers DEPENDS ${generated_headers} )

# set(SOURCES
#   z_simple.c
#   z_step-by-step.c
#   z_analyze.c
# )

# ## Require all PaStiX with multiple precision before...
# precisions_rules_py(generated_files
#    "${SOURCES}"
#    PRECISIONS "s;d;c;z")

# foreach(example_with_ext ${generated_files})
#   GET_FILENAME_COMPONENT(example ${example_with_ext} NAME_WE)
#   add_executable(${example} "${example_with_ext}" ${generated_headers})
#   target_link_libraries( ${example}
#     ${PASTIX_EXAMPLES_LIBRARIES})

#   if (PASTIX_WITH_MPI)
#     target_link_libraries( ${example}
#       ${MPI_LIBRARIES} )
#   endif()
# endforeach()

add_executable(analyze analyze.c)
target_link_libraries(analyze ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(simple simple.c)
target_link_libraries(simple ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(compare_tolerance compare_tolerance.c)
target_link_libraries(compare_tolerance ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(spm_convert_tests spm_convert_tests.c)
target_link_libraries(spm_convert_tests ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(spm_matvec_tests spm_matvec_tests.c)
target_link_libraries(spm_matvec_tests spm_test ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(spm_norm_tests spm_norm_tests.c)
target_link_libraries(spm_norm_tests spm_test ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(bcsc_norm_tests bcsc_norm_tests.c)
target_link_libraries(bcsc_norm_tests spm_test ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(bcsc_matvec_tests bcsc_matvec_tests.c)
target_link_libraries(bcsc_matvec_tests spm_test ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(compress_bfacto_test compress_bfacto_test.c)
target_link_libraries(compress_bfacto_test ${PASTIX_EXAMPLES_LIBRARIES})

add_executable(compress_afacto_test compress_afacto_test.c)
target_link_libraries(compress_afacto_test ${PASTIX_EXAMPLES_LIBRARIES})

#
#-- CTests
#
enable_testing()
include(CTest)

foreach(example simple step-by-step simple_dist simple_fortran)
  foreach(prec s d c z)
    set(test_args
      ${test_args}
      ${prec}_${example} "-lap 1000"
      )
  endforeach()
endforeach()

foreach(example murge_2Dlap murge_3Dlap)
  foreach(prec s d c z)
    set(test_args
      ${test_args}
      ${prec}${example} "-n 10")
  endforeach()
endforeach()

foreach(prec s d c z)
  set(test_args
    ${test_args}
    ${prec}_murge_fortran "100 1")
  set(test_args
    ${test_args}
    ${prec}_murge_fortran "100 4")
endforeach()

array2d_begin_loop( advanced "${test_args}" 2 "test;args" )
while( advanced )
  set(name "${test}_${args}")
  string(REPLACE " " "_" name ${name})
  string(REPLACE " " ";" args_list ${args})
  add_test(${name} ${test} ${args_list})
  set_tests_properties(${name} PROPERTIES TIMEOUT 4)
  if (PASTIX_WITH_MPI)
    add_test("mpirun_np_2_${name}" ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${test} ${MPIEXEC_POSTFLAGS} ${args_list})
    set_tests_properties(${name} PROPERTIES TIMEOUT 4)
  endif(PASTIX_WITH_MPI)
  array2d_advance()
endwhile()
