include(RulesPrecisions)
include(array2d)

# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}/")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/")
link_directories(
  ${CBLAS_LIBRARY_DIRS}
  ${BLAS_LIBRARY_DIRS}
)

set(PASTIX_EXAMPLES_LIBRARIES pastix)


set(PASTIX_EXAMPLES_LIBRARIES
  ${PASTIX_EXAMPLES_LIBRARIES}
  read_matrix
  pastix_murge
  pastix
  pastix_sopalin
  pastix_libcsc
  pastix_csc
  pastix
  pastix_blend)
if(PASTIX_WITH_STARPU)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    pastix_starpu)
endif(PASTIX_WITH_STARPU)
set(PASTIX_EXAMPLES_LIBRARIES
  ${PASTIX_EXAMPLES_LIBRARIES}
  pastix_libcsc
  pastix
  pastix_kernels
  pastix_common
  ${CBLAS_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  ${EXTRA_LIBS}
  m
)

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

set(HEADERS z_utils.h)
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")
add_custom_target(examples_headers DEPENDS ${generated_headers} )

#Generate examples utils library
set(SOURCES  z_utils.F90)
precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")
add_library(pastix_examples ${generated_files} mem_trace.c)
add_dependencies(pastix_examples generated_headers)
set(PASTIX_EXAMPLES_LIBRARIES ${PASTIX_EXAMPLES_LIBRARIES} pastix_examples read_matrix)
set(generated_files "")
message(STATUS PASTIX_EXAMPLES_LIBRARIES :  ${PASTIX_EXAMPLES_LIBRARIES})



set(SOURCES z_simple.c z_step-by-step.c z_analyze.c zmurge_2Dlap.c zmurge_3Dlap.c
  z_murge_fortran.F90 z_simple_fortran.F90)
if(PASTIX_DISTRIBUTED)
  set(SOURCES ${SOURCES} z_simple_dist.c)
endif(PASTIX_DISTRIBUTED)
precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

foreach(example_with_ext ${generated_files})
  GET_FILENAME_COMPONENT(example ${example_with_ext} NAME_WE)
  GET_FILENAME_COMPONENT(ext ${example_with_ext} EXT)
  if (PASTIX_WITH_CUDA)
    cuda_add_executable(${example} "${example_with_ext}" ${generated_headers})
    target_link_libraries(${example} cuda)
  else(PASTIX_WITH_CUDA)
    add_executable(${example} "${example_with_ext}" ${generated_headers})
  endif(PASTIX_WITH_CUDA)
  target_link_libraries( ${example}
    ${PASTIX_EXAMPLES_LIBRARIES})

  if (PASTIX_WITH_MPI)
    if (${ext} STREQUAL ".c")
      message(STATUS "${example}: ${MPI_C_LIBRARIES}")
      target_link_libraries( ${example}
        ${MPI_C_LIBRARIES} )
    elseif(${ext} STREQUAL ".F90")
      message(STATUS "${example}: ${MPI_Fortran_LIBRARIES}")
      set_property(TARGET ${example} PROPERTY LINKER_LANGUAGE Fortran) 
      target_link_libraries( ${example}
        ${MPI_Fortran_LIBRARIES} )
    else()
      message(STATUS "${example} : _${ext}_")
    endif()
  endif()
endforeach()

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
