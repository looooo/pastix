include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
link_directories(
  ${CBLAS_LIBRARY_DIRS}
  ${BLAS_LIBRARY_DIRS}
)

set(PASTIX_EXAMPLES_LIBRARIES pastix)

if(PASTIX_WITH_STARPU)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    pastix_starpu)
endif(PASTIX_WITH_STARPU)

set(PASTIX_EXAMPLES_LIBRARIES
  ${PASTIX_EXAMPLES_LIBRARIES}
  pastix_kernels
  pastix_libcsc
  pastix_blend
  pastix # pastix_blend requires dofInit from symbol which is in libpastix.a
  read_matrix
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
  endif()

  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    ${SCOTCH_LIBRARIES}
    )
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

set(HEADERS zutils.h)
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")
message(STATUS ${SOURCES})
message(STATUS ${generated_headers})
add_custom_target(examples_headers ALL SOURCES
  ${generated_headers} )
set(SOURCES zsimple.c zstep-by-step.c)
## Require all PaStiX with multiple precision before...
##precisions_rules_py(generated_files
##   "${SOURCES}"
##   PRECISIONS "s;d;c;z")
message(STATUS ${generated_files})
foreach(example_with_ext ${generated_files} simple.c step-by-step.c analyse.c)
  message(STATUS "${example_with_ext}")
  GET_FILENAME_COMPONENT(example ${example_with_ext} NAME_WE)
  add_executable(${example} "${example_with_ext}")
  add_dependencies(${example} ${example_with_ext})
  target_link_libraries( ${example}
    ${PASTIX_EXAMPLES_LIBRARIES} )

  if (PASTIX_WITH_MPI)
    target_link_libraries( ${example}
      ${MPI_LIBRARIES} )
  endif()
endforeach()
