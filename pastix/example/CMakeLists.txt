include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
link_directories(
  ${CBLAS_LIBRARY_DIRS}
  ${BLAS_LIBRARY_DIRS}
)

set(PASTIX_EXAMPLES_LIBRARIES
  pastix
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

add_executable(simple
  simple.c
)

add_executable(analyse
  analyse.c
)

target_link_libraries( simple
  ${PASTIX_EXAMPLES_LIBRARIES} )

if (PASTIX_WITH_MPI)
  target_link_libraries( simple
    ${MPI_LIBRARIES} )
endif()

target_link_libraries( analyse
  ${PASTIX_EXAMPLES_LIBRARIES} )

if (PASTIX_WITH_MPI)
  target_link_libraries( analyse
    ${MPI_LIBRARIES} )
endif()
