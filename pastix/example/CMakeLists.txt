include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
link_directories(
  ${BLAS_LIBRARY_DIRS}
)

set(PASTIX_EXAMPLES_LIBRARIES
  pastix
  read_matrix
  hwloc
  ${BLAS_LIBRARIES}
  pthread
  rt
  m
)

if (PASTIX_ORDERING_SCOTCH)
  link_directories(${SCOTCH_DIR}/lib)

  if (PASTIX_ORDERING_PTSCOTCH)
    set(PASTIX_EXAMPLES_LIBRARIES
      ${PASTIX_EXAMPLES_LIBRARIES}
      ptscotch
      ptscotcherr
      )
  endif()

  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    scotch
    scotcherr
    )
endif()


if (PASTIX_ORDERING_METIS)
  link_directories(${METIS_DIR}/lib)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    metis
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
