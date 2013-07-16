include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
link_directories(
  /opt/plasma-extlibs-dbg/lib
)

set(PASTIX_EXAMPLES_LIBRARIES
  pastix
  read_matrix
  refblas
  hwloc
  pthread
  m
)

if (PASTIX_ORDERING_SCOTCH)
  link_directories(/opt/scotch/int64/lib)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    scotch
    scotcherr
    )
endif()

if (PASTIX_ORDERING_METIS)
  link_directories(/opt/metis/lib)
  set(PASTIX_EXAMPLES_LIBRARIES
    ${PASTIX_EXAMPLES_LIBRARIES}
    metis
    )
endif()

add_executable(simple
  simple.c
)

target_link_libraries( simple
  ${PASTIX_EXAMPLES_LIBRARIES} )

if (PASTIX_WITH_MPI)
  target_link_libraries( simple
    ${MPI_LIBRARIES} )
endif()