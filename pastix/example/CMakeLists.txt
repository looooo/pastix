include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
link_directories(
  /opt/plasma-extlibs-dbg/lib
  /opt/scotch/int64/lib
)

add_executable(simple
  simple.c
)


target_link_libraries( simple
    pastix read_matrix refblas scotch scotcherr hwloc pthread m )
