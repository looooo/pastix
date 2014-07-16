include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../sopalin")

# api_str_to_int
add_custom_command(
  OUTPUT  api_str_to_int.c
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../tools/genapi.sh ${CMAKE_CURRENT_SOURCE_DIR}/../include/pastix/api.h api_str_to_int.c
  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/../tools/genapi.sh
  ${CMAKE_CURRENT_SOURCE_DIR}/../include/pastix/api.h
)

set(MTXDRIVER_LIB_SRCS
  skitf.f
  iohb.c
  mmio.c
  common_drivers.c
  get_options.c
  read_matrix.c
  rsaread.c
  hbread.c
  mmread.c
  mmdread.c
  petscread.c
  cccread.c
  olafread.c
  chbread.c
  cscdread.c
  peerread.c
  threefilesread.c
  fdupread.c
  laplacian.c
  api_str_to_int.c
)

set_source_files_properties(mmio.c PROPERTIES COMPILE_FLAGS "-std=c99" )
set_source_files_properties(read_matrix.c PROPERTIES COMPILE_FLAGS "-std=c99" )

add_library(read_matrix
  ${MTXDRIVER_LIB_SRCS}
)

add_dependencies(read_matrix
  headers)

install(TARGETS read_matrix
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

