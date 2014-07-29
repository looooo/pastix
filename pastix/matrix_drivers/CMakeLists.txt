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

set(HEADERS
  z_cccread.h
  z_chbread.h
  z_cscdread.h
  #z_csparse.h
  z_fdupread.h
  z_hbread.h
  z_laplacian.h
  z_mmdread.h
  z_mmread.h
  #z_mtx.h
  z_olafread.h
  z_peerread.h
  z_petscread.h
  z_read_matrix.h
  z_rsaread.h
  z_threefilesread.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(matrix_drivers_headers
  DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_cccread.c
  z_chbread.c
  z_cscdread.c
  #z_csparse.c
  z_fdupread.c
  z_hbread.c
  z_laplacian.c
  z_mmdread.c
  z_mmread.c
  #z_mtx.c
  z_olafread.c
  z_peerread.c
  z_petscread.c
  z_read_matrix.c
  z_rsaread.c
  z_threefilesread.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")
set(MTXDRIVER_LIB_SRCS
  skitf.f
  iohb.c
  mmio.c
  common_drivers.c
  get_options.c
  api_str_to_int.c
  ${generated_files}
)

set_source_files_properties(mmio.c PROPERTIES COMPILE_FLAGS "-std=c99" )
set_source_files_properties(read_matrix.c PROPERTIES COMPILE_FLAGS "-std=c99" )

add_library(read_matrix
  ${MTXDRIVER_LIB_SRCS}
)

add_dependencies(read_matrix matrix_drivers_headers)
add_dependencies(read_matrix include_headers)
add_dependencies(read_matrix common_headers)

install(TARGETS read_matrix
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

