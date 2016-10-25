include(RulesPrecisions)

### reset variables
set(generated_sources "")
set(generated_headers "")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/drivers")

set(HEADERS
  z_spm.h
)

### Generate the headers in all precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "p;s;d;c;z")

set(SOURCES
  z_spm.c
  z_spm_2dense.c
  z_spm_convert_to_csc.c
  z_spm_convert_to_csr.c
  z_spm_convert_to_ijv.c
  z_spm_dofs2flat.c
  z_spm_expand.c
  z_spm_genrhs.c
  z_spm_integer.c
  z_spm_laplacian.c
  z_spm_matrixvector.c
  z_spm_norm.c
  z_spm_print.c
  )

set(spm_headers
  ${generated_headers}
  spm.h
  spm_drivers.h
  )

### Generate the sources in all precisions
precisions_rules_py(generated_sources
                    "${SOURCES}"
                    PRECISIONS "p;s;d;c;z")

set(spm_sources
  ${generated_sources}
  spm.c
  spm_io.c
  spm_integers.c
  spm_expand.c
  spm_read_driver.c
  drivers/skitf.f
  drivers/iohb.c
  drivers/mmio.c
  drivers/laplacian.c
  drivers/readhb.c
  drivers/readijv.c
  drivers/readmm.c
  drivers/readrsa.c
  )

add_custom_target(spm_headers_tgt
  DEPENDS ${spm_headers} )

add_library(pastix_spm
  ${spm_sources}
  )

add_dependencies(pastix_spm
  spm_headers_tgt
)

### Generate the lib
if (MPI_FOUND)
  set_target_properties(pastix_spm PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif (MPI_FOUND)

install(TARGETS pastix_spm
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
