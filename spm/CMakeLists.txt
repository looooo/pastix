include(RulesPrecisions)

### reset variables
set(generated_sources "")
set(generated_headers "")

set(HEADERS
  z_spm.h
)

### Generate the headers in all precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "p;s;d;c;z")

set(SOURCES
  z_spm.c
  z_spm_convert_to_csc.c
  z_spm_convert_to_csr.c
  z_spm_convert_to_ijv.c
  z_spm_genrhs.c
  z_spm_matrixvector.c
  z_spm_norm.c
  )

set(spm_headers
  ${generated_headers}
  spm.h
  )

### Generate the sources in all precisions
precisions_rules_py(generated_sources
                    "${SOURCES}"
                    PRECISIONS "p;s;d;c;z")

set(spm_sources
  ${generated_sources}
  spm.c
  spm_io.c
  )

add_custom_target(spm_headers_tgt
  DEPENDS ${spm_headers} )

set_source_files_properties(
  spm.c
  PROPERTIES DEPENDS spm_headers_tgt
  )

set(PASTIX_LIB_SRCS
  ${PASTIX_LIB_SRCS}
  ${spm_sources}
  )
