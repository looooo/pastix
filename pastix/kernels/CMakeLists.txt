include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  pastix_zcores.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  # extra kernel
  core_zgemdm.c
  core_zgetro.c
  core_zgeadd.c
  # Kernels
  core_zsytrfsp.c
  core_zhetrfsp.c
  core_zpotrfsp.c
  core_zgetrfsp.c
)
# precisions_rules_py(generated_files
#   "${SOURCES}"
#   PRECISIONS "s;d;c;z")

### Generate the lib
add_library(dsparse_cores
  ${SOURCES}
#  ${generated_files}
)
install(TARGETS dsparse_cores ARCHIVE DESTINATION lib)

