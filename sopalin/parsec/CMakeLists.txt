include(RulesPrecisions)
include(RulesJDF)
# reset variables
set(generated_files "")
set(generated_jdf "")
set(generated_headers "")

set(EXTRA_SOURCES
  sparse-matrix.c
)

set(HEADERS
  sparse-matrix.h
)

# precisions_rules_py(generated_headers
#                  "${HEADERS}"
#                  PRECISIONS "z;c;d;s")

### Generate .c files from .jdf for all required precisions
set(JDF
  # zcsc2cblk.jdf
  zpotrf_sp1dplus.jdf
  zgetrf_sp1dplus.jdf
  zhetrf_sp1dplus.jdf
  zsytrf_sp1dplus.jdf
  #
  zpotrf_sp2d.jdf
  zgetrf_sp2d.jdf
  # # Solvers
  # zpotrs_sp1dplus.jdf
)
precisions_rules_py(generated_jdf
  "${JDF}"
  PRECISIONS "s;d;c;z")

jdf_rules(generated_files "${generated_jdf}")

### Generate the pastix_parsec wrappers for all required precisions
set(SOURCES
  #zcsc2cblk_wrapper.c
  zsytrf_sp1dplus_wrapper.c
  zhetrf_sp1dplus_wrapper.c
  zpotrf_sp1dplus_wrapper.c
  zgetrf_sp1dplus_wrapper.c
  #zsytrf_sp2d_wrapper.c
  #zhetrf_sp2d_wrapper.c
  #zpotrf_sp2d_wrapper.c
  zgetrf_sp2d_wrapper.c
  #zpotrs_sp1dplus_wrapper.c
)

precisions_rules_py(generated_files
  "${SOURCES}"
  PRECISIONS "s;d;c;z")

add_library(pastix_parsec
  ${generated_headers}
  ${generated_files}
  ${EXTRA_SOURCES})

target_link_libraries(pastix_parsec pastix_kernels)

add_dependencies(pastix_parsec
  sopalin_headers
)

### Generate the lib
if (MPI_FOUND)
  set_target_properties(pastix_parsec PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif (MPI_FOUND)

install(TARGETS pastix_parsec
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
