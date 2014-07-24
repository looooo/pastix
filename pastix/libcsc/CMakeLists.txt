include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  z_csc_utils.h
  z_cscd_utils.h
  z_cscd_utils_intern.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(libcsc_headers DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_csc_utils.c
  z_cscd_utils.c
  z_cscd_utils_fortran.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

### Generate the lib
add_library(pastix_libcsc ${generated_files} )

add_dependencies(pastix_libcsc
  libcsc_headers
)

install(TARGETS pastix_libcsc
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


