include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  z_csc_intern_io.h
  z_csc.h
  z_csc_intern_compute.h
  z_csc_intern_updown.h
  z_csc_intern_build.h
  z_csc_intern_solve.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(bcsc_headers 
  DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_csc_intern_compute.c
  z_csc_intern_compute_ge.c
  z_csc_intern_compute_po.c
  z_csc_intern_compute_sy.c
  z_csc_intern_compute_he.c
  z_csc_intern_build.c
  z_csc_intern_io.c
  z_csc_intern_solve.c
  z_csc_intern_updown.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

##list(REMOVE_ITEM generated_files
##  s_csc_intern_compute.c
##  d_csc_intern_compute.c
##  c_csc_intern_compute.c
##  z_csc_intern_compute.c
##)

### Generate the lib
add_library(pastix_bcsc ${generated_files} )

add_dependencies(pastix_bcsc
  bcsc_headers
)

install(TARGETS pastix_bcsc
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


