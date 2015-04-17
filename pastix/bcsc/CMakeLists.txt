include(RulesPrecisions)

# reset variables
set(generated_files   "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

# set(HEADERS
#   z_csc_intern_io.h
#   z_csc.h
#   z_csc_intern_compute.h
#   z_csc_intern_updown.h
#   z_csc_intern_build.h
#   z_csc_intern_solve.h
#   # This C file is included in all other _he,_sy,_po,_ge
#   z_csc_intern_compute.c
# )

# ### generate the dsparse_cores headers for all possible precisions
# precisions_rules_py(generated_headers
#                     "${HEADERS}"
#                     PRECISIONS "s;d;c;z")

# add_custom_target(bcsc_headers
#   DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  bcsc_zinit.c # Not finished
)

precisions_rules_py(generated_files
  "${SOURCES}"
  PRECISIONS "s;d;c;z"
  )

### Generate the lib
add_library(pastix_bcsc
  ${generated_files}
  bcsc.c        # Not finished
  bcsc_fake.c   # Not finished
)

install(TARGETS pastix_bcsc
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)


