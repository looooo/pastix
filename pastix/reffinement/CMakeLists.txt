include(RulesPrecisions)

# reset variables
set(generated_files   "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  z_raff_functions.h
)

### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "p;s;d;c;z")

add_custom_target(reffinement_headers
  DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
z_raff_functions.c
# z_raff_functions_ge.c
# z_raff_functions_he.c
# z_raff_functions_po.c
# z_raff_functions_sy.c
z_raff_gmres.c
# z_raff_grad.c
# z_raff_pivot.c
# z_raff_bicgstab.c
)

# precisions_rules_py(generated_files
#   "${SOURCES}"
#   PRECISIONS "s;d;c;z"
#   )

# set_source_files_properties( ${generated_files} PROPERTIES DEPENDS reffinement_headers )


