include(RulesPrecisions)
include(AddSourceFiles)

### reset variables
set(generated_sources "")
set(generated_headers "")

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

### Generate the headers in all precisions
set(HEADERS
  z_raff_functions.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "p;s;d;c;z")

set(refinement_headers
  ${generated_headers}
  )

add_custom_target(refinement_headers_tgt
  DEPENDS ${refinement_headers} )

### Generate the sources in all precisions
set(SOURCES
  z_raff_functions.c
  z_raff_gmres.c
  z_raff_grad.c
  z_raff_pivot.c
  z_raff_bicgstab.c
)

precisions_rules_py(generated_sources
  "${SOURCES}"
  PRECISIONS "p;s;d;c;z"
  TARGETDIR "${CMAKE_CURRENT_BINARY_DIR}")

set(refinement_sources
  ${generated_sources}
  pastix_task_raff.c
)

# ### Add documented files to the global property
# add_documented_files(
#   DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#   ${generated_headers}
#   ${generated_sources}
#   )

add_documented_files(
  # Headers
  # Source files
  pastix_task_raff.c
  )
