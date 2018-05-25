###
#
#  @copyright 2013-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @date 2013-06-24
#
###
include(RulesPrecisions)
include(AddSourceFiles)

# reset variables
set(generated_files   "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

### Generate the headers in all precisions
set(HEADERS
  z_bcsc.h
)

### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "p;s;d;c;z")

set(bcsc_headers
  ${generated_headers}
  bcsc.h
  )

add_custom_target(bcsc_headers_tgt
  DEPENDS ${bcsc_headers} )

### Generate the sources in all precisions
set(SOURCES
  bcsc_zinit.c          # Not finished
  bcsc_zcompute.c       # Not finished
  z_bcsc_norm.c         # Not finished
  z_bcsc_matrixvector.c # Not finished
)

precisions_rules_py(generated_sources
  "${SOURCES}"
  PRECISIONS "s;d;c;z")

### Generate the lib
add_library(pastix_bcsc
  ${generated_sources}
  bcsc.c        # Not finished
)

add_dependencies( pastix_bcsc
  bcsc_headers_tgt )

install(TARGETS pastix_bcsc
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

# Not commented yet
#
# ### Add documented files to the global property
# add_documented_files(
#   DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#   ${generated_headers}
#   ${generated_sources}
#   )
#
# add_documented_files(
#   # Headers
#   bcsc.h
#   # Source files
#   bcsc.c
#   bcsc_fake.c
#   )
