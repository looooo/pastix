include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  z_pastixdata.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(common_headers
  DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_integer.c
)
precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")



### Generate the lib
add_library(pastix_common ${generated_files} )

#install(TARGETS pastix_common
#  RUNTIME DESTINATION bin
#  ARCHIVE DESTINATION lib
#  LIBRARY DESTINATION lib)


