include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

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


