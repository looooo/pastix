include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  z_pastix.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(include_headers
  DEPENDS ${generated_headers} )


#install(TARGETS 
#  RUNTIME DESTINATION bin
#  ARCHIVE DESTINATION lib
#  LIBRARY DESTINATION lib)


