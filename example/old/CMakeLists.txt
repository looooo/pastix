set (EXAMPLES
  simple.c
  step-by-step.c
  schur.c
  )

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable("old_${_name_we}" ${_file})
  target_link_libraries("old_${_name_we}" pastix)

  install(TARGETS old_${_name_we} RUNTIME DESTINATION examples/old )
  install(FILES   ${_file}            DESTINATION examples/old )

endforeach()
