add_subdirectory(old)
set (OLD_EXAMPLES
  analyze.c
  simple.c
  reentrant.c
  bench_facto.c
  step-by-step.c
  )

foreach (_file ${OLD_EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix)
endforeach()
