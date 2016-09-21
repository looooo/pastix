set (EXAMPLES
  analyze.c
  simple.c
  bench_facto.c)

foreach (_file ${EXAMPLES})
  get_filename_component(_name_we ${_file} NAME_WE)
  add_executable(${_name_we} ${_file})
  target_link_libraries(${_name_we} pastix)
endforeach()
