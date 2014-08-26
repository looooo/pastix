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


if(PASTIX_INT64)
  set(INTSIZE 8)
else(PASTIX_INT64)
  set(INTSIZE 4)
endif(PASTIX_INT64)

# generate pastix_fortran.h file
# insert defines, add API variables using APIc2f.pl
add_custom_command(OUTPUT pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define PASTIX_INT_KIND    ${INTSIZE}" > pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define pastix_int_t       INTEGER(kind=${INTSIZE})" >> pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define pastix_uint_t      unsigned INTEGER(kind=${INTSIZE})" >> pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define pastix_data_ptr_t  INTEGER(kind=${CMAKE_SIZEOF_VOID_P})" >> pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define MPI_PASTIX_INT     MPI_INTEGER${INTSIZE}" >> pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define PASTIX_MAJOR_VERSION  ${PASTIX_VERSION_MAJOR}" >> pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define PASTIX_VERSION_MAJOR  ${PASTIX_VERSION_MAJOR}" >> pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define PASTIX_VERSION_MINOR ${PASTIX_VERSION_MINOR}" >> pastix_fortran.h
  COMMAND ${CMAKE_COMMAND}
  ARGS -E echo "#define PASTIX_VERSION_MICRO ${PASTIX_VERSION_MICRO}" >> pastix_fortran.h
  COMMAND ${PERL_EXECUTABLE}
  ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../tools/APIc2f.pl -f ${CMAKE_CURRENT_SOURCE_DIR}/../include/pastix/api.h >> pastix_fortran.h
  DEPENDS ../tools/APIc2f.pl pastix/api.h
  VERBATIM
  )
add_custom_target(pastix_fortran_h ALL DEPENDS pastix_fortran.h)
install(FILES
  pastix_fortran.h
  DESTINATION include)


#install(TARGETS
#  RUNTIME DESTINATION bin
#  ARCHIVE DESTINATION lib
#  LIBRARY DESTINATION lib)


