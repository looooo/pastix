include(RulesPrecisions)
include(array2d)

# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  zmurge.h
  zmurge_pastix.h
  zmurge_defines.h
  ZMURGE_GetLocalElementList.c
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(murge_headers
  DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  zmurge.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

if(PASTIX_INT64)
  set(INTSIZE 8)
  set(INTTYPE int64_t)
else(PASTIX_INT64)
  set(INTSIZE 4)
  set(INTTYPE int32_t)
endif(PASTIX_INT64)

set(name_prec_coefsize_cplx
  "murge"        "z" 8 1
  "murge"        "d" 8 0
  "murge"        "c" 4 1
  "murge"        "s" 4 0
  "murge_pastix" "z" 8 1
  "murge_pastix" "d" 8 0
  "murge_pastix" "c" 4 1
  "murge_pastix" "s" 4 0
  )
array2d_begin_loop( advanced "${name_prec_coefsize_cplx}"
  4 "name;prec;coefsize;cplx;" )
if (NOT advanced)
  message(FATAL_ERROR "advanced: ${advanced}")
endif(NOT advanced)
while( advanced )
  # genrate fortran interface file
  # include common.h then generate using perl script
  add_custom_command(OUTPUT ${prec}${name}_fortran.c
    COMMAND ${CMAKE_COMMAND}
    ARGS -E echo "#define PRECISION_${prec}" > ${prec}${name}_fortran.c
    COMMAND ${CMAKE_COMMAND}
    ARGS -E echo "#include \"common.h\"" >> ${prec}${name}_fortran.c
    COMMAND ${PERL_EXECUTABLE}
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/geninterface.pl -n -f ${CMAKE_CURRENT_BINARY_DIR}/${prec}${name}.h  >> ${prec}${name}_fortran.c
    DEPENDS geninterface.pl ${CMAKE_CURRENT_BINARY_DIR}/${prec}${name}.h
    VERBATIM)
  add_custom_target(${prec}${name}_fortran ALL DEPENDS ${prec}${name}_fortran.c)
  set(generated_files ${generated_files} ${prec}${name}_fortran.c)

  # generate .inc file
  if ("murge_pastix" STREQUAL ${name})
    # generate interface using genfort.pl, add API variables using APIc2f.pl
    # and remove MURGE_XXXX_KIND definition which are already in murge.inc
    add_custom_command(OUTPUT ${prec}${name}.inc
      COMMAND ${PERL_EXECUTABLE}
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/genfort.pl -f  ${CMAKE_CURRENT_BINARY_DIR}/${prec}${name}.h -s ${INTSIZE} -l ${INTSIZE} -r ${coefsize} -c ${cplx} > ${prec}${name}.inc.tmp
      COMMAND ${PERL_EXECUTABLE}
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../tools/APIc2f.pl -m -f ${CMAKE_CURRENT_SOURCE_DIR}/../include/pastix/api.h >> ${prec}${name}.inc.tmp
      COMMAND ${PERL_EXECUTABLE}
      ARGS -ne 'print if !(/ "PARAMETER.*_KIND"/)' ${prec}${name}.inc.tmp > ${prec}${name}.inc
      COMMAND ${CMAKE_COMMAND}
      ARGS -E remove ${prec}${name}.inc.tmp
      DEPENDS
      ../tools/APIc2f.pl genfort.pl
      ${CMAKE_CURRENT_BINARY_DIR}/${prec}${name}.h
      ${CMAKE_CURRENT_SOURCE_DIR}/../include/pastix/api.h
      )
  else()
    # generate interface using genfort.pl
    add_custom_command(OUTPUT ${prec}${name}.inc
      COMMAND ${PERL_EXECUTABLE}
      ARGS ${CMAKE_CURRENT_SOURCE_DIR}/genfort.pl -f  ${CMAKE_CURRENT_BINARY_DIR}/${prec}${name}.h -s ${INTSIZE} -l ${INTSIZE} -r ${coefsize} -c ${cplx} > ${prec}${name}.inc
      DEPENDS genfort.pl ${CMAKE_CURRENT_BINARY_DIR}/${prec}${name}.h
      )
  endif()
  add_custom_target(${prec}${name} ALL DEPENDS ${prec}${name}.inc)

  install(FILES
    ${prec}${name}.inc
    ${prec}${name}.h
    DESTINATION include)

  array2d_advance()
endwhile(advanced)


### Generate the lib
add_library(pastix_murge ${generated_files})

install(TARGETS pastix_murge
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


