include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../common")

add_executable(genheader
  genheader.c)
add_dependencies(genheader include_headers)
add_dependencies(genheader common_headers)
set(INCDIR
  ${CMAKE_CURRENT_SOURCE_DIR}/../include )

# add_custom_command(
#   OUTPUT  ${INCDIR}/pastix.h
#   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/genheader ${INCDIR}/pastix.h unused unused unused unused unused unused C
#   COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../sopalin/pastix.h >> ${INCDIR}/pastix.h
#   COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../common/api.h     >> ${INCDIR}/pastix.h
#   DEPENDS
#   genheader
#   ${CMAKE_CURRENT_SOURCE_DIR}/../sopalin/pastix.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/../common/api.h
# )

# add_custom_command(
#   OUTPUT  ${INCDIR}/pastix_fortran.h
#   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/genheader unused ${INCDIR}/pastix_fortran.h unused unused unused unused unused Fortran
#   COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../common/pastix_fortran.h >> ${INCDIR}/pastix_fortran.h
#   DEPENDS
#   genheader
#   ${CMAKE_CURRENT_SOURCE_DIR}/../common/pastix_fortran.h
# )

# add_custom_command(
#   OUTPUT  ${INCDIR}/murge.inc
#   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/genheader unused unused ${INCDIR}/murge.inc ${CMAKE_CURRENT_SOURCE_DIR}/../sopalin/murge_pastix.inc ${CMAKE_CURRENT_SOURCE_DIR}/../murge/include/murge.h ${CMAKE_CURRENT_SOURCE_DIR}/../murge/scripts/genfort.pl ${CMAKE_CURRENT_SOURCE_DIR}/insert-fortran-inc.sh Murge
#   COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../common/pastix_fortran.inc >> ${INCDIR}/murge.inc
#   DEPENDS
#   genheader
#   ${CMAKE_CURRENT_SOURCE_DIR}/../common/pastix_fortran.inc
#   ${CMAKE_CURRENT_SOURCE_DIR}/../sopalin/murge_pastix.inc
#   ${CMAKE_CURRENT_SOURCE_DIR}/../murge/include/murge.h
#   ${CMAKE_CURRENT_SOURCE_DIR}/../murge/scripts/genfort.pl
#   insert-fortran-inc.sh
# )

# add_custom_command(
#   OUTPUT  ${INCDIR}/pastix_nompi.h
#   COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/../common/nompi.h ${INCDIR}/pastix_nompi.h
#   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../common/nompi.h
# )

# add_custom_target(headers ALL
#   DEPENDS
#   ${INCDIR}/pastix.h
#   ${INCDIR}/pastix_fortran.h
#   ${INCDIR}/pastix_nompi.h
# #  ${INCDIR}/murge.inc
# )
