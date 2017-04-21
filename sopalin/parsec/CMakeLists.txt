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
include(RulesJDF)
# reset variables
set(generated_sources "")
set(generated_jdf     "")
set(generated_headers "")

set(EXTRA_SOURCES
  parsec_sparse_matrix.c
)

set(HEADERS
  pastix_parsec.h
)

### Generate the pastix_parsec wrappers for all required precisions
set(SOURCES
  parsec_zpotrf.c
  parsec_zgetrf.c
)

precisions_rules_py(generated_sources
  "${SOURCES}"
  PRECISIONS "s;d;c;z")

### Add generated documented files to the global property before generating JDFs
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  ${generated_sources}
  )

### Generate .c files from .jdf for all required precisions
set(JDF
  zpotrf_sp1dplus.jdf
  zgetrf_sp1dplus.jdf
  #
  zpotrf_sp2d.jdf
  zgetrf_sp2d.jdf
)
precisions_rules_py(generated_jdf
  "${JDF}"
  PRECISIONS "s;d;c;z")

jdf_rules(generated_sources
	 "${generated_jdf}"
	)

add_library(pastix_parsec
  ${generated_headers}
  ${generated_sources}
  ${EXTRA_SOURCES})

target_link_libraries(pastix_parsec
  pastix_kernels
  ${PARSEC_LIBRARIES}
  parsec_distribution
  parsec_distribution_matrix
  ${PARSEC_EXTRA_LIBRARIES}
  ${PARSEC_LIBRARIES}
  )

add_dependencies(pastix_parsec
  sopalin_headers
)

### Generate the lib
if (MPI_FOUND)
  set_target_properties(pastix_parsec PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif (MPI_FOUND)

install(TARGETS pastix_parsec
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

### Add non-generated documented files to the global property
add_documented_files(
  # Headers
  pastix_parsec.h
  # Source files
  parsec_sparse_matrix.c
  )
