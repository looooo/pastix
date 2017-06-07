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
  starpu_sparse_matrix.c
)

set(HEADERS
  pastix_starpu.h
)

### Generate the pastix_parsec wrappers for all required precisions
set(SOURCES
  codelet_zgemmsp.c
  codelet_zpotrfsp.c
  starpu_zpotrf.c
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

add_library(pastix_starpu
  ${generated_headers}
  ${generated_sources}
  ${EXTRA_SOURCES})

target_link_libraries(pastix_starpu
  pastix_kernels
  ${STARPU_LIBRARIES}
  ${STARPU_EXTRA_LIBRARIES}
  )

add_dependencies(pastix_starpu
  sopalin_headers
)

### Generate the lib
if (MPI_FOUND)
  set_target_properties(pastix_starpu PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif (MPI_FOUND)

install(TARGETS pastix_starpu
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

### Add non-generated documented files to the global property
add_documented_files(
  # Headers
  pastix_starpu.h
  # Source files
  starpu_sparse_matrix.c
  )
