###
#
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.1
#  @author Mathieu Faverge
#  @author Ian Masliah
#  @author Pierre Ramet
#  @author Tony Delarue
#  @author Vincent Bridonneau
#  @author Alycia Lisito
#  @author Florent Pruvost
#  @author Nolan Bredel
#  @author Tom Moenne-Loccoz
#  @date 2023-10-25
#
###
include(RulesPrecisions)

# reset variables
set(generated_sources "")
set(generated_headers "")

set(EXTRA_SOURCES
  starpu.c
  starpu_tags.c
  starpu_sparse_matrix.c
  starpu_dense_matrix.c
  starpu_profile.c
  pastix_starpu_interface.c
)

### Generate the headers in all precisions
set(HEADERS
  pastix_zstarpu.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "s;d;c;z")

set(starpu_headers
  ${generated_headers}
  )

add_custom_target(starpu_headers_tgt
  DEPENDS ${starpu_headers} pastix_starpu.h )

### Generate the pastix_parsec wrappers for all required precisions
set(SOURCES
  codelet_blok_zgemmsp.c
  codelet_blok_zgetrfsp.c
  codelet_blok_zhetrfsp.c
  codelet_blok_zpotrfsp.c
  codelet_blok_zpxtrfsp.c
  codelet_blok_zscalo.c
  codelet_blok_zsytrfsp.c
  codelet_blok_ztrsmsp.c
  codelet_cblk_zgemmsp.c
  codelet_cblk_zgetrfsp.c
  codelet_cblk_zhetrfsp.c
  codelet_cblk_zpotrfsp.c
  codelet_cblk_zpxtrfsp.c
  codelet_cblk_zsytrfsp.c
  codelet_solve_zdiag.c
  codelet_solve_zgemm.c
  codelet_solve_ztrsm.c
  starpu_zdiag.c
  starpu_zgetrf.c
  starpu_zhetrf.c
  starpu_zpotrf.c
  starpu_zpxtrf.c
  starpu_zsytrf.c
  starpu_ztrsm.c
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

target_include_directories( pastix_starpu PRIVATE
  $<BUILD_INTERFACE:${PASTIX_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PASTIX_BINARY_DIR}>
  $<BUILD_INTERFACE:${PASTIX_SOURCE_DIR}/common> )

target_include_directories( pastix_starpu PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

target_link_libraries( pastix_starpu PUBLIC
  pastix_kernels
  MORSE::STARPU
  )

add_dependencies(pastix_starpu
  sopalin_headers
  starpu_headers_tgt
)

### Generate the lib
if (MPI_C_FOUND)
  set_target_properties(pastix_starpu PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif (MPI_C_FOUND)

# export target pastix_starpu
# ---------------------------
install(EXPORT pastix_starpuTargets
  NAMESPACE PASTIX::
  DESTINATION lib/cmake/pastix
  )

# Installation
# ------------
install(TARGETS pastix_starpu
        EXPORT pastix_starpuTargets
        ARCHIVE       DESTINATION lib
        LIBRARY       DESTINATION lib
        RUNTIME       DESTINATION bin
        PUBLIC_HEADER DESTINATION include)

### Add non-generated documented files to the global property
add_documented_files(
  # Headers
  pastix_starpu.h
  # Source files
  starpu.c
  starpu_sparse_matrix.c
  starpu_tags.c
  starpu_dense_matrix.c
  pastix_starpu_interface.c
  )
