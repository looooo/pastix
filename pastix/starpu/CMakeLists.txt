include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  starpu_zkernels.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(starpu_headers ALL SOURCES
  ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  starpu_zgetrfsp_cpu.c
  starpu_zpotrfsp_cpu.c
  starpu_zsytrfsp_cpu.c
  starpu_zhetrfsp_cpu.c
 )

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

### Generate the lib
add_library(pastix_starpu
  ${generated_files}
)

add_dependencies(pastix_starpu
  kernel_headers
)

include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
link_directories(
  ${CBLAS_LIBRARY_DIRS}
  ${BLAS_LIBRARY_DIRS}
)

set(TESTING_KERNELS_LIBRARIES
  pastix_starpu
  ${CBLAS_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  ${EXTRA_LIBS}
  m
)

