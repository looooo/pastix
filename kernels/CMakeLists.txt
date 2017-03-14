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
include(AddSourceFiles)

# reset variables
set(generated_sources "")
set(generated_headers "")

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

### Generate the headers in all precisions
set(HEADERS
  pastix_zcores.h
  z_nan_check.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "s;d;c;z")

set(kernels_headers
  ${generated_headers}
  )

add_custom_target(kernels_headers_tgt
  DEPENDS ${kernels_headers} )

### Generate the sources in all precisions
set(SOURCES
  # CPU extra kernel
  core_zgemdm.c
  core_zgetro.c
  core_zgeadd.c
  core_zplrnt.c
  core_zgemmsp.c
  core_ztrsmsp.c
  # Kernels
  core_zsytrfsp.c
  core_zhetrfsp.c
  core_zpotrfsp.c
  core_zgetrfsp.c
  core_zgelrops.c
  core_zgelrops_SVD.c
  core_zgelrops_RRQR.c
)

precisions_rules_py(generated_sources
  "${SOURCES}"
  PRECISIONS "s;d;c;z")

set(kernels_sources
  ${generated_sources}
  )

if(PASTIX_WITH_CUDA)
   set( CUDA_SOURCES
     # GPU kernel
     gpu_zgemmsp.c
    )
   set(CUDA_generated_files "")
   precisions_rules_py(CUDA_generated_files
        "${CUDA_SOURCES}"
        PRECISIONS "s;d;c;z")

   set(kernels_sources
       ${kernels_sources}
       ${CUDA_generated_files})

   add_subdirectory(gpus/kepler)

endif(PASTIX_WITH_CUDA)

add_library(pastix_kernels
  ${kernels_sources}
  )

add_dependencies(pastix_kernels
  kernels_headers_tgt
)

target_link_libraries( pastix_kernels
  ${LAPACKE_LIBRARIES}
  ${LAPACK_SEQ_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  ${EXTRA_LIBS}
)

if(PASTIX_WITH_CUDA)
  target_link_libraries( pastix_kernels
    pastix_cucores_sm35
    )
endif()


### Add documented files to the global property
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  ${generated_sources}
  )

# add_documented_files(
#   # Headers
#   # Source files
#   )
