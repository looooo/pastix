include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  starpu_zdefines.h
  starpu_zkernels.h
  starpu_zsubmit.h
  starpu_zregister_data.h
  #starpu_zkernels.h
  starpu_zsubmit_tasks.h
  starpu_zupdo.h
  starpu_zupdo_kernels.h
  # included in _ge, _he, _sy, _po files
  starpu_zkernels.c
  starpu_zsubmit_tasks.c
  starpu_zupdo.c
  starpu_zupdo_kernels.c
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(starpu_headers DEPENDS ${generated_headers})

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  starpu_zgetrfsp_cpu.c
  starpu_zpotrfsp_cpu.c
  starpu_zsytrfsp_cpu.c
  starpu_zhetrfsp_cpu.c
  starpu_zregister_data.c
  starpu_zsubmit.c
  starpu_zkernels_ge.c
  starpu_zkernels_he.c
  starpu_zkernels_sy.c
  starpu_zkernels_po.c
  starpu_zsubmit_tasks_ge.c
  starpu_zsubmit_tasks_he.c
  starpu_zsubmit_tasks_sy.c
  starpu_zsubmit_tasks_po.c
  starpu_zupdo_ge.c
  starpu_zupdo_he.c
  starpu_zupdo_sy.c
  starpu_zupdo_po.c
  starpu_zupdo_kernels_ge.c
  starpu_zupdo_kernels_he.c
  starpu_zupdo_kernels_sy.c
  starpu_zupdo_kernels_po.c
 )

if (PASTIX_WITH_CUDA)
  list(APPEND SOURCES starpu_zsytrfsp_cuda.c)
endif(PASTIX_WITH_CUDA)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

### Generate the lib
add_library(pastix_starpu
  ${generated_files}
  starpu_callback.c
)

add_dependencies(pastix_starpu starpu_headers)
add_dependencies(pastix_starpu include_headers)
add_dependencies(pastix_starpu common_headers)

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

