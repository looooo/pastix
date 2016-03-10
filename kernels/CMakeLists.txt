include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

set(HEADERS
  pastix_zcores.h
)

### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(kernels_headers DEPENDS ${generated_headers} )

### Generate the CPU kernels for all required precisions
set(SOURCES
  # CPU extra kernel
  core_zgemdm.c
  core_zgetro.c
  core_zgeadd.c
  core_zplrnt.c
  core_zgemmsp.c
  # Kernels
  core_zsytrfsp.c
  core_zhetrfsp.c
  core_zpotrfsp.c
  core_zgetrfsp.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

if(PASTIX_WITH_CUDA)
   set( CUDA_SOURCES
     # GPU kernel
     gpu_zgemmsp.c
    )
   set(CUDA_generated_files "")
   precisions_rules_py(CUDA_generated_files
        "${CUDA_SOURCES}"
        PRECISIONS "s;d;c;z")

   set(generated_files
       ${generated_files}
       ${CUDA_generated_files})

   add_subdirectory(gpus/kepler)

endif(PASTIX_WITH_CUDA)

add_library(pastix_kernels ${generated_files} )
add_dependencies(pastix_kernels kernels_headers)

target_link_libraries( pastix_kernels pastix_cucores_sm35 )
