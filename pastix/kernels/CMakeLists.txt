include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  pastix_zcores.h
  sparse_zgemm_fermi.h
  gemm_stencil.h
  gemdm_stencil.h
  zgemm_fermi_define.h
  zgemm_fermi_define_bottom.h
  zgemm_fermi_define_corner.h
  zgemm_fermi_define_right.h
  zgemdm_fermi_define.h
  zgemdm_fermi_define_bottom.h
  zgemdm_fermi_define_corner.h
  zgemdm_fermi_define_right.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

list(APPEND generated_headers
  gemm_stencil_bottom.h
  gemdm_stencil_bottom.h
  gemm_stencil_right.h
  gemdm_stencil_right.h
  gemm_stencil_corner.h
  gemdm_stencil_corner.h)

add_custom_target(kernels_headers DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  # extra kernel
  core_zgemdm.c
  core_zgetro.c
  core_zgeadd.c
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
    # CUDA
    sparse_zgemm_fermi.cu
    sparse_zgemdm_fermi.cu
    )
  precisions_rules_py(CUDA_generated_files
    "${CUDA_SOURCES}"
    PRECISIONS "s;d;c;z")
endif(PASTIX_WITH_CUDA)

if(PASTIX_WITH_CUDA)
  foreach( _cudafile ${CUDA_generated_files} )
    get_source_file_property(_IsInBinaryDir ${_cudafile} IS_IN_BINARY_DIR )
    if( _IsInBinaryDir )
      string(REGEX REPLACE 
        ${_cudafile}
        ${CMAKE_CURRENT_BINARY_DIR}/${_cudafile}
        CUDA_generated_files "${CUDA_generated_files}")
    endif( _IsInBinaryDir )
    MESSAGE(STATUS "CUDA_generated_files ${CUDA_generated_files}")
  endforeach()
  list(APPEND generated_files ${CUDA_generated_files})
  list(APPEND CUDA_NVCC_FLAGS "-DCUDA_SM_VERSION=20  -arch sm_20")
  ### Generate the lib
  cuda_add_library(pastix_kernels ${generated_files} )
else(PASTIX_WITH_CUDA)
  ### Generate the lib
  add_library(pastix_kernels ${generated_files} )
endif(PASTIX_WITH_CUDA)


add_dependencies(pastix_kernels kernels_headers)
add_dependencies(pastix_kernels include_headers)
add_dependencies(pastix_kernels common_headers)

include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
link_directories(
  ${CBLAS_LIBRARY_DIRS}
  ${BLAS_LIBRARY_DIRS}
)

set(TESTING_KERNELS_LIBRARIES
  pastix
  pastix_kernels
  ${CBLAS_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  ${EXTRA_LIBS}
  m
)

set(SOURCES test_zgetrfsp.c)
if (PASTIX_WITH_CUDA)
  list(APPEND SOURCES test_sparse_zgemm.c)
endif(PASTIX_WITH_CUDA)
precisions_rules_py(generated_tests
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

foreach(example ${generated_tests})
  GET_FILENAME_COMPONENT(output ${example} NAME_WE)
  add_executable(${output} ${example})
  target_link_libraries( ${output}
    ${TESTING_KERNELS_LIBRARIES} )
  if (PASTIX_WITH_MPI)
    target_link_libraries( ${output} ${MPI_LIBRARIES} )
  endif(PASTIX_WITH_MPI)
  if (PASTIX_WITH_CUDA)
    target_link_libraries(${output} ${CUDA_LIBRARIES})
    target_link_libraries(${output} ${CUDA_CUBLAS_LIBRARIES})
  endif(PASTIX_WITH_CUDA)
endforeach()

add_test(test_sgetrfsp test_sgetrfsp)
add_test(test_dgetrfsp test_dgetrfsp)
add_test(test_cgetrfsp test_cgetrfsp)
add_test(test_zgetrfsp test_zgetrfsp)
if (PASTIX_WITH_CUDA)
  add_test(test_sparse_sgemm_2000_2000_200 test_sparse_sgemm 2000 2000 20)
  add_test(test_sparse_dgemm_2000_2000_200 test_sparse_dgemm 2000 2000 20)
  add_test(test_sparse_cgemm_2000_2000_200 test_sparse_cgemm 2000 2000 20)
  add_test(test_sparse_zgemm_2000_2000_200 test_sparse_zgemm 2000 2000 20)
endif(PASTIX_WITH_CUDA)
install(TARGETS pastix_kernels
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


