###
#
#  @copyright 2013-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2013-06-24
#
###
include(RulesPrecisions)
include(AddSourceFiles)

### reset variables
set(generated_sources "")
set(generated_headers "")

include_directories("${CMAKE_CURRENT_BINARY_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/drivers")

### Generate the headers in all precisions
set(HEADERS
  z_spm.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "p;s;d;c;z")

set(spm_headers
  ${generated_headers}
  spm.h
  spm_drivers.h
  )

add_custom_target(spm_headers_tgt
  DEPENDS ${spm_headers} )

### Generate the sources in all precisions
set(SOURCES
  z_spm.c
  z_spm_2dense.c
  z_spm_dof_extend.c
  z_spm_norm.c
  z_spm_scal.c

  z_spm_convert_to_csc.c
  z_spm_convert_to_csr.c
  z_spm_convert_to_ijv.c
  z_spm_expand.c
  z_spm_genrhs.c
  z_spm_integer.c
  z_spm_laplacian.c
  z_spm_matrixvector.c
  z_spm_print.c
  )

precisions_rules_py(generated_sources
  "${SOURCES}"
  PRECISIONS "p;s;d;c;z")

set(spm_sources
  ${generated_sources}
  spm.c
  spm_io.c
  spm_integers.c
  spm_dof_extend.c
  spm_read_driver.c
  spm_gen_fake_values.c
  drivers/skitf.f
  drivers/iohb.c
  drivers/mmio.c
  drivers/laplacian.c
  drivers/readhb.c
  drivers/readijv.c
  drivers/readmm.c
  drivers/readrsa.c
  )

add_library(pastix_spm
  ${spm_sources}
  )

target_link_libraries(pastix_spm
  ${LAPACKE_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${CBLAS_LIBRARIES}
  )

add_dependencies(pastix_spm
  spm_headers_tgt
)

### Generate the lib
if (MPI_FOUND)
  set_target_properties(pastix_spm PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif (MPI_FOUND)

install(TARGETS pastix_spm
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES spm.h
  DESTINATION include )

### Add documented files to the global property
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  ${generated_sources}
  )

add_documented_files(
  # Headers
  spm.h
  #spm_drivers.h
  # Source files
  spm.c
  spm_io.c
  spm_read_driver.c
  spm_dof_extend.c
  spm_integers.c
  )
