include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  z_debug_dump.h
  z_pastix_internal.h
  z_dynsched.h
  z_sopalin3d.h
  z_sopalin_init.h
  z_sopalin_option.h
  z_tools.h
  z_coefinit.h
  z_ooc.h
  # This C files is included in all other _he,_sy,_po,_ge
  z_sopalin3d.c
  z_raff_functions.c
  z_raff_functions.h
  z_raff_bicgstab.c
  z_raff.c
  z_raff_gmres.c
  z_raff_grad.c
  # included in pastix.c
  z_scaling.h
  z_scaling.c
  # Will be removed when we'll use kernels from kernels directory.
  z_compute_diag.c
  z_compute_diag.h
  z_compute_gemdm.c
  z_compute_trsm.c
  z_compute_trsm.h
  z_contrib.c
  z_sopalin_compute.c
  z_sopalin_compute.h
  # included in sopalin3d.c
  z_sopalin_sendrecv.c
  #z_dynsched.c
  z_contrib.c
  z_updo.c
  z_updo_sendrecv.c
  z_updo_sendrecv.h
  z_raff_pivot.c
  z_raff_gmres.c
  z_raff_grad.c
  z_raff_bicgstab.c
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(sopalin_headers
  DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_raff_functions_he.c
  z_raff_functions_sy.c
  z_raff_functions_po.c
  z_raff_functions_ge.c
  z_sopalin3d_ge.c
  z_sopalin3d_he.c
  z_sopalin3d_sy.c
  z_sopalin3d_po.c
#  sopalin/compute_context_nbr.c
  z_coefinit.c
  z_debug_dump.c
#  sopalin/ooc.c
  z_pastix.c
  z_pastix_fortran.c
  z_sopalin_init.c
  z_sopalin_option.c
#  sopalin/sparse_gemm_cpu.c
  z_tools.c
#  sopalin/murge.c
#  sopalin/csc_intern_old.c
#  sopalin/compute_diag.c
#  sopalin/compute_gemdm.c
#  sopalin/compute_trsm.c
#  sopalin/contrib.c
#  sopalin/MURGE_GetLocalElementList.c
#  sopalin/murge_pastix_fortran.c
#  sopalin/raff.c
#  sopalin/scaling.c
#  sopalin/sopalin_compute.c
#  sopalin/sopalin_sendrecv.c
#  sopalin/updo.c
#  sopalin/updo_sendrecv.c
#  sopalin/variable_csc.c
#  sopalin/zgetrf_stapiv_gpu.c
#  sopalin/zhetrf_stapiv_gpu.c
#  sopalin/zpotrf_stapiv_gpu.c
#  sopalin/zsytrf_stapiv_gpu.c
#
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")


### Generate the lib
add_library(pastix_sopalin ${generated_files}
  sopalin_thread.c)

add_dependencies(pastix_sopalin sopalin_headers)
add_dependencies(pastix_sopalin include_headers)
add_dependencies(pastix_sopalin common_headers)

install(TARGETS pastix_sopalin
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


