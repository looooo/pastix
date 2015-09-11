include(RulesPrecisions)

# reset variables
set(generated_files "")
set(generated_headers "")

set(HEADERS
   z_tools.h
#  z_debug_dump.h
#  z_pastix_internal.h
#  z_dynsched.h
#  z_sopalin3d.h
#  z_sopalin_init.h
#  z_sopalin_option.h
#  z_coefinit.h
#  z_ooc.h
#  # This C files is included in all other _he,_sy,_po,_ge
#  z_sopalin3d.c
#  z_raff_functions.c
#  z_raff_functions.h
#  z_raff_bicgstab.c
#  z_raff.c
#  z_raff_gmres.c
#  z_raff_grad.c
#  # included in pastix.c
#  z_scaling.h
#  z_scaling.c
#  # Will be removed when we'll use kernels from kernels directory.
#  z_compute_diag.c
#  z_compute_diag.h
#  z_compute_gemdm.c
#  z_compute_trsm.c
#  z_compute_trsm.h
#  z_contrib.c
#  z_sopalin_compute.c
#  z_sopalin_compute.h
#  # included in sopalin3d.c
#  z_sopalin_sendrecv.c
#  #z_dynsched.c
#  z_contrib.c
#  z_updo.c
#  z_updo_sendrecv.c
#  z_updo_sendrecv.h
#  z_raff_pivot.c
#  z_raff_gmres.c
#  z_raff_grad.c
#  z_raff_bicgstab.c
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(sopalin_headers
  DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
#set(SOURCES
#)
#
#precisions_rules_py(generated_files
#   "${SOURCES}"
#   PRECISIONS "s;d;c;z")


### Generate the lib
#add_library(pastix_sopalin ${generated_files}
#  sopalin_thread.c)

#add_dependencies(pastix_sopalin sopalin_headers)
#add_dependencies(pastix_sopalin include_headers)
#add_dependencies(pastix_sopalin common_headers)
#
#install(TARGETS pastix_sopalin
#  RUNTIME DESTINATION bin
#  ARCHIVE DESTINATION lib
#  LIBRARY DESTINATION lib)


