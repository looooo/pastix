include(RulesPrecisions)
# reset variables
set(generated_files "")
set(generated_headers "")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  z_solver.h
  z_updown.h
  z_ftgt.h
)
### generate the dsparse_cores headers for all possible precisions
precisions_rules_py(generated_headers
                    "${HEADERS}"
                    PRECISIONS "s;d;c;z")

add_custom_target(blend_headers 
  DEPENDS ${generated_headers} )

### Generate the dsparse wrappers for all required precisions
set(SOURCES
  z_pastix_task_analyze.c
  z_solver.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

list(APPEND generated_files
  cand.c
  elimin_graph.c
  elimin_tree.c
  propmap.c
  extracblk.c
  splitsymbol.c
  blendctrl.c

  blend.c
  blend_symbol_cost.c
  bulles.c
  cost.c
  costfunc.c
  simu_run.c
#  distribPart.c
  extendVector.c
  fanboth2.c
  queue.c
  simu.c
#  partbuild.c
#  smart_cblk_split.c
  solverMatrixGen.c
  solverRealloc.c
  solver_check.c
  solver_io.c
#  splitpartlocal.c
  task.c
  blend_distributeOnGPU.c
)

### Generate the lib
add_library(pastix_blend ${generated_files} )

add_dependencies(pastix_blend
  blend_headers
)

install(TARGETS pastix_blend
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)


