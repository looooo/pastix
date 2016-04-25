include(RulesPrecisions)
include(array2d)

# reset variables
set(generated_files "")
set(generated_headers "")

# Headers
set(HEADERS
  testing_zmain.h
  )
precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "s;d;c;z")

# Add testing files
set(SOURCES
  testing_zrradd.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

add_library(testing ${generated_files} ${generated_headers})

include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}/")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/")

link_directories(
  ${LAPACKE_LIBRARY_DIRS}
  ${LAPACK_LIBRARY_DIRS}
  ${CBLAS_LIBRARY_DIRS}
  ${BLAS_LIBRARY_DIRS}
)

message(STATUS "LAPACKE_LIBRARIES:    ${LAPACKE_LIBRARIES}")
message(STATUS "LAPACK_SEQ_LIBRARIES: ${LAPACK_SEQ_LIBRARIES}")
message(STATUS "CBLAS_LIBRARIES:      ${CBLAS_LIBRARIES}")
message(STATUS "BLAS_SEQ_LIBRARIES:   ${BLAS_SEQ_LIBRARIES}")

set(PASTIX_TESTINGS_LIBRARIES
  testing
  pastix_kernels
  pastix
  ${LAPACKE_LIBRARIES}
  ${LAPACK_SEQ_LIBRARIES}
  ${BLAS_SEQ_LIBRARIES}
  ${EXTRA_LIBS}
  m
)

# Build four precision executables
set_source_files_properties(
  testing_zmain.c
  PROPERTIES DEPENDS testing_headers
  )

set(EXECUTABLE
  testing_zmain.c
)

precisions_rules_py(generated_files
  "${EXECUTABLE}"
  PRECISIONS "s;d;c;z")

foreach(prec s d c z)
  add_executable(testing_${prec}main testing_${prec}main.c)
  target_link_libraries(testing_${prec}main ${PASTIX_TESTINGS_LIBRARIES})
endforeach()
