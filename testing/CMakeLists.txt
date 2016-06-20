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
  testing_zgradd.c
  testing_zlrm2.c
  testing_zlrmm.c
  testing_zlrmge.c
  testing_zge2lr.c
)

precisions_rules_py(generated_files
   "${SOURCES}"
   PRECISIONS "s;d;c;z")

add_library(pastix_testings ${generated_files} ${generated_headers})

find_package(TMG)
if(TMG_FOUND)
  include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
  include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}/")
  include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/")

  link_directories(
    ${LAPACKE_LIBRARY_DIRS}
    ${LAPACK_LIBRARY_DIRS}
    ${TMG_LIBRARY_DIRS}
    ${CBLAS_LIBRARY_DIRS}
    ${BLAS_LIBRARY_DIRS}
    )

  message(STATUS "LAPACKE_LIBRARIES:    ${LAPACKE_LIBRARIES}")
  message(STATUS "TMG_LIBRARIES: ${TMG_LIBRARIES}")
  message(STATUS "LAPACK_SEQ_LIBRARIES: ${LAPACK_SEQ_LIBRARIES}")
  message(STATUS "CBLAS_LIBRARIES:      ${CBLAS_LIBRARIES}")
  message(STATUS "BLAS_SEQ_LIBRARIES:   ${BLAS_SEQ_LIBRARIES}")

  set(PASTIX_TESTINGS_LIBRARIES
    pastix_testings
    pastix_kernels
    pastix
    ${LAPACKE_LIBRARIES}
    ${TMG_LIBRARIES}
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
    testing_zrun.c
    )

  precisions_rules_py(generated_files
    "${EXECUTABLE}"
    PRECISIONS "s;d;c;z")

  foreach(prec s d c z)
    add_executable(testing_${prec}main testing_${prec}main.c)
    target_link_libraries(testing_${prec}main ${PASTIX_TESTINGS_LIBRARIES})
    add_executable(testing_${prec}run testing_${prec}run.c)
    target_link_libraries(testing_${prec}run ${PASTIX_TESTINGS_LIBRARIES})
  endforeach()
else()
  message(STATUS "libtmg not found -- Can't compile testings")
endif()
