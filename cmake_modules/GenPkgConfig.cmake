###
#
# @copyright (c) 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                          Univ. Bordeaux. All rights reserved.
#
###
#
#  @file GenPkgConfig.cmake
#
#  @project PaStiX
#
#  @version 6.0.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 10-11-2014
#
###

###
#
# convert_incstyle_to_pkgconfig():
# convert a libraries list to follow the pkg-config style used in
# CLEAN_LIB_LIST
#
###
macro(convert_incstyle_to_pkgconfig _inclist)
  set(${_inclist}_CPY "${${_inclist}}")
  set(${_inclist} "")
  foreach(_dep ${${_inclist}_CPY})
    if (${_dep} MATCHES "^-D")
      list(APPEND ${_inclist} ${_dep})
    else()
      list(APPEND ${_inclist} "-I${_dep}")
    endif()
  endforeach()
endmacro(convert_incstyle_to_pkgconfig)

###
#
# convert_libstyle_to_pkgconfig():
# convert a libraries list to follow the pkg-config style used in
# CLEAN_LIB_LIST
#
###
macro(convert_libstyle_to_pkgconfig _liblist)
  set(${_liblist}_CPY "${${_liblist}}")
  set(${_liblist} "")
  foreach(_dep ${${_liblist}_CPY})
    if (${_dep} MATCHES "^/")
      get_filename_component(dep_libname ${_dep} NAME)
      get_filename_component(dep_libdir  ${_dep} PATH)
      string(REPLACE ".so"    "" dep_libname "${dep_libname}")
      string(REPLACE ".a"     "" dep_libname "${dep_libname}")
      string(REPLACE ".dylib" "" dep_libname "${dep_libname}")
      string(REPLACE ".dll"   "" dep_libname "${dep_libname}")
      string(REPLACE "lib"    "" dep_libname "${dep_libname}")
      list(APPEND ${_liblist} -L${dep_libdir} -l${dep_libname})
    elseif(NOT ${_dep} MATCHES "^-")
      list(APPEND ${_liblist} "-l${_dep}")
    else()
      list(APPEND ${_liblist} ${_dep})
    endif()
  endforeach()
endmacro(convert_libstyle_to_pkgconfig)

###
#
# clean_lib_list: clean libraries lists to follow the pkg-config style
#                 used in GENERATE_PKGCONFIG_FILE
#
###
macro(clean_lib_list _package)
    if ( ${_package}_PKGCONFIG_INCS )
      list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_INCS)
      convert_incstyle_to_pkgconfig(${_package}_PKGCONFIG_INCS)
      string(REPLACE ";" " " ${_package}_PKGCONFIG_INCS "${${_package}_PKGCONFIG_INCS}")
    endif()
    if ( ${_package}_PKGCONFIG_LIBS )
      list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_LIBS)
      convert_libstyle_to_pkgconfig(${_package}_PKGCONFIG_LIBS)
      string(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS "${${_package}_PKGCONFIG_LIBS}")
    endif()
    if ( ${_package}_PKGCONFIG_LIBS_PRIVATE )
      list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_LIBS_PRIVATE)
      convert_libstyle_to_pkgconfig(${_package}_PKGCONFIG_LIBS_PRIVATE)
      string(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS_PRIVATE "${${_package}_PKGCONFIG_LIBS_PRIVATE}")
    endif()
    if ( ${_package}_PKGCONFIG_REQUIRED )
      list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED)
      string(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED "${${_package}_PKGCONFIG_REQUIRED}")
    endif()
    if ( ${_package}_PKGCONFIG_REQUIRED_PRIVATE )
      list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED_PRIVATE)
      string(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED_PRIVATE "${${_package}_PKGCONFIG_REQUIRED_PRIVATE}")
    endif()
endmacro(clean_lib_list)

###
#
# generate_pkgconfig_file: generate file pastix.pc
#
###
macro(generate_pkgconfig_file)

    # The link flags specific to this package and any required libraries
    # that don't support PkgConfig
    set(PASTIX_PKGCONFIG_LIBS pastix pastix_bcsc pastix_spm pastix_kernels)

    # The link flags for private libraries required by this package but not
    # exposed to applications
    set(PASTIX_PKGCONFIG_LIBS_PRIVATE "")

    # A list of packages required by this package
    set(PASTIX_PKGCONFIG_REQUIRED "")

    # A list of private packages required by this package but not exposed to
    # applications
    set(PASTIX_PKGCONFIG_REQUIRED_PRIVATE "")

    if(PASTIX_WITH_STARPU)
      list(APPEND PASTIX_PKGCONFIG_LIBS pastix_starpu)
      if ( STARPU_FOUND_WITH_PKGCONFIG )
        if ( PASTIX_WITH_MPI )
          list(APPEND PASTIX_PKGCONFIG_REQUIRED starpumpi-${PASTIX_STARPU_VERSION})
        else()
          list(APPEND PASTIX_PKGCONFIG_REQUIRED starpu-${PASTIX_STARPU_VERSION})
        endif()
      else()
        list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE ${STARPU_LIBRARIES})
      endif()
    endif()

    # PaRSEC
    if(PASTIX_WITH_PARSEC)
      list(APPEND PASTIX_PKGCONFIG_LIBS pastix_parsec)
      if ( PARSEC_FOUND_WITH_PKGCONFIG )
        list(APPEND PASTIX_PKGCONFIG_REQUIRED parsec)
      else()
        list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE ${PARSEC_LIBRARIES})
      endif()
    endif()

    if(PASTIX_WITH_CUDA)
      list(APPEND PASTIX_PKGCONFIG_LIBS pastix_kernels_cuda)
    endif()

    if(PASTIX_WITH_EZTRACE)
      list(APPEND PASTIX_PKGCONFIG_REQUIRED eztrace litl)
    endif()

    list(APPEND PASTIX_PKGCONFIG_INCS
      )
    list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE
      ${LAPACKE_LIBRARIES}
      ${LAPACK_SEQ_LIBRARIES}
      ${CBLAS_LIBRARIES}
      ${BLAS_SEQ_LIBRARIES}
      )
    list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE
      ${EXTRA_LIBRARIES}
      )

    # HwLoc
    if ( HWLOC_FOUND_WITH_PKGCONFIG )
      list(APPEND PASTIX_PKGCONFIG_REQUIRED hwloc)
    else()
      list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE ${HWLOC_LIBRARIES})
    endif()

    # Define required package
    # -----------------------
    clean_lib_list(PASTIX)

    # Create .pc file
    # ---------------
    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/pastix.pc.in"
      "${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/pastix.pc" @ONLY)

    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/pastixf.pc.in"
      "${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/pastixf.pc" @ONLY)

    # installation
    # ------------
    install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/pastix.pc"
      "${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/pastixf.pc"
      DESTINATION lib/pkgconfig)

endmacro(generate_pkgconfig_file)

###
#
# generate_env_file: generate files pastix.pc
#
###
macro(generate_env_file)

    # Create .sh file
    # ---------------
    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/pastix_env.sh.in"
      "${CMAKE_CURRENT_BINARY_DIR}/bin/pastix_env.sh" @ONLY)

    # installation
    # ------------
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/bin/pastix_env.sh"
      DESTINATION bin)

endmacro(generate_env_file)

##
## @end file GenPkgConfig.cmake
##
