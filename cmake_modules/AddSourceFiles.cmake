#
# See https://cmake.org/pipermail/cmake/2010-March/035992.html
#
function(add_source_files)
  set(oneValueArgs PROJECT_NAME DIRECTORY)
  cmake_parse_arguments(add_source_files "" "${oneValueArgs}" "" ${ARGN})

  get_property(is_defined GLOBAL PROPERTY ${PROJECT_NAME}_LIB_SRCS2 DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY ${PROJECT_NAME}_LIB_SRCS2
      BRIEF_DOCS "List of source files included in the project ${PROJECT_NAME}"
      FULL_DOCS "List of source files included in the project ${PROJECT_NAME}")
  endif()

  if(NOT ${add_source_files_DIRECTORY} STREQUAL "")
    set(_dir "${add_source_files_DIRECTORY}/")
  else()
    set(_dir "")
  endif()

  set(ARGN ${add_source_files_UNPARSED_ARGUMENTS})
  # make absolute paths
  set(SRCS)
  foreach(s IN LISTS ARGN)
    set(s "${_dir}${s}")
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY ${PROJECT_NAME}_LIB_SRCS2 "${SRCS}")
endfunction(add_source_files)

#
# See https://cmake.org/pipermail/cmake/2010-March/035992.html
#
function(add_documented_files)
  set(oneValueArgs PROJECT_NAME DIRECTORY)
  cmake_parse_arguments(add_documented_files "" "${oneValueArgs}" "" ${ARGN})

  get_property(is_defined GLOBAL PROPERTY ${PROJECT_NAME}_DOX_SRCS DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY ${PROJECT_NAME}_DOX_SRCS
      BRIEF_DOCS "List of documented source files"
      FULL_DOCS "List of source files to be included into the in-code documentation")
  endif()

  if(NOT ${add_documented_files_DIRECTORY} STREQUAL "")
    set(_dir "${add_documented_files_DIRECTORY}/")
  else()
    set(_dir "")
  endif()

  set(ARGN ${add_documented_files_UNPARSED_ARGUMENTS})
  # make absolute paths
  set(SRCS)
  foreach(s IN LISTS ARGN)
    set(s "${_dir}${s}")
    if(NOT IS_ABSOLUTE "${s}")
      get_filename_component(s "${s}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY ${PROJECT_NAME}_DOX_SRCS "${SRCS}")
endfunction(add_documented_files)


#
# See https://cmake.org/pipermail/cmake/2010-March/035992.html
#
function(add_binary_to_completion)
  cmake_parse_arguments(add_binary_to_completion "" "${oneValueArgs}" "" ${ARGN})

  get_property(is_defined GLOBAL PROPERTY ${PROJECT_NAME}_COMP_BINARIES DEFINED)
  if(NOT is_defined)
    define_property(GLOBAL PROPERTY ${PROJECT_NAME}_COMP_BINARIES
      BRIEF_DOCS "List of binaries candidate for bash auto-completion"
      FULL_DOCS "List of binaries candidate for bash auto-completion")
  endif()

  set(ARGN ${add_binary_to_completion_UNPARSED_ARGUMENTS})
  # make absolute paths
  set(SRCS)
  foreach(s IN LISTS ARGN)
    get_filename_component(s "${s}" NAME_WE)
    list(APPEND SRCS "${s}")
  endforeach()
  # append to global list
  set_property(GLOBAL APPEND PROPERTY ${PROJECT_NAME}_COMP_BINARIES "${SRCS}")
endfunction(add_binary_to_completion)
