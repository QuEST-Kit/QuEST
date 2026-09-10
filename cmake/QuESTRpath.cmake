# Each executable can live at a different depth, e.g. bin/examples/extended.
function(setup_quest_rpath target destination)
  if(APPLE)
    set(_origin "@loader_path")
  elseif(UNIX)
    set(_origin "$ORIGIN")
  else()
    return()
  endif()
  cmake_path(ABSOLUTE_PATH destination BASE_DIRECTORY "${CMAKE_INSTALL_PREFIX}" OUTPUT_VARIABLE _from)
  set(_libdir "${CMAKE_INSTALL_LIBDIR}")
  cmake_path(ABSOLUTE_PATH _libdir BASE_DIRECTORY "${CMAKE_INSTALL_PREFIX}" OUTPUT_VARIABLE _to)
  file(RELATIVE_PATH _relative "${_from}" "${_to}")
  set_target_properties(${target} PROPERTIES
    BUILD_RPATH_USE_ORIGIN TRUE
    INSTALL_REMOVE_ENVIRONMENT_RPATH TRUE
    INSTALL_RPATH "${_origin}/${_relative}"
    INSTALL_RPATH_USE_LINK_PATH FALSE)
endfunction()
