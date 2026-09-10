cmake_minimum_required(VERSION 3.28)
project(QuESTRpathTests NONE)

enable_testing()
if(NOT UNIX)
  return()
endif()
get_filename_component(QUEST_MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake" ABSOLUTE)

foreach(case IN ITEMS default explicit_rpath use_link_path no_link_path
    cache_false target_false cache_true target_true custom_dirs skip_rpath
    skip_install_rpath build_with_install_rpath)
  add_test(NAME rpath.${case} COMMAND "${CMAKE_COMMAND}"
    --fresh -S "${CMAKE_CURRENT_SOURCE_DIR}/fixture" -B "${CMAKE_CURRENT_BINARY_DIR}/${case}"
    "-DQUEST_MODULE_DIR=${QUEST_MODULE_DIR}" "-DCASE=${case}")
endforeach()

if(APPLE)
  find_program(RPATH_INSPECTOR NAMES otool REQUIRED)
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
  find_program(RPATH_INSPECTOR NAMES readelf llvm-readelf REQUIRED)
else()
  return()
endif()

foreach(case IN ITEMS default explicit_rpath use_link_path no_link_path
    skip_rpath skip_install_rpath build_with_install_rpath)
  add_test(NAME rpath.runtime.${case} COMMAND "${CMAKE_COMMAND}"
    "-DQUEST_MODULE_DIR=${QUEST_MODULE_DIR}"
    "-DFIXTURE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/runtime ${case}"
    "-DRPATH_INSPECTOR=${RPATH_INSPECTOR}"
    "-DGENERATOR=${CMAKE_GENERATOR}" "-DCONFIG=$<CONFIG>"
    "-DCASE=${case}" -P "${CMAKE_CURRENT_SOURCE_DIR}/Runtime.cmake")
endforeach()
