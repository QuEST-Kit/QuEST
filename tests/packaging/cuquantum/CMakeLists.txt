cmake_minimum_required(VERSION 3.28)
project(QuESTCuQuantumFinderTests NONE)
enable_testing()
get_filename_component(QUEST_MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake" ABSOLUTE)
foreach(case IN ITEMS partial missing_header missing_library static_only optional unknown_required unknown_optional default_all repeated versions legacy_root prefix_root env_root tensor_dependencies missing_cutensor missing_cuda recovery required_missing sdk_version)
  add_test(NAME cuquantum.${case} COMMAND "${CMAKE_COMMAND}"
    --fresh -S "${CMAKE_CURRENT_SOURCE_DIR}/fixture" -B "${CMAKE_CURRENT_BINARY_DIR}/${case}"
    "-DQUEST_MODULE_DIR=${QUEST_MODULE_DIR}" "-DCASE=${case}")
  if(case STREQUAL "required_missing")
    set_tests_properties(cuquantum.${case} PROPERTIES WILL_FAIL TRUE)
  endif()
endforeach()
