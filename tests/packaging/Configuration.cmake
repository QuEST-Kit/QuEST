include("${SETTINGS}")
include("${CMAKE_CURRENT_LIST_DIR}/Helpers.cmake")
set(work "${TEST_BINARY_DIR}/configuration")
file(REMOVE_RECURSE "${work}")
file(MAKE_DIRECTORY "${work}/parent")
file(WRITE "${work}/parent/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.28)
project(Parent LANGUAGES CXX)
set(CMAKE_BUILD_TYPE "" CACHE STRING "" FORCE)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS OFF)
set(QUEST_ENABLE_OMP OFF CACHE BOOL "")
add_custom_target(min_example)
add_custom_target(package)
add_subdirectory("${QUEST_SOURCE_DIR}" quest)
if(QUEST_ENABLE_INSTALL OR QUEST_ENABLE_PACKAGING OR QUEST_BUILD_MIN_EXAMPLE)
  message(FATAL_ERROR "Embedding QuEST unexpectedly enabled standalone facilities")
endif()
if(CMAKE_BUILD_TYPE OR CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS)
  message(FATAL_ERROR "QuEST changed parent build defaults")
endif()
if(NOT TARGET QuEST::QuEST)
  message(FATAL_ERROR "Missing subproject alias")
endif()
]=])
run_checked("${CMAKE_COMMAND}" -S "${work}/parent" -B "${work}/build"
  "-DQUEST_SOURCE_DIR=${QUEST_SOURCE_DIR}" "-DCMAKE_CXX_COMPILER=${QUEST_CXX_COMPILER}")

# A disabled installation must not generate package configs or install exports.
if(EXISTS "${work}/build/quest/QuESTConfig.cmake" OR EXISTS "${work}/build/quest/CPackConfig.cmake")
  message(FATAL_ERROR "A subproject emitted install/package configuration by default")
endif()

# ADIOS2 rejection must be deterministic even on machines with an installed SDK
# or no MPI installation. The MPI fixture supplies the discovery target only;
# these negative tests never compile QuEST or call an MPI function.
file(MAKE_DIRECTORY "${work}/adios2/modules" "${work}/adios2/fallback")
file(WRITE "${work}/adios2/modules/FindMPI.cmake" [=[
set(MPI_FOUND TRUE)
set(MPI_CXX_FOUND TRUE)
if(NOT TARGET MPI::MPI_CXX)
  add_library(MPI::MPI_CXX INTERFACE IMPORTED)
endif()
]=])
file(WRITE "${work}/adios2/fallback/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.28)
project(ForbiddenADIOS2Fallback LANGUAGES NONE)
file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/entered" "Fallback was entered")
message(FATAL_ERROR "FORBIDDEN_ADIOS2_FALLBACK: incompatible installed package must be rejected before fetching")
]=])

function(expect_adios2_rejection case expected)
  execute_process(COMMAND "${CMAKE_COMMAND}"
    -S "${QUEST_SOURCE_DIR}" -B "${work}/adios2/${case}/build"
    "-DCMAKE_C_COMPILER=${QUEST_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${QUEST_CXX_COMPILER}"
    "-DCMAKE_MODULE_PATH=${work}/adios2/modules"
    -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF
    -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF
    -DQUEST_BUILD_MIN_EXAMPLE=OFF -DQUEST_ENABLE_PACKAGING=OFF
    -DQUEST_ENABLE_OMP=OFF -DQUEST_ENABLE_ADIOS2=ON -DQUEST_DOWNLOAD_ADIOS2=ON
    "-DFETCHCONTENT_SOURCE_DIR_ADIOS2=${work}/adios2/fallback"
    ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
  set(output "${out}\n${err}")
  if(result EQUAL 0)
    message(FATAL_ERROR "ADIOS2 ${case}: incompatible configuration unexpectedly succeeded")
  endif()
  if(NOT output MATCHES "${expected}")
    message(FATAL_ERROR "ADIOS2 ${case}: missing expected rejection '${expected}':\n${output}")
  endif()
  if(output MATCHES "fetching ADIOS2|FORBIDDEN_ADIOS2_FALLBACK|already exists"
      OR EXISTS "${work}/adios2/fallback/entered")
    message(FATAL_ERROR "ADIOS2 ${case}: attempted fallback or conflicting imports:\n${output}")
  endif()
endfunction()

expect_adios2_rejection(missing_external "Installable QuEST requires an external ADIOS2 package"
  -DQUEST_ENABLE_INSTALL=ON -DCMAKE_DISABLE_FIND_PACKAGE_adios2=ON)

foreach(interface IN ITEMS serial mpi)
  set(config_dir "${work}/adios2/missing_${interface}/package")
  file(MAKE_DIRECTORY "${config_dir}")
  if(interface STREQUAL "serial")
    set(other_target adios2::cxx_mpi)
    set(expected_target adios2::cxx)
    set(enable_mpi OFF)
  else()
    set(other_target adios2::cxx)
    set(expected_target adios2::cxx_mpi)
    set(enable_mpi ON)
  endif()
  file(WRITE "${config_dir}/adios2-config.cmake"
    "add_library(${other_target} INTERFACE IMPORTED)\nset(adios2_FOUND TRUE)\n")
  expect_adios2_rejection("missing_${interface}" "does not provide ${expected_target}"
    -DQUEST_ENABLE_INSTALL=ON "-DQUEST_ENABLE_MPI=${enable_mpi}"
    "-Dadios2_DIR=${config_dir}")
endforeach()

set(config_dir "${work}/adios2/partially_imported/package")
file(MAKE_DIRECTORY "${config_dir}")
# Deliberately unguarded: a second find_package call would produce a duplicate
# target error, which must not replace QuEST's useful incompatibility diagnostic.
file(WRITE "${config_dir}/adios2-config.cmake" [=[
add_library(adios2::cxx INTERFACE IMPORTED)
set(adios2_FOUND FALSE)
set(adios2_NOT_FOUND_MESSAGE "fixture package is installed but incompatible")
]=])
expect_adios2_rejection(partially_imported "external ADIOS2 configuration was found but is unusable"
  -DQUEST_ENABLE_INSTALL=OFF "-Dadios2_DIR=${config_dir}")
