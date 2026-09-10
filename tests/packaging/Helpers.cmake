function(run_checked)
  execute_process(COMMAND ${ARGV} RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Command failed (${result}): ${ARGV}\n${out}\n${err}")
  endif()
endfunction()

function(consume prefix binary_dir)
  file(MAKE_DIRECTORY "${binary_dir}")
  # An initial cache keeps list-valued prefixes and paths with spaces intact.
  file(WRITE "${binary_dir}/initial.cmake"
    "set(CMAKE_PREFIX_PATH [==[${prefix};${QUEST_DEPENDENCY_PREFIXES}]==] CACHE STRING \"\")\n")
  foreach(pair IN ITEMS "CMAKE_C_COMPILER|QUEST_C_COMPILER" "CMAKE_CXX_COMPILER|QUEST_CXX_COMPILER"
      "CMAKE_TOOLCHAIN_FILE|QUEST_TOOLCHAIN" "CUDAToolkit_ROOT|QUEST_CUDA_ROOT"
      "CUQUANTUM_ROOT|QUEST_CUQUANTUM_ROOT" "CUQUANTUM_DIR|QUEST_CUQUANTUM_DIR"
      "adios2_DIR|QUEST_ADIOS2_DIR" "NUMA_ROOT|QUEST_NUMA_ROOT"
      "HIP_DIR|QUEST_HIP_DIR" "MPI_CXX_COMPILER|QUEST_MPI_COMPILER")
    string(REPLACE "|" ";" fields "${pair}")
    list(GET fields 0 name)
    list(GET fields 1 source)
    if(NOT "${${source}}" STREQUAL "")
      file(APPEND "${binary_dir}/initial.cmake" "set(${name} [==[${${source}}]==] CACHE STRING \"\")\n")
    endif()
  endforeach()
  # Conflicting consumer feature options must not change the installed graph.
  run_checked("${CMAKE_COMMAND}" -S "${QUEST_SOURCE_DIR}/tests/packaging/consumer"
    -B "${binary_dir}" -C "${binary_dir}/initial.cmake"
    -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF
    -DQUEST_ENABLE_OMP=OFF -DQUEST_ENABLE_MPI=OFF -DBUILD_SHARED_LIBS=ON)
  run_checked("${CMAKE_COMMAND}" --build "${binary_dir}" --config "${CONFIG}" --parallel 2)
  run_checked("${CMAKE_CTEST_COMMAND}" --test-dir "${binary_dir}" -C "${CONFIG}" --output-on-failure)
endfunction()

function(check_relocation prefix)
  file(GLOB_RECURSE exports "${prefix}/*QuEST*.cmake")
  if(NOT exports)
    message(FATAL_ERROR "No installed QuEST CMake package")
  endif()
  foreach(export IN LISTS exports)
    file(READ "${export}" content)
    foreach(forbidden IN ITEMS "${QUEST_SOURCE_DIR}" "${QUEST_BUILD_DIR}")
      string(FIND "${content}" "${forbidden}" position)
      if(NOT position EQUAL -1)
        message(FATAL_ERROR "Producer path leaked into ${export}: ${forbidden}")
      endif()
    endforeach()
  endforeach()
  if(NOT EXISTS "${prefix}/${QUEST_INSTALL_INCLUDEDIR}/quest.h")
    message(FATAL_ERROR "Missing installed umbrella header")
  endif()
endfunction()
