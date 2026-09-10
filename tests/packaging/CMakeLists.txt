if(NOT QUEST_ENABLE_INSTALL)
  message(FATAL_ERROR "QUEST_BUILD_PACKAGING_TESTS requires QUEST_ENABLE_INSTALL")
endif()
option(QUEST_TEST_ARCHIVES "Consume CPack binary archives in packaging tests" ${QUEST_ENABLE_PACKAGING})
configure_file(Settings.cmake.in Settings.cmake @ONLY)
add_test(NAME packaging.install COMMAND "${CMAKE_COMMAND}"
  "-DSETTINGS=${CMAKE_CURRENT_BINARY_DIR}/Settings.cmake"
  "-DCONFIG=$<CONFIG>" -P "${CMAKE_CURRENT_SOURCE_DIR}/InstallConsumer.cmake")
set_tests_properties(packaging.install PROPERTIES LABELS packaging RUN_SERIAL TRUE FIXTURES_SETUP quest_installed)
if(QUEST_TEST_ARCHIVES AND QUEST_ENABLE_PACKAGING)
  add_test(NAME packaging.archives COMMAND "${CMAKE_COMMAND}"
    "-DSETTINGS=${CMAKE_CURRENT_BINARY_DIR}/Settings.cmake"
    "-DCONFIG=$<CONFIG>" -P "${CMAKE_CURRENT_SOURCE_DIR}/ArchiveConsumer.cmake")
  set_tests_properties(packaging.archives PROPERTIES LABELS packaging RUN_SERIAL TRUE)
endif()
add_subdirectory(cuquantum)
get_property(_finder_tests DIRECTORY cuquantum PROPERTY TESTS)
set_tests_properties(${_finder_tests} DIRECTORY cuquantum PROPERTIES LABELS packaging)
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME packaging.cpack_policy COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpack/test_cpack.py")
  set_tests_properties(packaging.cpack_policy PROPERTIES LABELS packaging)
endif()

add_test(NAME packaging.behavior COMMAND "${CMAKE_COMMAND}"
  "-DSETTINGS=${CMAKE_CURRENT_BINARY_DIR}/Settings.cmake"
  "-DCONFIG=$<CONFIG>" -P "${CMAKE_CURRENT_SOURCE_DIR}/PackageBehavior.cmake")
set_tests_properties(packaging.behavior PROPERTIES LABELS packaging FIXTURES_REQUIRED quest_installed)

add_test(NAME packaging.configuration COMMAND "${CMAKE_COMMAND}"
  "-DSETTINGS=${CMAKE_CURRENT_BINARY_DIR}/Settings.cmake"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/Configuration.cmake")
set_tests_properties(packaging.configuration PROPERTIES LABELS packaging)

if(QUEST_ENABLE_CUQUANTUM AND NOT QUEST_BUILT_SHARED AND UNIX)
  add_test(NAME packaging.relocated_sdk COMMAND "${CMAKE_COMMAND}"
    "-DSETTINGS=${CMAKE_CURRENT_BINARY_DIR}/Settings.cmake"
    "-DCONFIG=$<CONFIG>"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/RelocatedSDK.cmake")
  set_tests_properties(packaging.relocated_sdk PROPERTIES LABELS packaging FIXTURES_REQUIRED quest_installed)
endif()
