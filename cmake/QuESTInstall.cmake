include(CMakePackageConfigHelpers)
set(quest_install_config_dir "${CMAKE_INSTALL_LIBDIR}/cmake/QuEST")
install(TARGETS QuEST EXPORT QuESTTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT Runtime NAMELINK_COMPONENT Development
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT Development
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT Runtime
  FILE_SET umbrella_header DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT Development
  FILE_SET api_headers DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT Development
  FILE_SET config_header DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT Development)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/QuESTConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}" COMPATIBILITY SameMajorVersion)
configure_package_config_file("${CMAKE_CURRENT_LIST_DIR}/QuESTConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfig.cmake"
  INSTALL_DESTINATION "${quest_install_config_dir}")
configure_file("${CMAKE_CURRENT_LIST_DIR}/QuESTConfigDependencies.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfigDependencies.cmake" @ONLY)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfigDependencies.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfigVersion.cmake"
  DESTINATION "${quest_install_config_dir}" COMPONENT Development)
install(FILES "${CMAKE_CURRENT_LIST_DIR}/FindNUMA.cmake"
  "${CMAKE_CURRENT_LIST_DIR}/FindCUQUANTUM.cmake" "${CMAKE_CURRENT_LIST_DIR}/FindCUTENSOR.cmake"
  DESTINATION "${quest_install_config_dir}/modules" COMPONENT Development)
install(EXPORT QuESTTargets FILE QuESTTargets.cmake NAMESPACE QuEST::
  DESTINATION "${quest_install_config_dir}" COMPONENT Development)
# Each independently installable package carries its license.
install(FILES "${PROJECT_SOURCE_DIR}/LICENCE.txt" "${PROJECT_SOURCE_DIR}/AUTHORS.txt"
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/licenses/QuEST/Development" COMPONENT Development)
if(QUEST_BUILT_SHARED)
  install(FILES "${PROJECT_SOURCE_DIR}/LICENCE.txt" "${PROJECT_SOURCE_DIR}/AUTHORS.txt"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/licenses/QuEST/Runtime" COMPONENT Runtime)
endif()
if(QUEST_HAVE_INSTALLABLE_EXAMPLES)
  install(FILES "${PROJECT_SOURCE_DIR}/LICENCE.txt" "${PROJECT_SOURCE_DIR}/AUTHORS.txt"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/licenses/QuEST/Examples" COMPONENT Examples)
endif()
