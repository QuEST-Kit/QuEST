#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "QuEST::QuEST" for configuration "Release"
set_property(TARGET QuEST::QuEST APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(QuEST::QuEST PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libQuEST.so.4.1.0"
  IMPORTED_SONAME_RELEASE "libQuEST.so.4"
  )

list(APPEND _cmake_import_check_targets QuEST::QuEST )
list(APPEND _cmake_import_check_files_for_QuEST::QuEST "${_IMPORT_PREFIX}/lib/libQuEST.so.4.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
