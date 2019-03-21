# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(QuEST_Project)
message(${PROJECT_BINARY_DIR})
# -----------------------------------------------------------------------------
# ----- USER OPTIONS ----------------------------------------------------------
# -----------------------------------------------------------------------------

set(USER_SOURCE  "tutorial_example.c"  CACHE STRING "Source to build with QuEST library")
set(OUTPUT_EXE   "demo"  CACHE STRING "Executable to compile to")

option(TESTING "Enable test suite functionality -- requires test suite available" ON)
set(QuEST_TEST_DIR "tests" CACHE STRING
  "Name of the directory containing the QuEST library tests. It must be located in the same directory as the root CMakeLists.txt")

set(UPDATE_QUEST_LIB_LOC OFF CACHE BOOL "Set default QuESTPy library location to this compilation")

# -----------------------------------------------------------------------------
# ----- QuEST LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------

# Build the QuEST library if the path to libQuEST.so is not specified
if (NOT DEFINED ${QuEST_LIB_PATH})
    # Build libQuEST.so
    set(QuEST_DIR "QuEST" CACHE STRING 
        "Name of the directory containing the QuEST library sources. It must be located in the same directory as the root CMakeLists.txt")
    add_subdirectory(${QuEST_DIR})
    set(QuEST_LIB_PATH "${CMAKE_CURRENT_BINARY_DIR}/${QuEST_DIR}")
    set(QuEST_LIB_EXACT "${QuEST_LIB_PATH}/libQuEST.so")
endif()


# -----------------------------------------------------------------------------
# ----- USER EXECUTABLE -------------------------------------------------------
# -----------------------------------------------------------------------------

message(STATUS "Compiling ${USER_SOURCE} to executable ${OUTPUT_EXE}")

# Create user executable
add_executable(${OUTPUT_EXE} ${USER_SOURCE})

# Link libraries to user executable, including QuEST library
target_link_libraries(${OUTPUT_EXE} QuEST m)

# -----------------------------------------------------------------------------
# ----- TESTS -----------------------------------------------------------------
# -----------------------------------------------------------------------------

if (${TESTING})
  set(Python_ADDITIONAL_VERSIONS 3.4;3.5;3.6;3.7;3.8)
  find_package(PythonInterp)
  if ((NOT PYTHONINTERP_FOUND) OR (PYTHON_VERSION_MAJOR LESS 3)  OR (PYTHON_VERSION_MINOR LESS 4))
    set(TESTING OFF CACHE BOOL  "Enable test suite functionality -- requires test suite available" FORCE)
    message(WARNING "Python 3 not found on system -- Disabling testing")
  endif()
endif()

message(STATUS "Testing is ${TESTING}")

if (${TESTING})
  enable_testing()
  add_subdirectory(${QuEST_TEST_DIR})
endif()

# -----------------------------------------------------------------------------
# ----- QuEST Python Interface ------------------------------------------------
# -----------------------------------------------------------------------------

if (${UPDATE_QUEST_LIB_LOC})
  set(QUESTPY_LOC "${PROJECT_SOURCE_DIR}/${QuEST_TEST_DIR}/QuESTPy/QuESTLibDir.py")
  message(STATUS "Updating default QuEST Python library (${QUESTPY_LOC}) to ${QuEST_LIB_EXACT}")
  set(UPDATE_QUEST_LIB_LOC OFF CACHE BOOL "Set default QuESTPy library location to this compilation" FORCE)
  configure_file(${PROJECT_SOURCE_DIR}/cmake/QuESTPy.in ${QUESTPY_LOC})
endif()
