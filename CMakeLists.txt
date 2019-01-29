# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(QuEST_Project)

# -----------------------------------------------------------------------------
# ----- QuEST LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------

# Build the QuEST library if the path to libQuEST.so is not specified
if (NOT DEFINED ${QuEST_LIB_PATH})
    # Build libQuEST.so
    # The directory containing the QuEST library sources can be renamed but must be 
    # located in the same directory as this CMakeLists.txt
    if (NOT DEFINED ${QuEST_DIR})
        set(QuEST_DIR "QuEST")
    endif()
    add_subdirectory(${QuEST_DIR})
    set(QuEST_LIB_PATH "${CMAKE_CURRENT_BINARY_DIR}/${QuEST_DIR}")
endif()


# -----------------------------------------------------------------------------
# ----- USER EXECUTABLE -------------------------------------------------------
# -----------------------------------------------------------------------------

# Create user executable
add_executable(demo runTests.c)

# Link libraries to user executable, including QuEST library
target_link_libraries(demo QuEST m)


# -----------------------------------------------------------------------------
# ----- TESTS -----------------------------------------------------------------
# -----------------------------------------------------------------------------

# The directory containing the QuEST library tests can be renamed but must be 
# located in the same directory as this CMakeLists.txt
set(QuEST_TEST_DIR "tests")
add_subdirectory(${QuEST_TEST_DIR})

