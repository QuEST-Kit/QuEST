# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(QuEST_Project)

# The location of the QuEST library sources. Will build the QuEST_lib object if required
set(QuEST_DIR "QuEST")
add_subdirectory(${QuEST_DIR})


# Create user executable
add_executable(demo tutorial_example.c)

# Link libraries to user executable, including QuEST library
target_link_libraries(demo QuEST m)


set(QuEST_TEST_DIR "tests")
add_subdirectory(${QuEST_TEST_DIR})

