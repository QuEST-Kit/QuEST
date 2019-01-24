# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(QuEST_Project)

# The location of the QuEST library sources
add_subdirectory(QuEST)

# Create user executable
add_executable(demo tutorial_example.c)

# Link libraries to user executable, including QuEST library
target_link_libraries(demo QuEST_lib m)
