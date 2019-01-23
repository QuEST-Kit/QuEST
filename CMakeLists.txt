# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(QuEST_Project)

add_subdirectory(QuEST)

add_executable(demo tutorial_example.c)
target_link_libraries(demo QuEST_lib m)
