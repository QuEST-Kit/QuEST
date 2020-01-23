# author: Ania Brown 
# author: Jacob Wilkins
# author: Balint Koczor (Windows compatibility)

# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(QuEST_Project)

# -----------------------------------------------------------------------------
# ----- USER OPTIONS ----------------------------------------------------------
# -----------------------------------------------------------------------------

set(USER_SOURCE  "examples/tutorial_example.c"  CACHE STRING "Source to build with QuEST library")
set(OUTPUT_EXE   "demo"  CACHE STRING "Executable to compile to")


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
if (WIN32)
    target_link_libraries(${OUTPUT_EXE} QuEST)
else ()
    target_link_libraries(${OUTPUT_EXE} QuEST m)
endif()
