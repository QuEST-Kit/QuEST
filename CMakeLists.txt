# author: Ania Brown 
# author: Jacob Wilkins
# author: Balint Koczor (Windows compatibility)
# author: Tyson Jones (testing)

# CMake initialisation.
cmake_minimum_required(VERSION 3.7)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

message(STATUS "CMake version: ${CMAKE_VERSION}")
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")

# Project name
project(QuEST_Project)

# -----------------------------------------------------------------------------
# ----- USER OPTIONS ----------------------------------------------------------
# -----------------------------------------------------------------------------

set(USER_SOURCE  "examples/tutorial_example.c"  CACHE STRING "Source to build with QuEST library")
set(OUTPUT_EXE   "demo"  CACHE STRING "Executable to compile to")
set(WITH_EXAMPLES "Whether to build the examples" 0)

option(TESTING "Enable unit testing. This disables compiling USER_SOURCE" OFF)

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

if (NOT TESTING)
    if (WITH_EXAMPLES)
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
        file( GLOB APP_SOURCES examples/*.c )
        message(STATUS "Compiling examples")
        foreach( examplesourcefile ${APP_SOURCES} )
            string( REPLACE ".c" "" sname ${examplesourcefile} )
            get_filename_component(ename ${sname} NAME )
            add_executable( ${ename} ${examplesourcefile} )
            if (WIN32)
                target_link_libraries( ${ename} QuEST )
            else ()
                target_link_libraries( ${ename} QuEST m)
            endif()
            install(TARGETS ${ename} DESTINATION examples )
        endforeach( examplesourcefile ${APP_SOURCES} )
    endif()
endif()

# -----------------------------------------------------------------------------
# ----- TESTS -----------------------------------------------------------------
# -----------------------------------------------------------------------------

if (TESTING)
    message(STATUS "Compiling unit tests")

    enable_testing()
    add_subdirectory(tests)
endif()
