cmake_minimum_required(VERSION 3.28)
project(QuESTConsumer LANGUAGES C CXX)

# Deliberately do not find any of QuEST's dependencies here.
find_package(QuEST 4 CONFIG REQUIRED)
find_package(QuEST 4 CONFIG REQUIRED)
if(CMAKE_CUDA_COMPILER_LOADED OR CMAKE_HIP_COMPILER_LOADED)
  message(FATAL_ERROR "QuEST package discovery must not enable GPU languages")
endif()
foreach(lang IN ITEMS c cpp)
  add_executable(consumer_${lang} main.${lang})
  target_link_libraries(consumer_${lang} PRIVATE QuEST::QuEST)
  set_target_properties(consumer_${lang} PROPERTIES
    C_STANDARD 11 C_STANDARD_REQUIRED YES
    CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)
endforeach()
enable_testing()
add_test(NAME consumer_c COMMAND consumer_c)
add_test(NAME consumer_cpp COMMAND consumer_cpp)

if(WIN32)
  set_tests_properties(consumer_c consumer_cpp PROPERTIES
    ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:$<TARGET_FILE_DIR:QuEST::QuEST>")
endif()
