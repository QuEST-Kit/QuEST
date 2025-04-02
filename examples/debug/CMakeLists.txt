# @author Tyson Jones (by duplicating files from...)
# @author Oliver Brown
# @author Erich Essmann


# C example

add_executable(debug_c_errorhandling
  errorhandling.c
)

target_link_libraries(debug_c_errorhandling PUBLIC QuEST)

install(TARGETS debug_c_errorhandling
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/debug
)
set_target_properties(debug_c_errorhandling 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "c_errorhandling"
)


# C++ example

add_executable(debug_cpp_errorhandling
  errorhandling.cpp
)

target_link_libraries(debug_cpp_errorhandling PUBLIC QuEST)

install(TARGETS debug_cpp_errorhandling
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/debug
)
set_target_properties(debug_cpp_errorhandling 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "cpp_errorhandling"
)
