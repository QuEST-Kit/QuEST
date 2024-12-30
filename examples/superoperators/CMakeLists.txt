add_executable(superoperators_cpp_initialisation
  initialisation.cpp
)

target_link_libraries(superoperators_cpp_initialisation PUBLIC QuEST)

install(TARGETS superoperators_cpp_initialisation
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/superoperators
)
set_target_properties(superoperators_cpp_initialisation 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/quest/${CMAKE_INSTALL_LIBDIR}"
    OUTPUT_NAME "cpp_initialisation"
)

add_executable(superoperators_c_initialisation
  initialisation.cpp
)

target_link_libraries(superoperators_c_initialisation PUBLIC QuEST)

install(TARGETS superoperators_c_initialisation
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/superoperators
)
set_target_properties(superoperators_c_initialisation 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/quest/${CMAKE_INSTALL_LIBDIR}"
    OUTPUT_NAME c_initialisation
)