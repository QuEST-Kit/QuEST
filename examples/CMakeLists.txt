add_executable(hadamard
  hadamard.cpp
)

target_link_libraries(hadamard PUBLIC QuEST)

install(TARGETS hadamard
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples
)
set_target_properties(hadamard PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/quest/${CMAKE_INSTALL_LIBDIR}")