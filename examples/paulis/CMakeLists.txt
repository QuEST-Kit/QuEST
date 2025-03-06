# @author Oliver Thomson Brown
# @author Erich Essmann (updating rpath)

add_executable(paulis_cpp_initialisation
  initialisation.cpp
)

target_link_libraries(paulis_cpp_initialisation PUBLIC QuEST)

install(TARGETS paulis_cpp_initialisation
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/paulis
)
set_target_properties(paulis_cpp_initialisation 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "cpp_initialisation"
)

add_executable(paulis_c_initialisation
  initialisation.c
)

target_link_libraries(paulis_c_initialisation PUBLIC QuEST)

install(TARGETS paulis_c_initialisation
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/paulis
)
set_target_properties(paulis_c_initialisation 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME c_initialisation
)