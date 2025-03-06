# @author Oliver Thomson Brown
# @author Erich Essmann (updating rpath)

add_executable(matrices_cpp_initialisation
  initialisation.cpp
)

target_link_libraries(matrices_cpp_initialisation PUBLIC QuEST)

install(TARGETS matrices_cpp_initialisation
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/matrices
)
set_target_properties(matrices_cpp_initialisation 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "cpp_initialisation"
)

add_executable(matrices_c_initialisation
  initialisation.c
)

target_link_libraries(matrices_c_initialisation PUBLIC QuEST)

install(TARGETS matrices_c_initialisation
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/matrices
)
set_target_properties(matrices_c_initialisation 
  PROPERTIES 
  INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
  OUTPUT_NAME "c_initialisation"
)