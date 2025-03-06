# @author Oliver Thomson Brown
# @author Erich Essmann (updating rpath)

add_executable(krausmap_cpp_initialisation
  initialisation.cpp
)

target_link_libraries(krausmap_cpp_initialisation PUBLIC QuEST)

install(TARGETS krausmap_cpp_initialisation
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/krausmaps
)
set_target_properties(krausmap_cpp_initialisation 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "cpp_initialisation"
)

add_executable(krausmap_c_initialisation
  initialisation.c
)

target_link_libraries(krausmap_c_initialisation PUBLIC QuEST)

install(TARGETS krausmap_c_initialisation
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/krausmaps
)
set_target_properties(krausmap_c_initialisation 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "c_initialisation"
)