# @author Oliver Thomson Brown
# @author Erich Essmann (updating rpath)

add_executable(reporters_cpp_env
  env.cpp
)

target_link_libraries(reporters_cpp_env PUBLIC QuEST)

install(TARGETS reporters_cpp_env
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/reporters
)
set_target_properties(reporters_cpp_env 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "cpp_env"
)

add_executable(reporters_c_env
  env.c
)

target_link_libraries(reporters_c_env PUBLIC QuEST)

install(TARGETS reporters_c_env
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/reporters
)
set_target_properties(reporters_c_env 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}}"
    OUTPUT_NAME "c_env"
)

add_executable(reporters_cpp_matrices
  matrices.cpp
)

target_link_libraries(reporters_cpp_matrices PUBLIC QuEST)

install(TARGETS reporters_cpp_matrices
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/reporters
)
set_target_properties(reporters_cpp_matrices 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "cpp_matrices"
)

add_executable(reporters_c_matrices
  matrices.c
)

target_link_libraries(reporters_c_matrices PUBLIC QuEST)

install(TARGETS reporters_c_matrices
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/reporters
)
set_target_properties(reporters_c_matrices
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "c_matrices"
)

add_executable(reporters_cpp_paulis
  paulis.cpp
)

target_link_libraries(reporters_cpp_paulis PUBLIC QuEST)

install(TARGETS reporters_cpp_paulis
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/reporters
)
set_target_properties(reporters_cpp_paulis
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "cpp_paulis"
)

add_executable(reporters_c_paulis
  paulis.c
)

target_link_libraries(reporters_c_paulis PUBLIC QuEST)

install(TARGETS reporters_c_paulis
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/reporters
)
set_target_properties(reporters_c_paulis 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "c_paulis"
)

add_executable(reporters_cpp_qureg
  qureg.cpp
)

target_link_libraries(reporters_cpp_qureg PUBLIC QuEST)

install(TARGETS reporters_cpp_qureg
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/reporters
)
set_target_properties(reporters_cpp_qureg
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "cpp_qureg"
)

add_executable(reporters_c_qureg
  qureg.c
)

target_link_libraries(reporters_c_qureg PUBLIC QuEST)

install(TARGETS reporters_c_qureg
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/reporters
)
set_target_properties(reporters_c_qureg
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
    OUTPUT_NAME "c_qureg"
)