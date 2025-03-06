# @author Oliver Thomson Brown
# @author Erich Essmann (updating rpath)

add_executable(cpp_arithmetic
  arithmetic.cpp
)

target_link_libraries(cpp_arithmetic PUBLIC QuEST)

install(TARGETS cpp_arithmetic
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/numbers
)
set_target_properties(cpp_arithmetic 
  PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
)

add_executable(c_arithmetic
  arithmetic.c
)

target_link_libraries(c_arithmetic PUBLIC QuEST)

install(TARGETS c_arithmetic
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}/examples/numbers
)
set_target_properties(c_arithmetic 
  PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
)