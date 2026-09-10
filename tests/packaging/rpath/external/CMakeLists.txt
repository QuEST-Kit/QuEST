cmake_minimum_required(VERSION 3.28)
project(QuESTRpathExternal LANGUAGES C)

add_library(rpath_external SHARED external.c)
set_target_properties(rpath_external PROPERTIES INSTALL_NAME_DIR "@rpath")
install(TARGETS rpath_external EXPORT RpathExternal LIBRARY DESTINATION lib)
install(EXPORT RpathExternal FILE RpathExternal.cmake NAMESPACE Rpath:: DESTINATION cmake)
