# @author Oliver Thomson Brown
# @author Erich Essmann

configure_file(quest.h.in ../quest.h @ONLY)

file(GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" CONFIGURE_DEPENDS "*.h")
foreach(FILE IN LISTS HEADER_FILES)
	file(GENERATE OUTPUT ${FILE} INPUT ${FILE})
endforeach()
