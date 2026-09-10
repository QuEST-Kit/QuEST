# @author Oliver Thomson Brown
# @author Erich Essmann
# @author Luc Jaulmes (using config file)
# @author Tyson Jones (doc)

# Generate a header file which defines all configurable preprocessors
# needed by the QuEST source (e.g. QUEST_COMPILE_MPI), as informed by the
# user-set CMake options. This permits us to avoid passing any macros
# through compiler flags and the associated conflicts arising when
# installing QuEST. Note that config.h must be manually created when
# not compiling via CMake, e.g. when using a custom build script 
configure_file(config.h.in config.h @ONLY)

# File sets preserve the historical public include layout.
set(_quest_api_headers
  calculations.h channels.h debug.h decoherence.h deprecated.h environment.h
  experimental.h initialisations.h matrices.h modes.h multiplication.h
  operations.h paulis.h precision.h qureg.h trotterisation.h types.h wrappers.h)
list(TRANSFORM _quest_api_headers PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
target_sources(QuEST PUBLIC
  FILE_SET umbrella_header TYPE HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" FILES quest.h
  FILE_SET api_headers TYPE HEADERS
    BASE_DIRS "${PROJECT_SOURCE_DIR}" FILES ${_quest_api_headers}
  FILE_SET config_header TYPE HEADERS
    BASE_DIRS "${PROJECT_BINARY_DIR}" FILES "${CMAKE_CURRENT_BINARY_DIR}/config.h")
