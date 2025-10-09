# @author Oliver Thomson Brown
# @author Erich Essmann
# @author Luc Jaulmes (using config file)
# @author Tyson Jones (doc)

# Generate a header file which defines all configurable preprocessors
# needed by the QuEST source (e.g. COMPILE_MPI), as informed by the
# user-set CMake options. This permits us to avoid passing any macros
# through compiler flags and the associated conflicts arising when
# installing QuEST. Note that config.h must be manually created when
# not compiling via CMake, e.g. when using a custom build script 
configure_file(config.h.in config.h @ONLY)
