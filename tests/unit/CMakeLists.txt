# @author Oliver Thomson Brown

target_sources(tests
  PUBLIC
  bitwise.cpp
  calculations.cpp
  channels.cpp
  debug.cpp
  decoherence.cpp
  environment.cpp
  experimental.cpp
  initialisations.cpp
  matrices.cpp
  multiplication.cpp
  operations.cpp
  paulis.cpp
  qureg.cpp
  trotterisation.cpp
  types.cpp
)

# When the user opts in via QUEST_ENABLE_BMI2, compile only the issue-#717 bitwise test with -mbmi2 so
# it exercises the actual PEXT/PDEP path; otherwise (the default) it exercises the scalar fallback. The
# assertions hold identically either way. CMP0118: source properties are visible to the target's
# directory, so set it with explicit TARGET_DIRECTORY for the parent-scope 'tests' target.
if (QUEST_ENABLE_BMI2)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-mbmi2" QUEST_TEST_SUPPORTS_MBMI2)
  if (QUEST_TEST_SUPPORTS_MBMI2)
    set_source_files_properties(
      bitwise.cpp
      TARGET_DIRECTORY tests
      PROPERTIES COMPILE_OPTIONS "-mbmi2"
    )
  endif()
endif()