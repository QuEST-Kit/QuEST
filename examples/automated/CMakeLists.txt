# @author Tyson Jones

add_all_local_examples()

# The issue-#717 bitwise micro-benchmark builds with -mbmi2 (so its PEXT/PDEP path is enabled) only
# when the user opts in via QUEST_ENABLE_BMI2 — same switch the library uses. Without the opt-in it
# compiles the scalar fallback and prints "BMI2 fast path: INACTIVE" (never SIGILLs). add_example()
# names the target <filename>_<ext>; the flag is scoped to this one target.
if (QUEST_ENABLE_BMI2 AND TARGET benchmark_bitwise_bmi2_cpp)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-mbmi2" QUEST_EXAMPLE_SUPPORTS_MBMI2)
  if (QUEST_EXAMPLE_SUPPORTS_MBMI2)
    target_compile_options(benchmark_bitwise_bmi2_cpp PRIVATE -mbmi2)
  endif()
endif()
