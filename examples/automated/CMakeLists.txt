# @author Tyson Jones

add_all_local_examples()

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mbmi2" QUEST_COMPILER_SUPPORTS_MBMI2)
if (QUEST_COMPILER_SUPPORTS_MBMI2)
  target_compile_options(benchmark_bmi2_bitwise_cpp PRIVATE -mbmi2)
endif()
