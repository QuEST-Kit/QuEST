target_sources(QuEST
  PRIVATE
  gpu_config.cpp
  gpu_subroutines.cpp
)

if (ENABLE_GPU_ACCELERATION)
  if (ENABLE_HIP)
    set_source_files_properties(
      gpu_config.cpp
      gpu_subroutines.cpp
      TARGET_DIRECTORY QuEST
      PROPERTIES
      LANGUAGE HIP
    )
  else()
    set_source_files_properties(
      gpu_config.cpp
      gpu_subroutines.cpp
      TARGET_DIRECTORY QuEST
      PROPERTIES
      LANGUAGE CUDA
    )
  endif()
endif()
