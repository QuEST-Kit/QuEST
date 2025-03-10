# @author Oliver Thomson Brown

target_sources(QuEST
  PRIVATE
  gpu_config.cpp
  gpu_subroutines.cpp
)

if (ENABLE_CUDA)
  set_source_files_properties(
    gpu_config.cpp
    gpu_subroutines.cpp
    TARGET_DIRECTORY QuEST
    PROPERTIES
    LANGUAGE CUDA
  )
endif()

if (ENABLE_HIP)
  set_source_files_properties(
    gpu_config.cpp
    gpu_subroutines.cpp
    TARGET_DIRECTORY QuEST
    PROPERTIES
    LANGUAGE HIP
  )
endif()