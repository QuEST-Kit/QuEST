# @author Oliver Thomson Brown

target_sources(QuEST
  PRIVATE
  gpu_config.cpp
  gpu_subroutines.cpp
)

if (QUEST_ENABLE_CUDA)
  set_source_files_properties(
    gpu_config.cpp
    gpu_subroutines.cpp
    TARGET_DIRECTORY QuEST
    PROPERTIES
    LANGUAGE CUDA
  )
endif()

if (QUEST_ENABLE_HIP)
  set_source_files_properties(
    gpu_config.cpp
    gpu_subroutines.cpp
    TARGET_DIRECTORY QuEST
    PROPERTIES
    LANGUAGE HIP
  )
endif()