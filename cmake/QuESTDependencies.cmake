# Backend dependencies and target-local compilation requirements.

# OpenMP
if (QUEST_ENABLE_OMP)

  # find OpenMP, but fail gracefully...
  find_package(OpenMP QUIET COMPONENTS CXX)

  # so that we can customise the error message on MacOS
  if (NOT OpenMP_FOUND)
    set(ErrorMsg "Could not find OpenMP, necessary for enabling multithreading.")
    if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      string(APPEND ErrorMsg " Try first calling \n\tbrew install libomp\nthen\n\texport OpenMP_ROOT=$(brew --prefix)/opt/libomp")
    endif()
    message(FATAL_ERROR ${ErrorMsg})
  endif()

  target_link_libraries(QuEST
    PRIVATE
    OpenMP::OpenMP_CXX
  )

else()

  # suppress GCC "unknown pragma" warning when OpenMP disabled
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(QuEST PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wno-unknown-pragmas>)
  endif()

endif()


# NUMA is an optional enhancement, resolved once into the installed configuration.
if(QUEST_ENABLE_OMP AND QUEST_ENABLE_NUMA AND NOT WIN32)
  find_package(NUMA QUIET)
  if(NUMA_FOUND)
    target_link_libraries(QuEST PRIVATE NUMA::NUMA)
  else()
    message(WARNING "libnuma not found, QuEST will not be aware of NUMA locality")
    set(QUEST_ENABLE_NUMA OFF)
  endif()
else()
  set(QUEST_ENABLE_NUMA OFF)
endif()

# MPI
if (QUEST_ENABLE_MPI)
  find_package(MPI REQUIRED
    # Component CXX is the C api usable from C++
    # NOT the deprecated C++ API
    COMPONENTS CXX
  )

  if(QUEST_ENABLE_SUBCOMM)
    target_link_libraries(QuEST PUBLIC MPI::MPI_CXX)
  else()
    target_link_libraries(QuEST PRIVATE MPI::MPI_CXX)
  endif()
endif()


# CUDA
if (QUEST_ENABLE_CUDA)

  # make nvcc use user cxx-compiler as default host (before cuda-host is set below)
  if (NOT DEFINED CMAKE_CUDA_HOST_COMPILER)
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  endif()

  enable_language(CUDA)
  set_target_properties(QuEST PROPERTIES CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED YES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  find_package(CUDAToolkit REQUIRED)
  get_target_property(_quest_cuda_runtime QuEST CUDA_RUNTIME_LIBRARY)
  if(NOT _quest_cuda_runtime)
    if(CMAKE_CUDA_RUNTIME_LIBRARY_DEFAULT STREQUAL "SHARED")
      set(_quest_cuda_runtime Shared)
    else()
      set(_quest_cuda_runtime Static)
    endif()
    set_property(TARGET QuEST PROPERTY CUDA_RUNTIME_LIBRARY "${_quest_cuda_runtime}")
  endif()
  target_link_libraries(QuEST PRIVATE
    "$<$<STREQUAL:$<UPPER_CASE:${_quest_cuda_runtime}>,STATIC>:CUDA::cudart_static>"
    "$<$<STREQUAL:$<UPPER_CASE:${_quest_cuda_runtime}>,SHARED>:CUDA::cudart>")

  # force MSVC to use the modern preprocessor
  if (MSVC)
    target_compile_options(QuEST PRIVATE
      $<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>
      $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
    )
  endif()

endif()


# HIP
if (QUEST_ENABLE_HIP)

  # if generation fails (hip::amdhip64 not found), users can try setting
  # CMAKE_MODULE_PATH to '/opt/rocm/cmake' or '/opt/rocm/hip/lib/cmake/hip'
  # (suitable when shared library libamdhip64.so is located in /opt/rocm/lib/
  #  or /opt/rocm/hip/lib/ respectively). We avoid setting CMAKE_MODULE_PATH
  # pre-emptively since it made successful generation less likely in our tests!
  # example: list(APPEND CMAKE_MODULE_PATH "/opt/rocm/cmake"). Users should
  # also add '/opt/rocm/bin' or '/opt/rocm/hip/bin' to their $PATH env-var.

  enable_language(HIP)
  set_target_properties(QuEST PROPERTIES HIP_STANDARD 20 HIP_STANDARD_REQUIRED YES)

  find_package(HIP REQUIRED)
  message(STATUS "Found HIP: " ${HIP_VERSION})

  target_link_libraries(QuEST PRIVATE hip::host)

endif()


# cuQuantum
if (QUEST_ENABLE_CUQUANTUM)
  find_package(CUQUANTUM REQUIRED MODULE COMPONENTS cuStateVec)
  target_link_libraries(QuEST PRIVATE CUQUANTUM::cuStateVec)
endif()


# Checkpointing (ADIOS2)
if (QUEST_ENABLE_ADIOS2)

  set(_quest_adios2_components CXX)
  if(QUEST_ENABLE_MPI)
    list(APPEND _quest_adios2_components MPI)
  endif()
  find_package(adios2 CONFIG QUIET COMPONENTS ${_quest_adios2_components})
  if(QUEST_ENABLE_MPI)
    set(_quest_adios2_target adios2::cxx_mpi)
    set(_quest_adios2_legacy_target adios2::cxx11_mpi)
  else()
    set(_quest_adios2_target adios2::cxx)
    set(_quest_adios2_legacy_target adios2::cxx11)
  endif()
  # ADIOS2 2.9 (including Ubuntu 24.04) used the cxx11 target names.
  if(NOT TARGET ${_quest_adios2_target} AND TARGET ${_quest_adios2_legacy_target})
    set(_quest_adios2_target "${_quest_adios2_legacy_target}")
  endif()
  if(NOT adios2_FOUND AND (adios2_CONFIG OR TARGET adios2::core OR TARGET adios2::cxx OR TARGET adios2::cxx_mpi))
    message(FATAL_ERROR "The external ADIOS2 configuration was found but is unusable: ${adios2_NOT_FOUND_MESSAGE}. Select a compatible ADIOS2 installation; QuEST will not fetch over partially imported targets.")
  endif()
  if(adios2_FOUND AND NOT TARGET ${_quest_adios2_target})
    message(FATAL_ERROR "The installed ADIOS2 package does not provide ${_quest_adios2_target}. Select a compatible external ADIOS2 installation.")
  endif()
  if(NOT adios2_FOUND AND QUEST_ENABLE_INSTALL)
    message(FATAL_ERROR "Installable QuEST requires an external ADIOS2 package. Set adios2_DIR, or set QUEST_ENABLE_INSTALL=OFF and QUEST_ENABLE_PACKAGING=OFF for a developer build with downloaded ADIOS2.")
  endif()
  if(NOT adios2_FOUND AND QUEST_DOWNLOAD_ADIOS2)
    message(STATUS "fetching ADIOS2 via FetchContent")

    include(FetchContent)
    FetchContent_Declare(
      adios2
      GIT_REPOSITORY https://github.com/ornladios/ADIOS2.git
      GIT_TAG v2.12.1
    )

    # Match ADIOS2's MPI to QuEST's so distributed runs write per-rank slices
    # into one shared file. ADIOS2's CUDA support is deliberately left OFF:
    # checkpointing copies amps to host memory (syncQuregFromGpu/syncQuregToGpu)
    # before any I/O, so ADIOS2 never touches device pointers. Building it with
    # CUDA is unnecessary and stalls the Windows CUDA CI job.
    set(ADIOS2_USE_MPI  ${QUEST_ENABLE_MPI} CACHE BOOL "" FORCE)
    set(ADIOS2_USE_CUDA OFF CACHE BOOL "" FORCE)

    # Forego unused facilities
    set(ADIOS2_BUILD_TESTING  OFF CACHE BOOL "" FORCE)
    set(ADIOS2_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_SODIUM     OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_Fortran    OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_HDF5       OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_ZeroMQ     OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_SST        OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_DataMan    OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_SSC        OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_MHS        OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_DAOS       OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_MGARD      OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_BZip2      OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_Blosc      OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_Blosc2     OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_SZ         OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_ZFP        OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_PNG        OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_Profiling  OFF CACHE BOOL "" FORCE)
    set(ADIOS2_USE_Python     OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(adios2)

  else()
    # re-run non-QUIET so configuration fails with a clear error if the package
    # somehow became unavailable between the two calls
    find_package(adios2 CONFIG REQUIRED COMPONENTS ${_quest_adios2_components})
  endif()

  if(NOT TARGET ${_quest_adios2_target})
    message(FATAL_ERROR "ADIOS2 does not provide the required ${_quest_adios2_target} target")
  endif()

  # In distributed builds link ADIOS2's MPI-enabled C++ interface: it defines
  # ADIOS2_USE_MPI, which exposes the adios2::ADIOS(MPI_Comm) constructor used in
  # qureg.cpp for collective per-rank I/O. The serial target lacks it.
  target_link_libraries(QuEST PRIVATE ${_quest_adios2_target})
endif()
