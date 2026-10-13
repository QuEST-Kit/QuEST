#[=======================================================================[.rst:
FindCUTENSOR
------------
Find the shared cuTENSOR library and header using normal CMake search rules,
including CUTENSOR_ROOT and its environment equivalent. Defines CUTENSOR_FOUND,
CUTENSOR_INCLUDE_DIR, CUTENSOR_LIBRARY and CUTENSOR::cutensor. Does not enable
CUDA or accept the SDK's static archives.
#]=======================================================================]
include(FindPackageHandleStandardArgs)
function(_cutensor_find_artifacts)
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .dll.a)
  elseif(APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib .so)
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
  endif()
  find_path(CUTENSOR_INCLUDE_DIR NAMES cutensor.h PATH_SUFFIXES include)
  # NVIDIA archives separate CUDA major variants beneath lib/12 or lib/13.
  string(REGEX MATCH "^[0-9]+" _cuda_major "${CUDAToolkit_VERSION}")
  find_library(CUTENSOR_LIBRARY NAMES cutensor
    PATH_SUFFIXES "lib/${_cuda_major}" lib lib64)
endfunction()
find_package(CUDAToolkit QUIET)
_cutensor_find_artifacts()
set(_CUTENSOR_CUDA_FOUND FALSE)
if(TARGET CUDA::toolkit)
  set(_CUTENSOR_CUDA_FOUND TRUE)
endif()
set(_CUTENSOR_ARTIFACTS_VALID FALSE)
if(EXISTS "${CUTENSOR_INCLUDE_DIR}/cutensor.h" AND EXISTS "${CUTENSOR_LIBRARY}"
    AND NOT CUTENSOR_LIBRARY MATCHES "(_static\\.|\\.a$)")
  set(_CUTENSOR_ARTIFACTS_VALID TRUE)
endif()
find_package_handle_standard_args(CUTENSOR REQUIRED_VARS
  CUTENSOR_INCLUDE_DIR CUTENSOR_LIBRARY _CUTENSOR_CUDA_FOUND _CUTENSOR_ARTIFACTS_VALID)
mark_as_advanced(CUTENSOR_INCLUDE_DIR CUTENSOR_LIBRARY)
if(CUTENSOR_FOUND AND NOT TARGET CUTENSOR::cutensor)
  add_library(CUTENSOR::cutensor UNKNOWN IMPORTED)
  set_target_properties(CUTENSOR::cutensor PROPERTIES
    IMPORTED_LOCATION "${CUTENSOR_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CUTENSOR_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "CUDA::toolkit")
endif()
