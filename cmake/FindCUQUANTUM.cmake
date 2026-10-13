#[=======================================================================[.rst:
FindCUQUANTUM
------------
Find the shared NVIDIA cuQuantum libraries (CMake 3.28 or newer).

@author Oliver Thomson Brown

Components are ``cuStateVec``, ``cuTensorNet`` and ``cuDensityMat``; a call
without components requests all three, required. Targets retain these names
under ``CUQUANTUM::``. ``CUQUANTUM::cuQuantum`` aggregates successfully found
components across calls. No CUDA language is needed.

Searches use normal CMake roots, environment roots, CMAKE_PREFIX_PATH and
cross-compilation rules. An explicitly supplied ``CUQUANTUM_DIR`` is a legacy
prefix hint searched first (cached artifact overrides still take precedence).
The environment ``CUQUANTUM_DIR`` is a fallback hint, never copied into that
variable. Only shared libraries are supported, not the SDK's _static archives.

Results: ``CUQUANTUM_FOUND``, ``CUQUANTUM_<component>_FOUND``,
``CUQUANTUM_<component>_INCLUDE_DIR``, ``CUQUANTUM_<component>_LIBRARY`` and
``CUQUANTUM_<component>_VERSION``. Versions describe individual components,
not the SDK release; package-level version requests are unsupported.
Legacy ``CUQUANTUM_INCLUDE_PATH``, ``CUQUANTUM_INCLUDE_DIRS``,
``CUQUANTUM_LIBRARIES`` and ``CUQUANTUM_LIBRARY_DIRS`` contain only successfully
resolved requested components; LIBRARIES contains imported target names.

cuTensorNet and cuDensityMat require cuTENSOR (``CUTENSOR_ROOT`` is supported).
See https://docs.nvidia.com/cuda/cuquantum/latest/getting-started/index.html.
#]=======================================================================]
include(FindPackageHandleStandardArgs)

function(_cuquantum_artifacts component stem)
  # Keep suffix restrictions local, including when a parent prefers static libs.
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .dll.a)
  elseif(APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .dylib .so)
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
  endif()
  if(CUQUANTUM_DIR)
    find_path(CUQUANTUM_${component}_INCLUDE_DIR NAMES ${stem}.h
      PATHS "${CUQUANTUM_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
    find_library(CUQUANTUM_${component}_LIBRARY NAMES ${stem}
      PATHS "${CUQUANTUM_DIR}" PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
  endif()
  find_path(CUQUANTUM_${component}_INCLUDE_DIR NAMES ${stem}.h
    HINTS ENV CUQUANTUM_DIR PATH_SUFFIXES include)
  find_library(CUQUANTUM_${component}_LIBRARY NAMES ${stem}
    HINTS ENV CUQUANTUM_DIR PATH_SUFFIXES lib lib64)
  mark_as_advanced(CUQUANTUM_${component}_INCLUDE_DIR CUQUANTUM_${component}_LIBRARY)
  set(version "")
  if(EXISTS "${CUQUANTUM_${component}_INCLUDE_DIR}/${stem}.h")
    string(TOUPPER "${stem}" macro)
    if(component STREQUAL "cuStateVec")
      string(APPEND macro "_VER")
    endif()
    file(STRINGS "${CUQUANTUM_${component}_INCLUDE_DIR}/${stem}.h" lines
      REGEX "^#[ \t]*define[ \t]+${macro}_(MAJOR|MINOR|PATCH)[ \t]+[0-9]+")
    foreach(field IN ITEMS MAJOR MINOR PATCH)
      set(${field} "")
      foreach(line IN LISTS lines)
        if(line MATCHES "${macro}_${field}[ \t]+([0-9]+)")
          set(${field} "${CMAKE_MATCH_1}")
        endif()
      endforeach()
    endforeach()
    if(NOT MAJOR STREQUAL "" AND NOT MINOR STREQUAL "" AND NOT PATCH STREQUAL "")
      set(version "${MAJOR}.${MINOR}.${PATCH}")
    endif()
  endif()
  set(CUQUANTUM_${component}_VERSION "${version}" PARENT_SCOPE)
endfunction()

set(_cuquantum_known cuStateVec cuTensorNet cuDensityMat)
if(NOT CUQUANTUM_FIND_COMPONENTS)
  set(CUQUANTUM_FIND_COMPONENTS ${_cuquantum_known})
  foreach(_component IN LISTS _cuquantum_known)
    set(CUQUANTUM_FIND_REQUIRED_${_component} TRUE)
  endforeach()
endif()
set(_cuquantum_needed ${CUQUANTUM_FIND_COMPONENTS})
if("cuDensityMat" IN_LIST _cuquantum_needed)
  list(APPEND _cuquantum_needed cuTensorNet)
endif()
# Dependencies must precede dependents, regardless of caller ordering.
set(_cuquantum_reasons "")
set(_cuquantum_has_known FALSE)
foreach(_component IN LISTS _cuquantum_needed)
  set(CUQUANTUM_${_component}_FOUND FALSE)
  if(_component IN_LIST _cuquantum_known)
    set(_cuquantum_has_known TRUE)
  else()
    list(APPEND _cuquantum_reasons "Unknown component '${_component}'")
  endif()
endforeach()
if(_cuquantum_has_known)
  find_package(CUDAToolkit QUIET)
endif()
if("cuTensorNet" IN_LIST _cuquantum_needed)
  find_package(CUTENSOR QUIET MODULE)
endif()
foreach(_component IN LISTS _cuquantum_known)
  if(NOT _component IN_LIST _cuquantum_needed)
    continue()
  endif()
  string(TOLOWER "${_component}" _stem)
  _cuquantum_artifacts("${_component}" "${_stem}")
  set(_deps CUDA::toolkit CUDA::cublas)
  if(_component STREQUAL "cuStateVec")
    list(APPEND _deps CUDA::cublasLt)
  elseif(_component STREQUAL "cuTensorNet")
    list(APPEND _deps CUDA::cusolver CUTENSOR::cutensor)
  else()
    list(APPEND _deps CUDA::cusolver CUDA::cublasLt CUDA::curand CUDA::cusparse
      CUTENSOR::cutensor CUQUANTUM::cuTensorNet)
    if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 12)
      list(APPEND _deps CUDA::nvJitLink)
    endif()
  endif()
  set(_ready TRUE)
  foreach(_dep IN LISTS _deps)
    if(NOT TARGET "${_dep}")
      set(_ready FALSE)
      list(APPEND _cuquantum_reasons "${_component} requires ${_dep}")
    endif()
  endforeach()
  if(NOT CUQUANTUM_${_component}_INCLUDE_DIR OR NOT CUQUANTUM_${_component}_LIBRARY
      OR NOT EXISTS "${CUQUANTUM_${_component}_INCLUDE_DIR}/${_stem}.h"
      OR NOT EXISTS "${CUQUANTUM_${_component}_LIBRARY}"
      OR CUQUANTUM_${_component}_LIBRARY MATCHES "(_static\\.|\\.a$)")
    set(_ready FALSE)
    list(APPEND _cuquantum_reasons "${_component} requires its header and shared library")
  endif()
  if(_component STREQUAL "cuDensityMat" AND NOT CUQUANTUM_cuTensorNet_FOUND)
    set(_ready FALSE)
  endif()
  set(CUQUANTUM_${_component}_FOUND ${_ready})
  if(_ready AND NOT TARGET CUQUANTUM::${_component})
    add_library(CUQUANTUM::${_component} UNKNOWN IMPORTED)
    set_target_properties(CUQUANTUM::${_component} PROPERTIES
      IMPORTED_LOCATION "${CUQUANTUM_${_component}_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${CUQUANTUM_${_component}_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${_deps}")
  endif()
endforeach()

set(CUQUANTUM_INCLUDE_DIRS "")
set(CUQUANTUM_LIBRARIES "")
set(CUQUANTUM_LIBRARY_DIRS "")
foreach(_component IN LISTS CUQUANTUM_FIND_COMPONENTS)
  if(CUQUANTUM_${_component}_FOUND)
    list(APPEND CUQUANTUM_INCLUDE_DIRS "${CUQUANTUM_${_component}_INCLUDE_DIR}")
    list(APPEND CUQUANTUM_LIBRARIES "CUQUANTUM::${_component}")
    get_filename_component(_libdir "${CUQUANTUM_${_component}_LIBRARY}" DIRECTORY)
    list(APPEND CUQUANTUM_LIBRARY_DIRS "${_libdir}")
  endif()
endforeach()
list(REMOVE_DUPLICATES CUQUANTUM_INCLUDE_DIRS)
list(REMOVE_DUPLICATES CUQUANTUM_LIBRARY_DIRS)
set(CUQUANTUM_INCLUDE_PATH "${CUQUANTUM_INCLUDE_DIRS}")
set(_cuquantum_version_supported TRUE)
if(CUQUANTUM_FIND_VERSION)
  set(_cuquantum_version_supported FALSE)
  list(APPEND _cuquantum_reasons "No overall SDK version is available: inspect CUQUANTUM_<component>_VERSION")
endif()
list(JOIN _cuquantum_reasons ". " _cuquantum_reason)
find_package_handle_standard_args(CUQUANTUM HANDLE_COMPONENTS
  REQUIRED_VARS _cuquantum_version_supported
  REASON_FAILURE_MESSAGE "${_cuquantum_reason}. Set CUQUANTUM_ROOT (or legacy CUQUANTUM_DIR), CUDAToolkit_ROOT and, for tensor components, CUTENSOR_ROOT.")
if(CUQUANTUM_FOUND)
  if(NOT TARGET CUQUANTUM::cuQuantum)
    add_library(CUQUANTUM::cuQuantum INTERFACE IMPORTED)
  endif()
  get_target_property(_aggregate CUQUANTUM::cuQuantum INTERFACE_LINK_LIBRARIES)
  if(NOT _aggregate)
    set(_aggregate "")
  endif()
  list(APPEND _aggregate ${CUQUANTUM_LIBRARIES})
  list(REMOVE_DUPLICATES _aggregate)
  set_target_properties(CUQUANTUM::cuQuantum PROPERTIES INTERFACE_LINK_LIBRARIES "${_aggregate}")
endif()
