#[=======================================================================[.rst:
FindCuQuantum
-------------
 
Attempts to find NVIDIA's cuQuantum library. 
Use CUQUANTUM_ROOT or CUQUANTUM_DIR to specify the prefix path.
@author Oliver Thomson Brown


Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``CUQUANTUM_FOUND``
True if libcuquantum is found.
``CUQUANTUM_INCLUDE_DIRS``
Include directories needed to use cuQuantum.
``CUQUANTUM_LIBRARIES``
Libraries needed to link to cuQuantum.
``CUQUANTUM_LIBRARY_DIRS``
Location of libraries needed to link to cuQuantum.

#]=======================================================================]

include(FindPackageHandleStandardArgs)

find_package(PkgConfig QUIET)

# CUQUANTUM_DIR is the CMake standard, but cuQuantum uses CUQUANTUM_ROOT
# so we'll check if that's defined first
# A user supplied CUQUANTUM_DIR always takes precedence
if(NOT DEFINED CUQUANTUM_DIR)
  if(DEFINED ENV{CUQUANTUM_ROOT})
    set(CUQUANTUM_DIR $ENV{CUQUANTUM_ROOT})
  else()
    set(CUQUANTUM_DIR $ENV{CUQUANTUM_DIR})
  endif()
endif()

if (NOT CUQUANTUM_FOUND)
  # Until cuQuantum exports pkgconfig files or a CMake target
  # we're going to have to do this the hard way...

  # Look for custatevec.h in an include directory below CUQUANTUM_DIR
  # (or CUQUANTUM_ROOT)
  find_path(CUQUANTUM_INCLUDE_PATH
    NAMES 
      custatevec.h
      cudensitymat.h
      cutensornet.h
    PATHS
      ${CUQUANTUM_DIR}/include
  )

  set(CUQUANTUM_INCLUDE_DIRS "${CUQUANTUM_INCLUDE_PATH}")

  find_path(CUQUANTUM_LIBRARY_DIRS
    NAMES 
      libcustatevec.so
      libcudensitymat.so
      libcutensornet.so
    PATHS
      ${CUQUANTUM_DIR}/lib
      ${CUQUANTUM_DIR}/lib64
  )

  if(CUQUANTUM_LIBRARY_DIRS)
    set(CUQUANTUM_LIBRARIES "custatevec;cudensitymat;cutensornet")
  endif() 
endif()

find_package_handle_standard_args(CUQUANTUM
  REQUIRED_VARS
    CUQUANTUM_INCLUDE_DIRS
    CUQUANTUM_LIBRARIES
    CUQUANTUM_LIBRARY_DIRS
  REASON_FAILURE_MESSAGE
    "Try setting CUQUANTUM_DIR or CUQUANTUM_ROOT. Current values shown below.
      CUQUANTUM_DIR=${CUQUANTUM_DIR}
      CUQUANTUM_ROOT=${CUQUANTUM_ROOT}"
)

if(CUQUANTUM_FOUND AND NOT TARGET CUQUANTUM::cuQuantum)
  add_library(CUQUANTUM::cuQuantum INTERFACE IMPORTED)
  target_include_directories(CUQUANTUM::cuQuantum INTERFACE ${CUQUANTUM_INCLUDE_DIRS})
  target_link_libraries(CUQUANTUM::cuQuantum INTERFACE ${CUQUANTUM_LIBRARIES})
  target_link_directories(CUQUANTUM::cuQuantum INTERFACE ${CUQUANTUM_LIBRARY_DIRS})

  if(NOT TARGET CUQUANTUM::cuStateVec)
    add_library(CUQUANTUM::cuStateVec INTERFACE IMPORTED)
    target_include_directories(CUQUANTUM::cuStateVec INTERFACE ${CUQUANTUM_INCLUDE_DIRS})
    target_link_directories(CUQUANTUM::cuStateVec INTERFACE ${CUQUANTUM_LIBRARY_DIRS})
    target_link_libraries(CUQUANTUM::cuStateVec INTERFACE custatevec)
  endif()
    
  if(NOT TARGET CUQUANTUM::cuDensityMat)
    add_library(CUQUANTUM::cuDensityMat INTERFACE IMPORTED)
    target_include_directories(CUQUANTUM::cuDensityMat INTERFACE ${CUQUANTUM_INCLUDE_DIRS})
    target_link_directories(CUQUANTUM::cuDensityMat INTERFACE ${CUQUANTUM_LIBRARY_DIRS})
    target_link_libraries(CUQUANTUM::cuDensityMat INTERFACE cudensitymat)
  endif()
 
  if(NOT TARGET CUQUANTUM::cuTensorNet)
    add_library(CUQUANTUM::cuTensorNet INTERFACE IMPORTED)
    target_include_directories(CUQUANTUM::cuTensorNet INTERFACE ${CUQUANTUM_INCLUDE_DIRS})
    target_link_directories(CUQUANTUM::cuTensorNet INTERFACE ${CUQUANTUM_LIBRARY_DIRS})
    target_link_libraries(CUQUANTUM::cuTensorNet INTERFACE cutensornet)
  endif()
endif()
