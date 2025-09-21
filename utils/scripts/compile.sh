# A hacky Unix bash script which compiles QuEST alongside user
# source code, or the v4 tests, or the deprecated v3 unit tests, 
# into a single executable. This script should be copied to, and 
# run in, the root directory (i.e. that containing the quest/ direc). 
# Beware, this is NOT the recommended method of compiling QuEST, and
# users should instead use the CMake build. This script is mostly
# useful for experimenting with the flags passed to the unix build.
# 
# @author Tyson Jones



# USER SETTINGS

# numerical precision (1, 2, 4)
FLOAT_PRECISION=2

# deployments to compile (0, 1)
ENABLE_DISTRIBUTION=0       # MPI
ENABLE_MULTITHREADING=0     # OpenMP
ENABLE_CUDA=0               # NVIDIA GPU
ENABLE_HIP=0                # AMD GPU
ENABLE_CUQUANTUM=0          # NVIDIA cuStateVec
ENABLE_NUMA=0               # NUMA awareness

# other options (0, 1)
ENABLE_DEPRECATED_API=0
DISABLE_DEPRECATION_WARNINGS=0

# NVIDIA compute capability or AMD arch (e.g. 60 or gfx908)
GPU_ARCH=90

# backend compilers
TESTS_COMPILER=g++
BASE_COMPILER=g++
OMP_COMPILER=g++
MPI_COMPILER=mpic++
CUDA_COMPILER=nvcc
HIP_COMPILER=hipcc

# linker
LINKER=g++

# whether to compile the below user source files (0),
# or the unit tests (1), which when paired with above
# ENABLE_DEPRECATED_API=1, will use the v3 tests (which
# you should pair with DISABLE_DEPRECATION_WARNINGS=1)
COMPILE_TESTS=0

# name of the compiled test executable
TEST_EXEC_FILE="test"

# name of the compiled user-code executable
USER_EXEC_FILE="main"

# user files (.c or .cpp, intermixed)
USER_FILES=(
    "main.cpp"
)

# location of user files
USER_DIR="."

# user code language-specific compilers
USER_C_COMPILER=gcc
USER_CXX_COMPILER=g++

# user code language-specific compiler flags
USER_C_COMP_FLAGS='-std=c11'
USER_CXX_COMP_FLAGS='-std=c++14'

# user linker flags
USER_LINK_FLAGS='-lstdc++'

# whether to compile cuQuantum (consulted only when ENABLE_CUQUANTUM=1)
# in debug mode, which logs to below file with performance tips and errors
CUQUANTUM_LOG=0
CUQUANTUM_LOG_FN="./custatevec_log.txt"

# external library locations (replace with "." to default)
CUQUANTUM_LIB_DIR="${CUQUANTUM_ROOT}"
CUDA_LIB_DIR="/usr/local/cuda"
ROCM_LIB_DIR="/opt/rocm"
OMP_LIB_DIR="/opt/homebrew/opt/libomp"
MPI_LIB_DIR="/opt/homebrew/opt/openmpi"
CATCH_LIB_DIR="$(pwd)/catch"



# CONSTANTS

USER_OBJ_PREF="user_"
QUEST_OBJ_PREF="quest_"
TEST_OBJ_PREF='test_'

INDENT='  '

CATCH_VERSION="3.4.0"



# QUEST FILE LAYOUT

# directories
ROOT_DIR='quest'
INCLUDE_DIR="${ROOT_DIR}/include"

SRC_DIR="${ROOT_DIR}/src"
API_DIR="${SRC_DIR}/api"
OMP_DIR="${SRC_DIR}/cpu"
GPU_DIR="${SRC_DIR}/gpu"
MPI_DIR="${SRC_DIR}/comm"
CORE_DIR="${SRC_DIR}/core"

TEST_MAIN_DIR="tests"
TEST_UTIL_DIR="${TEST_MAIN_DIR}/utils"
TEST_UNIT_DIR="${TEST_MAIN_DIR}/unit"
TEST_DEPR_DIR="${TEST_MAIN_DIR}/deprecated"
TEST_DEPR_CATCH_DIR="${TEST_DEPR_DIR}/catch"

# files that require modification by this script
CONFIG_FILE_IN="${INCLUDE_DIR}/config.h.in"
CONFIG_FILE_OUT="${INCLUDE_DIR}/config.h"

# files in API_DIR
API_FILES=(
    "calculations"
    "channels"
    "debug"
    "decoherence"
    "environment"
    "initialisations"
    "matrices"
    "modes"
    "multiplication"
    "operations"
    "paulis"
    "qureg"
    "trotterisation"
    "types"
)

# files in CORE_DIR
CORE_FILES=(
    "accelerator"
    "autodeployer"
    "envvars"
    "errors"
    "localiser"
    "memory"
    "parser"
    "paulilogic"
    "printer"
    "randomiser"
    "utilities"
    "validation"
)

# files in GPU_DIR
GPU_FILES=(
    "gpu_config"
    "gpu_subroutines"
)

# files in OMP_DIR
OMP_FILES=(
    "cpu_config"
    "cpu_subroutines"
)

# files in MPI_DIR
MPI_FILES=(
    "comm_config"
    "comm_routines"
)

# files in TEST_MAIN_DIR
TEST_MAIN_FILES=(
    "main"
)

# files in TEST_UTIL_DIR
TEST_UTIL_FILES=(
    "cache"
    "compare"
    "config"
    "convert"
    "evolve"
    "linalg"
    "lists"
    "measure"
    "qmatrix"
    "qvector"
    "random"
)

# files in TEST_UNIT_DIR
TEST_UNIT_FILES=(
    "calculations"
    "channels"
    "debug"
    "decoherence"
    "environment"
    "initialisations"
    "matrices"
    "multiplication"
    "operations"
    "paulis"
    "qureg"
    "trotterisation"
    "types"
)

# files in TEST_DEPR_DIR
TEST_DEPR_FILES=(
    "test_main"
    "test_gates"
    "test_unitaries"
    "test_operators"
    "test_decoherence"
    "test_calculations"
    "test_data_structures"
    "test_state_initialisations"
)

# files in TEST_DEPR_DIR which use MPI
TEST_DEPR_MPI_FILES=(
    "test_utilities"
)



# COMPILER AND LINKER FLAG OPTIONS

# compiler flags given to all (non-deprecated) files
TEST_COMP_FLAGS="-std=c++20 -I${CATCH_LIB_DIR}/include"
TEST_LINK_FLAGS="-L${CATCH_LIB_DIR}/lib -lCatch2"

# compiler flags given to deprecated test files
TEST_DEPR_COMP_FLAGS="-std=c++17 -I${CATCH_LIB_DIR}/include"

# compiler flags given to all backend files
BACKEND_COMP_FLAGS='-std=c++17 -O3'

# warning flags which apply to all compiled and linked files including user's
WARNING_FLAGS='-Wall'

# CUDA specific flags
CUDA_COMP_FLAGS="-x cu -arch=sm_${GPU_ARCH} -I${CUDA_LIB_DIR}/include"
CUDA_LINK_FLAGS="-L${CUDA_LIB_DIR}/lib -L${CUDA_LIB_DIR}/lib64 -lcudart -lcuda"

if [ $ENABLE_CUQUANTUM == 1 ]
then
    # extend GPU flags if cuQuantum enabled
    CUDA_COMP_FLAGS+=" -I${CUQUANTUM_LIB_DIR}/include"
    CUDA_LINK_FLAGS+=" -L${CUQUANTUM_LIB_DIR}/lib -L${CUQUANTUM_LIB_DIR}/lib64 -lcustatevec"

    # optional debug logging - will slow down code
    if [ $CUQUANTUM_LOG == 1 ]
    then
        CUDA_COMP_FLAGS+=" -DCUSTATEVEC_LOG_LEVEL=5 -DCUSTATEVEC_LOG_FILE=${CUQUANTUM_LOG_FN}"
    fi
fi

# HIP specific flags
HIP_COMP_FLAGS="-x hip --offload-arch=${GPU_ARCH} -I${ROCM_LIB_DIR}/include"
HIP_LINK_FLAGS="-L${ROCM_LIB_DIR}/lib -lamdhip64"

# MPI-specific flags
MPI_COMP_FLAGS="-I${MPI_LIB_DIR}/include"
MPI_LINK_FLAGS="-L${MPI_LIB_DIR}/lib -lmpi"

if ! $LINKER --version | grep -iq "clang"
then
    MPI_LINK_FLAGS+=' -lmpi_cxx'
fi

# OpenMP specific flags (compiler dependent)
OMP_COMP_FLAGS="-I${OMP_LIB_DIR}/include"
OMP_LINK_FLAGS="-L${OMP_LIB_DIR}/lib"

if $OMP_COMPILER --version | grep -iq "clang"
then
    OMP_COMP_FLAGS+=' -Xclang -fopenmp'
else
    OMP_COMP_FLAGS+=' -fopenmp'
fi

if $LINKER --version | grep -iq "clang"
then
    OMP_LINK_FLAGS+=' -lomp'
else
    OMP_LINK_FLAGS+=' -fopenmp'
fi

if [ $ENABLE_NUMA == 1 ]
then
    OMP_LINK_FLAGS+=' -lnuma'
fi

# point compilers to QuEST src
HEADER_FLAGS="-I. -I${INCLUDE_DIR}"



# ASSEMBLE FLAGS

echo ""
echo "deployment modes:"

# choose linker flags (extended below)
ALL_LINK_FLAGS="${USER_LINK_FLAGS}"

# choose compiler and flags for CPU/OMP files
CPU_FILES_FLAGS='-Ofast -DCOMPLEX_OVERLOADS_PATCHED=1'

if [ $ENABLE_MULTITHREADING == 1 ]
then
    echo "${INDENT}(multithreading enabled)"
    echo "${INDENT}${INDENT}[compiling OpenMP]"
    if [ $ENABLE_NUMA == 1 ]
    then
        echo "${INDENT}${INDENT}[compiling NUMA]"
    fi
    CPU_FILES_COMPILER=$OMP_COMPILER
    CPU_FILES_FLAGS+=" ${OMP_COMP_FLAGS}"
    ALL_LINK_FLAGS+=" ${OMP_LINK_FLAGS}"
else
    echo "${INDENT}(multithreading disabled)"
    CPU_FILES_COMPILER=$BASE_COMPILER
fi

# choose compiler and flags for GPU files
if [ $ENABLE_CUDA == 1 ]
then
    echo "${INDENT}(GPU-acceleration enabled)"
    echo "${INDENT}${INDENT}[compiling CUDA]"
    GPU_FILES_COMPILER=$CUDA_COMPILER
    GPU_FILES_FLAGS=$CUDA_COMP_FLAGS
    ALL_LINK_FLAGS+=" ${CUDA_LINK_FLAGS}"
    GPU_WARNING_FLAGS="-Xcompiler ${WARNING_FLAGS}"
elif [ $ENABLE_HIP == 1 ]
then
    echo "${INDENT}(GPU-acceleration enabled)"
    echo "${INDENT}${INDENT}[compiling HIP]"
    GPU_FILES_COMPILER=$HIP_COMPILER
    GPU_FILES_FLAGS=$HIP_COMP_FLAGS
    ALL_LINK_FLAGS+=" ${HIP_LINK_FLAGS}"
    GPU_WARNING_FLAGS="-Xcompiler ${WARNING_FLAGS}"
else
    echo "${INDENT}(GPU-acceleration disabled)"
    GPU_FILES_COMPILER=$BASE_COMPILER
    GPU_FILES_FLAGS=''
    GPU_WARNING_FLAGS=$WARNING_FLAGS
fi

# merely report cuQuantum status
if [ $ENABLE_CUQUANTUM == 1 ]
then
    echo "${INDENT}(cuQuantum enabled)"
    echo "${INDENT}${INDENT}[compiling cuStateVec]"
else
    echo "${INDENT}(cuQuantum disabled)"
fi

# choose compiler and flags for communication files
if [ $ENABLE_DISTRIBUTION == 1 ]
then
    echo "${INDENT}(distribution enabled)"
    echo "${INDENT}${INDENT}[compiling MPI]"
    MPI_FILES_COMPILER=$MPI_COMPILER
    MPI_FILES_FLAGS=$MPI_COMP_FLAGS
    ALL_LINK_FLAGS+=" ${MPI_LINK_FLAGS}"
else
    echo "${INDENT}(distribution disabled)"
    MPI_FILES_COMPILER=$BASE_COMPILER
    MPI_FILES_FLAGS=''
fi

# choose linker warning flag (to avoid passing them to nvcc)
if [ "${LINKER}" = "nvcc" ]
then
    ALL_LINK_FLAGS+="-Xcompiler ${WARNING_FLAGS}"
else
    ALL_LINK_FLAGS+=" ${WARNING_FLAGS}"
fi

# test link flags
if [ "${COMPILE_TESTS}" -eq 1 ] 
then
    ALL_LINK_FLAGS+=" ${TEST_LINK_FLAGS}"
fi

# display precision
if [ $FLOAT_PRECISION == 1 ]; then
    echo "${INDENT}(single precision)"
elif [ $FLOAT_PRECISION == 2 ]; then
    echo "${INDENT}(double precision)"
elif [ $FLOAT_PRECISION == 4 ]; then
    echo "${INDENT}(quad precision)"
else
    echo ""
    echo "INVALID FLOAT_PRECISION (${FLOAT_PRECISION})"
    echo "Exiting..."
    exit
fi

echo ""



# REPORTING COMILERS FLAGS

echo "chosen compilers and flags..."

# user compilers
if (( $COMPILE_TESTS == 0 ))
then
    echo "${INDENT}user compilers and flags (for .c and .cpp files respectively):"
    echo "${INDENT}${INDENT}${USER_C_COMPILER} ${USER_C_COMP_FLAGS} ${WARNING_FLAGS}"
    echo "${INDENT}${INDENT}${USER_CXX_COMPILER} ${USER_CXX_COMP_FLAGS} ${WARNING_FLAGS}"
fi

# test compiler
if (( $COMPILE_TESTS == 1 && ENABLE_DEPRECATED_API == 0 ))
then
    echo "${INDENT}tests compiler and flags:"
    echo "${INDENT}${INDENT}${TESTS_COMPILER} ${TEST_COMP_FLAGS} ${WARNING_FLAGS}"
fi

# deprecated compiler
if (( $COMPILE_TESTS == 1 && ENABLE_DEPRECATED_API == 1 ))
then
    echo "${INDENT}deprecated tests compiler and flags:"
    echo "${INDENT}${INDENT}${TESTS_COMPILER} ${TEST_DEPR_COMP_FLAGS} ${WARNING_FLAGS}"
fi

# base compiler
echo "${INDENT}base files compiler and backend flags:"
echo "${INDENT}${INDENT}${BASE_COMPILER} ${BACKEND_COMP_FLAGS} ${WARNING_FLAGS}"

# OMP/CPU
echo "${INDENT}CPU files compiler and flags:"
echo "${INDENT}${INDENT}${CPU_FILES_COMPILER} ${CPU_FILES_FLAGS} ${WARNING_FLAGS}"

# GPU
echo "${INDENT}GPU files compiler and flags:"
echo "${INDENT}${INDENT}${GPU_FILES_COMPILER} ${GPU_FILES_FLAGS} ${GPU_WARNING_FLAGS}"

# MPI
echo "${INDENT}distributed files compiler and flags:"
echo "${INDENT}${INDENT}${MPI_FILES_COMPILER} ${MPI_FILES_FLAGS} ${WARNING_FLAGS}"

# linker 
echo "${INDENT}linker and all linker flags:"
echo "${INDENT}${INDENT}${LINKER} ${ALL_LINK_FLAGS}"

# globals
echo "${INDENT}header flags:"
echo "${INDENT}${INDENT}${HEADER_FLAGS}"
echo ""



# OPTIONALLY PREPARING CATCH2

if [ "${COMPILE_TESTS}" -eq 1 ] 
then
    echo "preparing Catch2:"

    if [ -d "${CATCH_LIB_DIR}" ]
    then
        echo "${INDENT}found at ${CATCH_LIB_DIR}"
    else
        echo "${INDENT}downloading to ${CATCH_LIB_DIR}..."
        git clone --quiet https://github.com/catchorg/Catch2.git "${CATCH_LIB_DIR}"

        ORIGINAL_DIR=$(pwd)

        echo "${INDENT}configuring..."
        cd "${CATCH_LIB_DIR}"
        git fetch --quiet --tags
        git checkout --quiet "v${CATCH_VERSION}"
        git submodule update --quiet --init --recursive

        echo "${INDENT}building..."
        mkdir build
        cd build
        cmake .. -DCMAKE_INSTALL_PREFIX="${CATCH_LIB_DIR}"
        cmake --build . --target install --parallel

        cd "${ORIGINAL_DIR}"
    fi

    echo ""
fi



# GENERATING CONFIG HEADER

echo "generating headers:"

# write user-options as macros to config.h (and set version info to -1)
sed \
  -e "s|#cmakedefine FLOAT_PRECISION @FLOAT_PRECISION@|#define FLOAT_PRECISION ${FLOAT_PRECISION}|" \
  -e "s|#cmakedefine01 INCLUDE_DEPRECATED_FUNCTIONS|#define INCLUDE_DEPRECATED_FUNCTIONS ${ENABLE_DEPRECATED_API}|" \
  -e "s|#cmakedefine01 DISABLE_DEPRECATION_WARNINGS|#define DISABLE_DEPRECATION_WARNINGS ${DISABLE_DEPRECATION_WARNINGS}|" \
  -e "s|#cmakedefine01 COMPILE_OPENMP|#define COMPILE_OPENMP ${ENABLE_MULTITHREADING}|" \
  -e "s|#cmakedefine01 COMPILE_MPI|#define COMPILE_MPI ${ENABLE_DISTRIBUTION}|" \
  -e "s|#cmakedefine01 COMPILE_CUDA|#define COMPILE_CUDA $(( ENABLE_CUDA || ENABLE_HIP ))|" \
  -e "s|#cmakedefine01 COMPILE_CUQUANTUM|#define COMPILE_CUQUANTUM ${ENABLE_CUQUANTUM}|" \
  -e "s|#cmakedefine01 COMPILE_HIP|#define COMPILE_HIP ${ENABLE_HIP}|" \
  -e "s|#cmakedefine01 NUMA_AWARE|#define NUMA_AWARE ${ENABLE_NUMA}|" \
  -e "s|@PROJECT_VERSION@|unknown (not populated by manual compilation)|" \
  -e "s|@PROJECT_VERSION_MAJOR@|-1|" \
  -e "s|@PROJECT_VERSION_MINOR@|-1|" \
  -e "s|@PROJECT_VERSION_PATCH@|-1|" \
  "${CONFIG_FILE_IN}" > "${CONFIG_FILE_OUT}"

echo "${INDENT}${CONFIG_FILE_OUT}"
echo ""



# COMPILING USER CODE

# abort script if any compilation fails
set -e

if (( $COMPILE_TESTS == 0 ))
then
    echo "compiling user files:"

    for fn in ${USER_FILES[@]}
    do
        # choose C or C++ compiler for each user file
        if [[ $fn == *.cpp ]]
        then
            echo "${INDENT}${fn} (C++) ..."
            COMP=$USER_CXX_COMPILER
            FLAG=$USER_CXX_COMP_FLAGS
        else
            echo "${INDENT}${fn} (C) ..."
            COMP=$USER_C_COMPILER
            FLAG=$USER_C_COMP_FLAGS
        fi

        # compile
        $COMP -c $USER_DIR/$fn -o ${USER_OBJ_PREF}${fn}.o $FLAG $HEADER_FLAGS $WARNING_FLAGS
    done

    echo ""
fi



# COMPILING TESTS

if (( $COMPILE_TESTS == 1 && $ENABLE_DEPRECATED_API == 0 ))
then

    echo "compiling unit test files:"

    echo "${INDENT}main:"

    for fn in ${TEST_MAIN_FILES[@]}
    do
        echo "${INDENT}${INDENT}${fn}.cpp ..."
        $TESTS_COMPILER -c $TEST_MAIN_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
    done

    echo "${INDENT}utils:"

    for fn in ${TEST_UTIL_FILES[@]}
    do
        echo "${INDENT}${INDENT}${fn}.cpp ..."
        $TESTS_COMPILER -c $TEST_UTIL_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
    done

    echo "${INDENT}unit:"

    for fn in ${TEST_UNIT_FILES[@]}
    do
        echo "${INDENT}${INDENT}${fn}.cpp ..."
        $TESTS_COMPILER -c $TEST_UNIT_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
    done

    echo ""
fi



# COMPILING DEPRECATED TESTS

if (( $COMPILE_TESTS == 1 && $ENABLE_DEPRECATED_API == 1 ))
then
    echo "compiling deprecated test files:"

    if (( $DISABLE_DEPRECATION_WARNINGS == 0 ))
    then
        echo "${INDENT}(beware deprecation warnings were not disabled)"
    fi

    for fn in ${TEST_DEPR_FILES[@]}
    do
        echo "${INDENT}${fn}.cpp ..."
        $TESTS_COMPILER -c $TEST_DEPR_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_DEPR_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
    done

    for fn in ${TEST_DEPR_MPI_FILES[@]}
    do
        echo "${INDENT}${fn}.cpp ..."
        $MPI_FILES_COMPILER -c $TEST_DEPR_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_DEPR_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
    done

    echo ""
fi



# COMPILING CORE

echo "compiling core files in C++"

for fn in ${CORE_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $BASE_COMPILER -c $CORE_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $BACKEND_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
done

echo ""



# COMPILING API

echo "compiling API files in C++:"

for fn in ${API_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $BASE_COMPILER -c $API_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $BACKEND_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
done

echo ""



# COMPILING OMP

echo "compiling CPU/OMP files..."

for fn in ${OMP_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $CPU_FILES_COMPILER -c $OMP_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $CPU_FILES_FLAGS $BACKEND_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
done

echo ""



# COMPILING GPU

echo "compiling GPU files..."

for fn in ${GPU_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $GPU_FILES_COMPILER -c $GPU_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $GPU_FILES_FLAGS $BACKEND_COMP_FLAGS $HEADER_FLAGS $GPU_WARNING_FLAGS
done

echo ""



# COMPILING MPI

echo "compiling communication/MPI files..."

for fn in ${MPI_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $MPI_FILES_COMPILER -c $MPI_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $MPI_FILES_FLAGS $BACKEND_COMP_FLAGS $HEADER_FLAGS $WARNING_FLAGS
done

echo ""



# LINKING

echo "linking all files..."

# collect list of all objects
OBJECTS=""
OBJECTS+=" $(printf " ${QUEST_OBJ_PREF}%s.o" "${CORE_FILES[@]}")"
OBJECTS+=" $(printf " ${QUEST_OBJ_PREF}%s.o" "${API_FILES[@]}")"
OBJECTS+=" $(printf " ${QUEST_OBJ_PREF}%s.o" "${GPU_FILES[@]}")"
OBJECTS+=" $(printf " ${QUEST_OBJ_PREF}%s.o" "${OMP_FILES[@]}")"
OBJECTS+=" $(printf " ${QUEST_OBJ_PREF}%s.o" "${MPI_FILES[@]}")"

if (( $COMPILE_TESTS == 0 ))
then
    OBJECTS+=" $(printf " ${USER_OBJ_PREF}%s.o" "${USER_FILES[@]}")"
elif (( $COMPILE_TESTS == 1 && $ENABLE_DEPRECATED_API == 0 ))
then
    OBJECTS+=" $(printf " ${TEST_OBJ_PREF}%s.o" "${TEST_MAIN_FILES[@]}")"
    OBJECTS+=" $(printf " ${TEST_OBJ_PREF}%s.o" "${TEST_UTIL_FILES[@]}")"
    OBJECTS+=" $(printf " ${TEST_OBJ_PREF}%s.o" "${TEST_UNIT_FILES[@]}")"
elif (( $COMPILE_TESTS == 1 && $ENABLE_DEPRECATED_API == 1 ))
then
    OBJECTS+=" $(printf " ${TEST_OBJ_PREF}%s.o" "${TEST_DEPR_FILES[@]}")"
    OBJECTS+=" $(printf " ${TEST_OBJ_PREF}%s.o" "${TEST_DEPR_MPI_FILES[@]}")"
fi

# decide executable name
if (( $COMPILE_TESTS == 0 ))
then
    EXEC_FN=$USER_EXEC_FILE
else
    EXEC_FN=$TEST_EXEC_FILE
fi

# link all objects
$LINKER $OBJECTS $ALL_LINK_FLAGS -o $EXEC_FN

echo "${INDENT}compiled executable '${EXEC_FN}'"
echo ""



# CLEAN UP

echo "deleting temporary object files..."

rm *.o

echo "${INDENT}done"
echo ""
