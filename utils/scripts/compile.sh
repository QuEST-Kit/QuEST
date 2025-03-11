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
COMPILE_MPI=0        # distribution
COMPILE_OPENMP=0     # multithreading
COMPILE_CUDA=0       # GPU acceleration
COMPILE_CUQUANTUM=0  # GPU + cuQuantum

# GPU compute capability
GPU_CC=90

# backend compilers
TESTS_COMPILER=g++
BASE_COMPILER=g++
OMP_COMPILER=g++
MPI_COMPILER=mpic++
GPU_COMPILER=nvcc

# linker
LINKER=g++

# whether to compile unit tests (1) or the below user files (0),
# or the v3 deprecated unit tests (2). when either tests are
# compiled, all user-source related settings are ignored. 
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

# whether to compile cuQuantum (consulted only when COMPILE_CUQUANTUM=1)
# in debug mode, which logs to below file with performance tips and errors
CUQUANTUM_LOG=0
CUQUANTUM_LOG_FN="./custatevec_log.txt"

# external library locations (replace with "." to default)
CUQUANTUM_LIB_DIR="${CUQUANTUM_ROOT}"
CUDA_LIB_DIR="/usr/local/cuda"
OMP_LIB_DIR="/opt/homebrew/opt/libomp"
MPI_LIB_DIR="/opt/homebrew/opt/openmpi"
CATCH_LIB_DIR="tests/deprecated/catch"

# TODO:
# use of 'CATCH_LIB_DIR' above will change when v4 tests 
# switch to using Catch2 as supplied by CMake, rather
# than the hacky use of deprecated v3's single-header  



# CONSTANTS

USER_OBJ_PREF="user_"
QUEST_OBJ_PREF="quest_"
TEST_OBJ_PREF='test_'

INDENT='  '



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
    "operations"
    "paulis"
    "qureg"
    "types"
)

# files in CORE_DIR
CORE_FILES=(
    "accelerator"
    "autodeployer"
    "errors"
    "localiser"
    "memory"
    "parser"
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
    "operations"
    "paulis"
    "qureg"
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
TEST_COMP_FLAGS="-std=c++20 -I${CATCH_LIB_DIR}"

# compiler flags given to deprecated test files
TEST_DEPR_COMP_FLAGS="-std=c++17 -I${TEST_DEPR_CATCH_DIR}"

# compiler flags given to all backend files
BACKEND_COMP_FLAGS='-std=c++17 -O3'

# warning flags which apply to all compiled and linked files including user's
WARNING_FLAGS='-Wall'

# GPU-specific flags
GPU_COMP_FLAGS="-x cu -arch=sm_${GPU_CC} -I${CUDA_LIB_DIR}/include"
GPU_LINK_FLAGS="-L${CUDA_LIB_DIR}/lib -L${CUDA_LIB_DIR}/lib64 -lcudart -lcuda"

if [ $COMPILE_CUQUANTUM == 1 ]
then
    # extend GPU flags if cuQuantum enabled
    GPU_COMP_FLAGS+=" -I${CUQUANTUM_LIB_DIR}/include"
    GPU_LINK_FLAGS+=" -L${CUQUANTUM_LIB_DIR}/lib -L${CUQUANTUM_LIB_DIR}/lib64 -lcustatevec"

    # optional debug logging - will slow down code
    if [ $CUQUANTUM_LOG == 1 ]
    then
        GPU_COMP_FLAGS+=" -DCUSTATEVEC_LOG_LEVEL=5 -DCUSTATEVEC_LOG_FILE=${CUQUANTUM_LOG_FN}"
    fi
fi

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

# define pre-processor macros to indicate deployment mode
MODE_FLAGS="-DCOMPILE_MPI=${COMPILE_MPI} "
MODE_FLAGS+="-DCOMPILE_OPENMP=${COMPILE_OPENMP} "
MODE_FLAGS+="-DCOMPILE_CUDA=${COMPILE_CUDA} "
MODE_FLAGS+="-DCOMPILE_CUQUANTUM=${COMPILE_CUQUANTUM}"

# define pre-processor macros to set qcomp precision
PREC_FLAG="-DFLOAT_PRECISION=${FLOAT_PRECISION}"

# point compilers to QuEST src
HEADER_FLAGS="-I. -I${INCLUDE_DIR}"



# ASSEMBLE FLAGS

echo ""
echo "deployment modes:"

# flags given to every compilation unit
GLOBAL_COMP_FLAGS="${HEADER_FLAGS} ${MODE_FLAGS} ${PREC_FLAG}"

# choose linker flags (extended below)
ALL_LINK_FLAGS="${USER_LINK_FLAGS}"

# choose compiler and flags for CPU/OMP files
if [ $COMPILE_OPENMP == 1 ]
then
    echo "${INDENT}(multithreading enabled)"
    echo "${INDENT}${INDENT}[compiling OpenMP]"
    CPU_FILES_COMPILER=$OMP_COMPILER
    CPU_FILES_FLAGS=$OMP_COMP_FLAGS
    ALL_LINK_FLAGS+=" ${OMP_LINK_FLAGS}"
else
    echo "${INDENT}(multithreading disabled)"
    CPU_FILES_COMPILER=$BASE_COMPILER
    CPU_FILES_FLAGS=''
fi

# choose compiler and flags for GPU files
if [ $COMPILE_CUDA == 1 ]
then
    echo "${INDENT}(GPU-acceleration enabled)"
    echo "${INDENT}${INDENT}[compiling CUDA]"
    GPU_FILES_COMPILER=$GPU_COMPILER
    GPU_FILES_FLAGS=$GPU_COMP_FLAGS
    ALL_LINK_FLAGS+=" ${GPU_LINK_FLAGS}"
    GPU_WARNING_FLAGS="-Xcompiler ${WARNING_FLAGS}"
else
    echo "${INDENT}(GPU-acceleration disabled)"
    GPU_FILES_COMPILER=$BASE_COMPILER
    GPU_FILES_FLAGS=''
    GPU_WARNING_FLAGS=$WARNING_FLAGS
fi

# merely report cuQuantum status
if [ $COMPILE_CUQUANTUM == 1 ]
then
    echo "${INDENT}(cuQuantum enabled)"
    echo "${INDENT}${INDENT}[compiling cuStateVec]"
else
    echo "${INDENT}(cuQuantum disabled)"
fi

# choose compiler and flags for communication files
if [ $COMPILE_MPI == 1 ]
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

# choose linker warning flag (to avoid pass them to nvcc)
if [ "${LINKER}" = "nvcc" ]
then
    ALL_LINK_FLAGS+="-Xcompiler ${WARNING_FLAGS}"
else
    ALL_LINK_FLAGS+=" ${WARNING_FLAGS}"
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
if (( $COMPILE_TESTS == 1 ))
then
    echo "${INDENT}tests compiler and flags:"
    echo "${INDENT}${INDENT}${TESTS_COMPILER} ${TEST_COMP_FLAGS} ${WARNING_FLAGS}"
fi

# deprecated compiler
if (( $COMPILE_TESTS == 2 ))
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
        $COMP -c $USER_DIR/$fn -o ${USER_OBJ_PREF}${fn}.o $FLAG $GLOBAL_COMP_FLAGS $WARNING_FLAGS
    done

    echo ""
fi



# COMPILING TESTS

if (( $COMPILE_TESTS == 1 ))
then

    echo "compiling unit test files:"

    echo "${INDENT}main:"

    for fn in ${TEST_MAIN_FILES[@]}
    do
        echo "${INDENT}${INDENT}${fn}.cpp ..."
        $TESTS_COMPILER -c $TEST_MAIN_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
    done

    echo "${INDENT}utils:"

    for fn in ${TEST_UTIL_FILES[@]}
    do
        echo "${INDENT}${INDENT}${fn}.cpp ..."
        $TESTS_COMPILER -c $TEST_UTIL_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
    done

    echo "${INDENT}unit:"

    for fn in ${TEST_UNIT_FILES[@]}
    do
        echo "${INDENT}${INDENT}${fn}.cpp ..."
        $TESTS_COMPILER -c $TEST_UNIT_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
    done

    echo ""
fi



# COMPILING DEPRECATED TESTS

if (( $COMPILE_TESTS == 2 ))
then
    echo "compiling deprecated test files:"

    for fn in ${TEST_DEPR_FILES[@]}
    do
        echo "${INDENT}${fn}.cpp ..."
        $TESTS_COMPILER -c $TEST_DEPR_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_DEPR_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
    done

    for fn in ${TEST_DEPR_MPI_FILES[@]}
    do
        echo "${INDENT}${fn}.cpp ..."
        $MPI_FILES_COMPILER -c $TEST_DEPR_DIR/$fn.cpp -o ${TEST_OBJ_PREF}${fn}.o $TEST_DEPR_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
    done

    echo ""
fi



# COMPILING CORE

echo "compiling core files in C++"

for fn in ${CORE_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $BASE_COMPILER -c $CORE_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $BACKEND_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
done

echo ""



# COMPILING API

echo "compiling API files in C++:"

for fn in ${API_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $BASE_COMPILER -c $API_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $BACKEND_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
done

echo ""



# COMPILING OMP

echo "compiling CPU/OMP files..."

for fn in ${OMP_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $CPU_FILES_COMPILER -c $OMP_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $CPU_FILES_FLAGS $BACKEND_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
done

echo ""



# COMPILING GPU

echo "compiling GPU files..."

for fn in ${GPU_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $GPU_FILES_COMPILER -c $GPU_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $GPU_FILES_FLAGS $BACKEND_COMP_FLAGS $GLOBAL_COMP_FLAGS $GPU_WARNING_FLAGS
done

echo ""



# COMPILING MPI

echo "compiling communication/MPI files..."

for fn in ${MPI_FILES[@]}
do
    echo "${INDENT}${fn}.cpp ..."
    $MPI_FILES_COMPILER -c $MPI_DIR/$fn.cpp -o ${QUEST_OBJ_PREF}${fn}.o $MPI_FILES_FLAGS $BACKEND_COMP_FLAGS $GLOBAL_COMP_FLAGS $WARNING_FLAGS
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
elif (( $COMPILE_TESTS == 1 ))
then
    OBJECTS+=" $(printf " ${TEST_OBJ_PREF}%s.o" "${TEST_MAIN_FILES[@]}")"
    OBJECTS+=" $(printf " ${TEST_OBJ_PREF}%s.o" "${TEST_UTIL_FILES[@]}")"
    OBJECTS+=" $(printf " ${TEST_OBJ_PREF}%s.o" "${TEST_UNIT_FILES[@]}")"
elif (( $COMPILE_TESTS == 2 ))
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
