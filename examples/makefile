# DEPRECATION WARNING!
# QuEST has switched to using the cmake build system. 
# See https://github.com/QuEST-Kit/QuEST/tree/master/examples#compiling
# for more information.
# This makefile is included as a backup in case of problems switching to the
# new system but is less rigorously tested than the cmake system. Use at your own risk. 

# This makefile builds the QuEST library and links/compiles user code
# While attempting to accomodate as many platforms and compilers,
# unforeseen problems are inevitable: please email quest@materials.ox.ac.uk
# about any errors or complications, or raise an issue on Github

# author: Tyson Jones, based on...
# author: Ania Brown's individual makefiles

#======================================================================#
#                                                                      #
#      User settings                                                   #
#                                                                      #
#======================================================================#

# name of the executable to create
EXE = demo

# space-separated names (no file type) of all user source files (.c or .cpp) in the root directory
SOURCES = tutorial_example

# path to QuEST library from root directory
QUEST_DIR = QuEST

# compiler to use, which should support both C and C++, to be wrapped by GPU/MPI compilers
COMPILER = gcc

# type of above compiler, one of {GNU, INTEL, CLANG, MSVC}, used for setting compiler flags
COMPILER_TYPE = GNU

# only for WINDOWS: whether OS is 32-bit (x86) or 64-bit (x64). Choose {32, 64}
WINDOWS_ARCH = 32

# hardwares to target: 1 means use, 0 means don't use
MULTITHREADED = 0
DISTRIBUTED = 0
GPUACCELERATED = 0

# GPU hardware dependent, lookup at https://developer.nvidia.com/cuda-gpus, write without fullstop
GPU_COMPUTE_CAPABILITY = 30

# whether to suppress the below warnings about compiler compatibility
SUPPRESS_WARNING = 0

# whether to use single, double or quad floating point precision in the state-vector {1,2,4}
PRECISION = 2

#======================================================================#
#                                                                      #
#      Checking user settings                                          #
#                                                                      #
#======================================================================#

# suppresses all non-gcc output, useful for calling scripts
SILENT = 0

# always allow cleaning without errors or warnings
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), veryclean)
ifneq ($(SILENT), 1)

    # check $COMPILER_TYPE is correct
    ifneq ($(COMPILER_TYPE), CLANG)
    ifneq ($(COMPILER_TYPE), GNU)
    ifneq ($(COMPILER_TYPE), INTEL)
	  ifneq ($(COMPILER_TYPE), MSVC)
        $(error COMPILER_TYPE must be one of CLANG, GNU or INTEL)
    endif
    endif
    endif
    endif

    # distributed GPU not supported
    ifeq ($(DISTRIBUTED), 1)
    ifeq ($(GPUACCELERATED), 1)
        $(error Distributed GPU acceleration not supported)
    endif
    endif

    # GPU doesn't use threading
    ifeq ($(MULTITHREADED), 1)
    ifeq ($(GPUACCELERATED), 1)
        $(warning GPU acceleration makes no use of multithreading. Disabling the latter...)
        override MULTITHREADED = 0
    endif
    endif
	
    # CLANG compilers don't support threading
    ifeq ($(MULTITHREADED), 1)
    ifeq ($(COMPILER_TYPE), CLANG)
        $(warning Clang does not support multithreading. Disabling...)
        override MULTITHREADED = 0
    endif
    endif
	
	# check PRECISION is valid
    ifneq ($(PRECISION), 1)
    ifneq ($(PRECISION), 2)
    ifneq ($(PRECISION), 4)
        $(error PRECISION must be set to 1, 2 or 4)
    endif
    endif
    endif
	
	# GPU does not support quad precision
    ifeq ($(PRECISION), 4)
    ifeq ($(GPUACCELERATED), 1)
    $(warning GPUs do not support quad precision. Setting PRECISION=2...)
    override PRECISION = 2	
    endif
    endif
	
    # NVCC doesn't support new CLANG compilers
    ifeq ($(GPUACCELERATED), 1)
    ifeq ($(COMPILER_TYPE), CLANG)
    ifeq ($(SUPPRESS_WARNING), 0)
        $(info Some versions of Clang are not NVIDIA-GPU compatible. If compilation fails, try Clang 3.7)
    endif
    endif
    endif

    # NVCC doesn't support GNU compilers on OSX
    ifeq ($(GPUACCELERATED), 1)
    ifeq ($(COMPILER_TYPE), GNU)
    ifeq ($(SUPPRESS_WARNING), 0)
        $(info On some platforms (e.g. OSX), NVIDIA-GPUs are not compatible with GNU compilers. If compilation fails, try an alternative compiler, like Clang 3.7)
    endif
    endif
    endif
    
    # Windows users must set WINDOWS_ARCH as {32, 64}
    ifeq ($(COMPILER_TYPE), MSVC)
    ifneq ($(WINDOWS_ARCH), 32)
    ifneq ($(WINDOWS_ARCH), 64)
        $(error When compiling with MSVC, WINDOWS_ARCH must be 32 or 64)
    endif
    endif
    endif

# end of allowed cleaning
endif
endif
endif


#======================================================================#
#                                                                      #
#     Compilation                                                      #
#                                                                      #
#======================================================================#


#
# --- libraries
#

ifeq ($(COMPILER_TYPE), MSVC)
    LIBS =
else
    LIBS = -lm
endif


#
# --- source and include paths
#

QUEST_INCLUDE_DIR = ${QUEST_DIR}/include
QUEST_SRC_DIR = ${QUEST_DIR}/src

QUEST_COMMON_DIR = $(QUEST_SRC_DIR)
ifeq ($(GPUACCELERATED), 1)
    QUEST_INNER_DIR = $(QUEST_SRC_DIR)/GPU
else
    QUEST_INNER_DIR = $(QUEST_SRC_DIR)/CPU
endif
QUEST_INCLUDE = -I${QUEST_INCLUDE_DIR} -I$(QUEST_INNER_DIR) -I$(QUEST_COMMON_DIR)


#
# --- wrapper compilers
#

CUDA_COMPILER = nvcc
MPI_COMPILER = mpicc



#
# --- compiler flags
#

# note:
#	several flag names depend not just on the compiler type, but also compiler version
#	the user should update these below. For example:
#
#	- GNU C++ compilers of version < 4.7 use -std=c++0x instead of -std=c++11
#	- INTEL compilers of version < ? use -openmp instead of -qopenmp
#	- INTEL compilers of version < ? won't recognise -diad-disable and -cpu-dispatch
#	- CLANG compilers don't support openmp (threading) at all


# threading flag
ifeq ($(MULTITHREADED), 1)
    ifeq ($(COMPILER_TYPE), GNU)
        THREAD_FLAGS = -fopenmp
    else ifeq ($(COMPILER_TYPE), INTEL)
        THREAD_FLAGS = -qopenmp
    else ifeq ($(COMPILER_TYPE), MSVC)
        THREAD_FLAGS = -openmp
    endif
else
    THREAD_FLAGS =
endif

# windows architecture flag
ifeq ($(WINDOWS_ARCH), 32)
    ARCH_FLAG = X86
else
    ARCH_FLAG = X64
endif

# c
C_CLANG_FLAGS = -O2 -std=c99 -mavx -Wall -DQuEST_PREC=$(PRECISION)
C_GNU_FLAGS = -O2 -std=c99 -mavx -Wall -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS)
C_INTEL_FLAGS = -O2 -std=c99 -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS)
C_MSVC_FLAGS = -O2 -EHs -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS) -nologo -DDWIN$(WINDOWS_ARCH) -D_WINDOWS -Fo$@

# c++
CPP_CLANG_FLAGS = -O2 -std=c++11 -mavx -Wall -DQuEST_PREC=$(PRECISION)
CPP_GNU_FLAGS = -O2 -std=c++11 -mavx -Wall -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS)
CPP_INTEL_FLAGS = -O2 -std=c++11 -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS)
CPP_MSVC_FLAGS = -O2 -EHs -std:c++latest -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS) -nologo -DDWIN$(WINDOWS_ARCH) -D_WINDOWS -Fo$@

# wrappers
CPP_CUDA_FLAGS = -O2 -arch=compute_$(GPU_COMPUTE_CAPABILITY) -code=sm_$(GPU_COMPUTE_CAPABILITY) -DQuEST_PREC=$(PRECISION)

# choose c/c++ flags based on compiler type
ifeq ($(COMPILER_TYPE), CLANG)
    C_FLAGS = $(C_CLANG_FLAGS)
    CPP_FLAGS = $(CPP_CLANG_FLAGS)
else ifeq ($(COMPILER_TYPE), GNU)
    C_FLAGS = $(C_GNU_FLAGS)
    CPP_FLAGS = $(CPP_GNU_FLAGS)
else ifeq ($(COMPILER_TYPE), INTEL)
    C_FLAGS = $(C_INTEL_FLAGS)
    CPP_FLAGS = $(CPP_INTEL_FLAGS)
else ifeq ($(COMPILER_TYPE), MSVC)
    C_FLAGS = $(C_MSVC_FLAGS)
    CPP_FLAGS = $(CPP_MSVC_FLAGS)
	
	# must specify machine type on Windows
    CPP_CUDA_FLAGS := $(CPP_CUDA_FLAGS) -m=$(WINDOWS_ARCH) -DDWIN$(WINDOWS_ARCH)
endif



#
# --- compiler mode and linker flags 
#

# format args based on compiler type
ifeq ($(COMPILER_TYPE), MSVC)
    C_MODE = 
    LINKER = link.exe
    LINK_FLAGS := -SUBSYSTEM:CONSOLE -nologo -MACHINE:$(ARCH_FLAG)
		# note MSVC linker does not receive thread flags
    	
    # must forward linker flags from NVCC to link.exe on Windows
    ifeq ($(GPUACCELERATED), 1)
        LINK_FLAGS := -o $(EXE).exe $(foreach option, $(LINK_FLAGS), -Xlinker $(option))
    else
        LINK_FLAGS := -out:$(EXE).exe $(LINK_FLAGS)
    endif
    
    MPI_VARS = 
else
    C_MODE = -x c
    LINKER = $(COMPILER)
    LINK_FLAGS := -o $(EXE) $(THREAD_FLAGS)
	
    MPI_VARS = I_MPI_CC=$(COMPILER) OMPI_CC=$(COMPILER) MPICH_CC=$(COMPILER)
endif

# prepare compiler + args shortcut
ifeq ($(DISTRIBUTED), 1)
    COMP_CMD = $(MPI_VARS) $(MPI_COMPILER)
    LINK_CMD = $(MPI_VARS) $(MPI_COMPILER)
else
    COMP_CMD = $(COMPILER)
    LINK_CMD = $(LINKER)
endif



#
# --- targets
#

OBJ = QuEST.o QuEST_validation.o QuEST_common.o QuEST_qasm.o mt19937ar.o
ifeq ($(GPUACCELERATED), 1)
    OBJ += QuEST_gpu.o
else ifeq ($(DISTRIBUTED), 1)
    OBJ += QuEST_cpu.o QuEST_cpu_distributed.o
else
    OBJ += QuEST_cpu.o QuEST_cpu_local.o
endif
OBJ += $(addsuffix .o, $(SOURCES))



#
# --- building
#

# notes:
#	- if $SOURCES appear as both c and c++ files, the c files will be compiled
#	- CUDA won't compile .c files ($COMPILER will), only .cpp and .cu
#	- MPICC will compile .c and .cpp files (wrapping $COMPILER)


# first, let NVCC compile any GPU sources
ifeq ($(GPUACCELERATED), 1)

  # final -o to force NVCC to use '.o' extension even on Windows
  %.o: %.cu
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) -ccbin $(COMPILER) $(QUEST_INCLUDE) -o $@ $<
  %.o: $(QUEST_INNER_DIR)/%.cu
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) -ccbin $(COMPILER) $(QUEST_INCLUDE) -o $@ $<

endif

# remaining files compiled for CPU
# C
%.o: %.c
	$(COMP_CMD) $(C_MODE) $(C_FLAGS) $(QUEST_INCLUDE) -c $<
%.o: $(QUEST_INNER_DIR)/%.c
	$(COMP_CMD) $(C_MODE) $(C_FLAGS) $(QUEST_INCLUDE) -c $<
%.o: $(QUEST_COMMON_DIR)/%.c
	$(COMP_CMD) $(C_MODE) $(C_FLAGS) $(QUEST_INCLUDE) -c $<
	
# C++
%.o: %.cpp
	$(COMP_CMD) $(CPP_FLAGS) $(QUEST_INCLUDE) -c $<
%.o: $(QUEST_INNER_DIR)/%.cpp
	$(COMP_CMD) $(CPP_FLAGS) -c $<
%.o: $(LINK_DIR)/%.cpp
	$(COMP_CMD) $(CPP_FLAGS) $(QUEST_INCLUDE) -c $<



#
# --- linking
#

# CUDA
ifeq ($(GPUACCELERATED), 1)

  # a dirty hack to silence cl when NVCC linking
  # (https://stackoverflow.com/questions/61178458/force-nvcc-straight-to-linking-phase)
  SHUTUP := 
  ifeq ($(COMPILER_TYPE), MSVC)
      SHUTUP := -Xcompiler 2>nul:
  endif

  all:	$(OBJ)
	$(CUDA_COMPILER) $(SHUTUP) $(CPP_CUDA_FLAGS) $(OBJ) $(LIBS) $(LINK_FLAGS)

# C and C++
else

  default:	$(EXE)
  $(EXE):	$(OBJ)
			$(LINK_CMD) $(OBJ) $(LIBS) $(LINK_FLAGS)

endif




#
# --- clean
#

# use OS delete command, inferred from compiler type
# (this incorrectly assumes intel compilers aren't run on Windows, 
#  despite them being far better there than MSVC! Fix this in Cmake)
ifeq ($(COMPILER_TYPE), MSVC)
    REM = del
    EXE_FN = $(EXE).exe
else
    REM = /bin/rm -f
    EXE_FN = $(EXE)
endif

.PHONY:		tidy clean veryclean
tidy:
			$(REM) *.o *.lib *.exp
clean:	tidy
			$(REM) $(EXE_FN)
veryclean:	clean
			$(REM) *.h~ *.c~ makefile~



#
# --- debug
#
	
print-%:
	@echo $*=$($*)

getvalue-%:
	@echo $($*)



