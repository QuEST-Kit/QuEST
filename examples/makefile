
# This makefile builds the QuEST library and links/compiles user code
# While attempting to accomodate as many platforms and compilers,
# unforeseen problems are inevitable: please email anna.brown@oerc.ox.ac.uk
# about any errors or complications, or raise an issue on Github

#======================================================================#
#                                                                      #
#      User Settings                                                   #
#                                                                      #
#======================================================================#

# name of the executable to create
EXE = demo

# space-separated names (no file type) of all user source files (.c or .cpp) in the root directory
SOURCES = tutorialExample

# path to QuEST library from root directory
QUEST_DIR = QuEST

# compiler to use, which should support both C and C++, to be wrapped by GPU/MPI compilers
COMPILER = gcc

# type of above compiler, one of {GNU, INTEL, CLANG}, used for setting compiler flags
COMPILER_TYPE = GNU

# hardwares to target: 1 means use, 0 means don't use
MULTITHREADED = 0
DISTRIBUTED = 0
GPUACCELERATED = 0

# GPU hardware dependent, lookup at https://developer.nvidia.com/cuda-gpus, write without fullstop
GPU_COMPUTE_CAPABILITY = 30

# whether to suppress the below warnings about compiler compatibility
SUPPRESS_WARNING = 0



#======================================================================#
#                                                                      #
#      Adjusting Settings                                              #
#                                                                      #
#======================================================================#

# always allow cleaning without errors or warnings
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), veryclean)

ifneq ($(COMPILER_TYPE), CLANG)
ifneq ($(COMPILER_TYPE), GNU)
ifneq ($(COMPILER_TYPE), INTEL)
    $(error COMPILER_TYPE must be one of CLANG, GNU or INTEL)
endif
endif
endif

ifeq ($(DISTRIBUTED), 1)
ifeq ($(GPUACCELERATED), 1)
    $(error Distributed GPU acceleration not supported)
endif
endif

ifeq ($(MULTITHREADED), 1)
ifeq ($(GPUACCELERATED), 1)
    $(warning GPU acceleration makes no use of multithreading. Disabling the latter...)
    MULTITHREADED = 0
endif
endif

ifeq ($(MULTITHREADED), 1)
ifeq ($(COMPILER_TYPE), CLANG)
    $(warning Clang does not support multithreading. Disabling...)
    MULTITHREADED = 0
endif
endif

ifeq ($(GPUACCELERATED), 1)
ifeq ($(COMPILER_TYPE), CLANG)
ifeq ($(SUPPRESS_WARNING), 0)
    $(info Some versions of Clang are not NVIDIA-GPU compatible. If compilation fails, try Clang 3.7)
endif
endif
endif

ifeq ($(GPUACCELERATED), 1)
ifeq ($(COMPILER_TYPE), GNU)
ifeq ($(SUPPRESS_WARNING), 0)
    $(info On some platforms (e.g. OSX), NVIDIA-GPUs are not compatible with GNU compilers. If compilation fails, try an alternative compiler, like Clang 3.7)
endif
endif
endif

# end of allowed cleaning
endif
endif



#======================================================================#
#                                                                      #
#      Automated Compilation                                           #
#                                                                      #
#======================================================================#

#
# --- libraries
#

LIBS = -lm

ifeq ($(GPUACCELERATED), 1)
    QUEST_INNER_DIR = $(QUEST_DIR)/GPU
else
    QUEST_INNER_DIR = $(QUEST_DIR)/CPU
endif
QUEST_INCLUDE = -I$(QUEST_INNER_DIR)



#
# --- wrapper compilers
#

CUDA_COMPILER = nvcc
MPI_COMPILER = mpicc



#
# --- compiler flags
#

# threading flag
ifeq ($(MULTITHREADED), 1)
    ifeq ($(COMPILER_TYPE), GNU)
        THREAD_FLAGS = -fopenmp
    else ifeq ($(COMPILER_TYPE), INTEL)
        THREAD_FLAGS = -qopenmp # -openmp   # intel make up your mind
    endif
else
    THREAD_FLAGS =
endif

# c
C_CLANG_FLAGS = -O2 -std=c99 -mavx -Wall
C_GNU_FLAGS = -O2 -std=c99 -mavx -Wall $(THREAD_FLAGS)
C_INTEL_FLAGS = -O2 -std=c99 -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch $(THREAD_FLAGS)

# c++
CPP_CLANG_FLAGS = -O2 -std=c++11 -mavx -Wall
CPP_GNU_FLAGS = -O2 -std=c++11 -mavx -Wall $(THREAD_FLAGS)
CPP_INTEL_FLAGS = -O2 -std=c++11 -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch $(THREAD_FLAGS)

# wrappers
CPP_CUDA_FLAGS = -O2 -arch=compute_$(GPU_COMPUTE_CAPABILITY) -code=sm_$(GPU_COMPUTE_CAPABILITY) -ccbin $(COMPILER)

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
endif

#
# --- targets
#

OBJ = QuEST.o mt19937ar.o
ifeq ($(GPUACCELERATED), 1)
    OBJ += QuEST_env_localGPU.o
else ifeq ($(DISTRIBUTED), 1)
    OBJ += QuEST_env_mpi.o
else
    OBJ += QuEST_env_local.o
endif
OBJ += $(addsuffix .o, $(SOURCES))



#
# --- rules
#

# GPU
ifeq ($(GPUACCELERATED), 1)

  %.o: %.c
	$(COMPILER) -x c $(C_FLAGS) $(QUEST_INCLUDE) -c $<
  %.o: $(QUEST_INNER_DIR)/%.c
	$(COMPILER) -x c $(C_FLAGS) -c $<

  %.o: %.cu
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $(QUEST_INCLUDE) $<
  %.o: $(QUEST_INNER_DIR)/%.cu
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $<
	
  %.o: %.cpp
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $(QUEST_INCLUDE) $<
  %.o: $(QUEST_INNER_DIR)/%.cpp
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $<

# distributed
else ifeq ($(DISTRIBUTED), 1)

  %.o: %.c
	OMPI_CC=$(COMPILER) $(MPI_COMPILER) -x c $(C_FLAGS) $(QUEST_INCLUDE) -c $<
  %.o: $(QUEST_INNER_DIR)/%.c
	OMPI_CC=$(COMPILER) $(MPI_COMPILER) -x c $(C_FLAGS) -c $<
	
  %.o: %.cpp
	OMPI_CC=$(COMPILER) $(MPI_COMPILER) $(CPP_FLAGS) $(QUEST_INCLUDE) -c $<
  %.o: $(QUEST_INNER_DIR)/%.cpp
	OMPI_CC=$(COMPILER) $(MPI_COMPILER) $(CPP_FLAGS) -c $<

# CPU
else

  %.o: %.c
	$(COMPILER) -x c $(C_FLAGS) $(QUEST_INCLUDE) -c $<
  %.o: $(QUEST_INNER_DIR)/%.c
	$(COMPILER) -x c $(C_FLAGS) -c $<
	
  %.o: %.cpp
	$(COMPILER) $(CPP_FLAGS) $(QUEST_INCLUDE) -c $<
  %.o: $(QUEST_INNER_DIR)/%.cpp
	$(COMPILER) $(CPP_FLAGS) -c $<

endif



#
# --- build
#

# CUDA
ifeq ($(GPUACCELERATED), 1)

  all:	$(OBJ)
		$(CUDA_COMPILER) $(CPP_CUDA_FLAGS) $(QUEST_INCLUDE) -o $(EXE) $(OBJ) $(LIBS)

# MPI
else ifeq ($(DISTRIBUTED), 1)

  default:	$(EXE)
  $(EXE):	$(OBJ)
			OMPI_CC=$(COMPILER) $(MPI_COMPILER) $(C_FLAGS) $(QUEST_INCLUDE) -o $(EXE) $(OBJ) $(LIBS) 

# C
else

  default:	$(EXE)
  $(EXE):	$(OBJ)
			$(COMPILER) $(C_FLAGS) $(QUEST_INCLUDE) -o $(EXE) $(OBJ) $(LIBS) 

endif



#
# --- clean
#

.PHONY:		clean veryclean
clean:
			/bin/rm -f *.o $(EXE)
veryclean:	clean
			/bin/rm -f *.h~ *.c~ makefile~



#
# --- debug
#
	
print-%:
	@echo $*=$($*)
