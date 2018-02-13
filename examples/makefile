#======================================================================#
#                                                                      #
#      Makefile -- build the qubit function library                    #
#                                                                      #
#======================================================================#

#
# --- COMMON CONFIG
#

# COMPILER options: GNU, INTEL
COMPILER = GNU
EXE = demo
MY_FILE_NAME = basicTemplate
USE_MPI=0
USE_OPENMP=0
QUEST_DIR = QuEST

#
# --- compiler
#

ifeq ($(COMPILER), GNU)
	CFLAGS_OMP=-fopenmp
else 
	CFLAGS_OMP=-openmp
endif

ifneq ($(USE_OPENMP), 1)
	# disable OPENMP
	CFLAGS_OMP=
endif

ifneq ($(USE_MPI), 1)
	ifeq ($(COMPILER), GNU)
		# GCC compilers
		CC         = gcc
		CFLAGS     = -O2 -std=c99 -mavx -Wall
	else ifeq ($(COMPILER), INTEL)
		# Intel compilers
		CC         = icc
		CFLAGS     = -O2 -fprotect-parens -std=c99 -Wall -xAVX -axCORE-AVX2 -restrict
	else 
		$(error " *** error: invalid compiler")
	endif
else 
	ifeq ($(COMPILER), GNU)
		# GCC compilers
		CC         = mpicc
		CFLAGS     = -O2 -std=c99 -mavx -Wall
	else ifeq ($(COMPILER), INTEL)
		# Mvapich2
		CC         = mpicc
		CFLAGS     = -O2 -fprotect-parens -std=c99
	else 
		$(error " *** error: invalid compiler")
	endif
endif

#
# --- libraries
#
LIBS = -lm


#
# --- targets
#
OBJ = $(MY_FILE_NAME).o qubits.o mt19937ar.o
ifneq ($(USE_MPI), 0)
	OBJ += qubits_env_mpi.o
else
	OBJ += qubits_env_local.o
endif

#
# --- rules
#
%.o: %.c
	$(CC) $(CFLAGS) $(CFLAGS_OMP) -c $<

%.o: $(QUEST_DIR)/%.c
	$(CC) $(CFLAGS) $(CFLAGS_OMP) -c $<


#
# --- build
#
default:	$(EXE)

$(EXE):		$(OBJ)
		$(CC) $(CFLAGS) $(CFLAGS_OMP) -o $(EXE) $(OBJ) $(LIBS)

.PHONY:		clean veryclean
clean:
		/bin/rm -f *.o $(EXE)
veryclean:	clean
		/bin/rm -f *.h~ *.c~ makefile~
