#======================================================================#
#                                                                      #
#      Makefile -- build the qubit function library                    #
#                                                                      #
#======================================================================#

#
# --- COMMON CONFIG
#

# COMPILER options: GCC, INTEL
COMPILER = GCC
EXE = demo
MY_FILE_NAME = basicTemplate
USE_MPI=0
QUEST_DIR = QUEST

#
# --- compiler
#

ifneq ($(USE_MPI), 1)
	ifeq ($(COMPILER), GCC)
		# GCC compilers
		CC         = gcc
		CFLAGS     = -O2 -std=c99 -mavx -Wall
		CFLAGS_OMP = -fopenmp
	else ifeq ($(COMPILER), INTEL)
		# Intel compilers
		CC         = icc
		CFLAGS     = -O2 -std=c99 -Wall -xAVX -axCORE-AVX2 -restrict
		CFLAGS_OMP = -qopenmp
	else 
		$(error " *** error: invalid compiler")
	endif
else 
	ifeq ($(COMPILER), GCC)
		# GCC compilers
		CC         = mpicc
		CFLAGS     = -O2 -std=c99 -mavx -Wall
		CFLAGS_OMP = -fopenmp
	else ifeq ($(COMPILER), INTEL)
		# Mvapich2
		CC         = mpicc
		CFLAGS     = -O2 -std=c99
		CFLAGS_OMP = -qopenmp
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
OBJ = $(MY_FILE_NAME).o qubits.o
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
default:	demo

demo:		$(OBJ)
		$(CC) $(CFLAGS) $(CFLAGS_OMP) -o $(EXE) $(OBJ) $(LIBS)

.PHONY:		clean veryclean
clean:
		/bin/rm -f *.o demo
veryclean:	clean
		/bin/rm -f *.h~ *.c~ makefile~
