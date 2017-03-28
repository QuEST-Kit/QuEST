#======================================================================#
#                                                                      #
#      Makefile -- build the qubit function library                    #
#                                                                      #
#======================================================================#

#
# --- common config
#

EXE = demo
# MODE options: OMP, OMPMPI, GPU
MODE = OMPMPI
# COMPILER options: GCC, INTEL, MPICC
COMPILER = MPICC
MYFILENAME = timingDemo
QUESTDIR = QUEST

#
# --- compiler
#

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
else ifeq ($(COMPILER), MPICC)
  	# Mvapich2
  	CC         = mpicc
  	CFLAGS     = -O2 -std=c99 -g
  	CFLAGS_OMP = -qopenmp
else 
    	$(error " *** error: invalid compiler")
endif

#
# --- libraries
#
LIBS = -lm


#
# --- targets
#
OBJ = $(MYFILENAME).o qubits.o
ifeq ($(MODE), OMPMPI)
	OBJ += qubits_mpi.o
endif

#
# --- rules
#
%.o: %.c
	$(CC) $(CFLAGS) $(CFLAGS_OMP) -c $<

%.o: $(QUESTDIR)/%.c
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
