#======================================================================#
#                                                                      #
#      Makefile -- build the qubit function library                    #
#                                                                      #
#======================================================================#

#
# --- compiler
#
# default is GCC
ifndef COMPILER
  COMPILER = GCC
endif

# GCC compilers
ifeq ($(COMPILER), GCC)
  CC         = gcc
  CFLAGS     = -O2 -std=c99 -mavx -Wall
  CFLAGS_OMP = -fopenmp
else
  # Intel compilers
  ifeq ($(COMPILER), INTEL)
    CC         = icc
    CFLAGS     = -O2 -std=c99 -Wall -xAVX -axCORE-AVX2 -restrict
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
EXE = demo
OBJ = timingDemo.o qubits.o 

#
# --- rules
#
.SUFFIXES:
.SUFFIXES: .c .h .o
.c.o:
		$(CC) $(CFLAGS) $(CFLAGS_OMP) -c $<


#
# --- build
#
default:	demo

help:
		@echo
		@echo " usage:    make"
		@echo "           make COMPILER=INTEL"
		@echo " defaults: COMPILER=GCC"
		@echo

demo:		$(OBJ)
		$(CC) $(CFLAGS) $(CFLAGS_OMP) -o $(EXE) $(OBJ) $(LIBS)

.PHONY:		clean veryclean
clean:
		/bin/rm -f *.o demo
veryclean:	clean
		/bin/rm -f *.h~ *.c~ makefile~
