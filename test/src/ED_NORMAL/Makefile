FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

GLOB_INC:=$(shell pkg-config --cflags scifor dmft_ed)

all:  NORMAL SUPERC #NONSU2

NORMAL SUPERC:
	make -C ./ED_$@

clean:
	rm -rf *.mod *.o 
