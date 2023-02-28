FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

GLOB_INC:=$(shell pkg-config --cflags scifor dmft_ed)

all:  NORMAL SUPERC #NONSU2

NORMAL SUPERC:
	@echo " 1 step inside, forward"
	make -C ./ED_$@
	@echo " 1 step inside, backward"
clean:
	rm -rf *.mod *.o 
