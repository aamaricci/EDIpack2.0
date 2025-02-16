FC=mpif90
FFLAG = -O2 -ffree-line-length-none -coverage
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

GLOB_INC:=$(shell pkg-config --cflags scifor )

OBJS=ASSERTING.o

all: $(OBJS) NORMAL REPLICA GENERAL # HYBRID
	@echo "all done"


NORMAL REPLICA GENERAL: #HYBRID
	make -C ./BATH_$@

.f90.o:
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}

clean:
	rm -rf *.mod *.o 

