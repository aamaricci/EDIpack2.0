FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

GLOB_INC:=$(shell pkg-config --cflags scifor )

OBJS=ASSERTING.o

all: $(OBJS) NORMAL REPLICA GENERAL # HYBRID
	@echo "all done"


NORMAL REPLICA GENERAL:
	@echo " pre $@ comp"
	make -C ./ED_$@
	@echo "post $@ comp"

#HYBRID:
#	make -C ./ED_$@

.f90.o:
	@echo "pre assert comp" #delete
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}
	@echo "post assert comp" #delete

clean:
	rm -rf *.mod *.o 

