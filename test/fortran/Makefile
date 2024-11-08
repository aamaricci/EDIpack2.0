FC=mpif90

INC=$(shell pkg-config --cflags edipack2 scifor)
LIB=$(shell pkg-config --libs   edipack2 scifor)

FFLAG = -O2 -ffree-line-length-none -cpp -D_MPI

# Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

all:
	@echo ""	
	$(FC) $(FFLAG)  hm_bethe.f90 -o hm_bethe.exe $(INC) $(LIB) 

clean: 
	@echo "Cleaning:"
	@rm -f *.exe *.mod *.o *~
