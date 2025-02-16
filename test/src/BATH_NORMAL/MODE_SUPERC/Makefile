EXE=ed_normal_superc
FC=mpif90
PLAT=gnu
TSTBIN=../../../bin/NORMAL/SUPERC
OBJS=../../ASSERTING.o

#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################
GLOB_INC:=$(shell pkg-config --cflags edipack2 scifor)
GLOB_LIB:=$(shell pkg-config --libs   edipack2 scifor | sed  "s/;/ /g"  | sed 's/\\/  /g')


ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz -coverage
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
FPPFLAG =-fpp
endif
ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none -coverage
DFLAG = -O2 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
FPPFLAG =-cpp -D_MPI
endif

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90 .o

all: FLAG:=${FFLAG} ${FPPFLAG}
all:
	$(FC) $(FLAG) $(OBJS) $(EXE).f90 -o $(TSTBIN)/$(EXE).x -I../../ ${GLOB_INC} ${GLOB_LIB}
	@cp *.check $(TSTBIN)/
	@cp *.restart $(TSTBIN)/
	@cp *.in    $(TSTBIN)/
	@echo ""

debug: FLAG:=${DFLAG} ${FPPFLAG}
debug:
	$(FC) $(FLAG) $(OBJS) $(EXE).f90 -o $(TSTBIN)/$(EXE).x -I../../ ${GLOB_INC} ${GLOB_LIB}
	@cp *.check $(TSTBIN)/
	@cp *.restart $(TSTBIN)/
	@cp *.in    $(TSTBIN)/
	@echo ""

clean: 
	@echo "Cleaning bath:NORMAL mode:SUPERC"
	@rm -rf *.mod *.o *~
	@rm -fv  $(TSTBIN)/$(EXE)



#########################################################################
