# ------------------------------------------------
# LSCSQ only requires a FORTRAN compiler and  the C and FORTAN netCDF libraries
# Compiler: either INTEL ifort, NVIDIA nvfortran, or gfortran
#
# One can specify the compiler by setting the F90 environment variable,
# especially if it is a different compiler than ifort, nvfortran, and gfortran:
#
#     F90=pgf90 make
#
# On the PPPL cluster: simply load the appropriate modules
#         module load intel netcdf-c netcdf-fortran


# Name of the program
EXE = lsc_driver
LIBNAME= lscsq
PREFIX := $(or $(PREFIX),$(shell pwd)/build)

ifndef F90
# If the variable F90 is not define, look for a valid Fortran compiler
# starting with Intel ifort, then Nvidia nvfortran, then gfortran
   ifneq (, $(shell which ifort))
      FC := ifort
      $(info "Using the INTEL ifort compiler")
   else ifeq (, $(shell which nvfortran))
      FC := nvfortran
      $(info "Using the NVIDIA nvfortran compiler")
   else ifeq (, $(shell which gfortran))
      FC := gfortran
      $(info "Using the GNU gfortran compiler")
   else
      $(error "!!! Cannot find a valid FORTRAN compiler")
   endif
else
   FC := $(F90)
endif

ifneq (,$(findstring ifort, $(FC)))
	CFLAGS  := -O -w #-mp  #-cxxlib-gcc -gcc-version=333
	DCFLAGS := -g  #-cxxlib-gcc -gcc-version=333

	FFLAGS  := -g -O3 -w -mcmodel small -restrict -align -debug
	DFFLAGS := -g -O0 -w -mcmodel small -restrict 
else ifneq (,$(findstring nvfortran, $(FC)))
        CFLAGS  := -O
        DCFLAGS := -g
        FFLAGS  := -g -O3
        DFFLAGS := -g -pg
else ifneq (,$(findstring gfortran, $(FC)))
	CSPECIAL= -fno-align-commons -fbounds-check -fbacktrace -fPIC -fcheck=pointer
        CFLAGS  := -O
        DCFLAGS := -g
        FFLAGS  := -g -std=f2003 -O3 -fdefault-real-8
        DFFLAGS := -g -Wall -std=f2003 -finit-real=zero -fcheckall -fbacktrace \
		   -ffpe-trap=invalid -ffpe-summary=invalid -fdebug-aux-vars
else
        CFLAGS  := -O
        DCFLAGS := -g
        FFLAGS  := -O
        DFFLAGS := -g
endif

ifdef DEBUG
	FFLAGS := $(DFFLAGS)
	CFLAGS := $(DCFLAGS)
endif


# compile all .f90 / .c files in directory
ALL_F := $(subst .f90,.o, $(wildcard *.f90))
ALL_C := $(subst .c,.o, $(wildcard *.c))
M_OBJ := driver.o
OBJECTS := lscsq_mod.o nml2lsc.o fits.o plasejv.o  \
          grids.o matr.o power.o raystore.o grap.o grapgrd.o ql.o raytrace.o cycle.o fe.o       \
          jrf.o rayini.o wr2tsc.o lscsq_main.o 


# Includes, F90 Use files
INCLUDES = -I./ -I`nf-config --includedir`

L_NETCDF  =  `nc-config --libs` `nf-config --flibs`


%.o: %.f90
	$(NICE) $(FC) $(FFLAGS) -c $(INCLUDES) -o $@ $<

%.o: %.f
	$(FC) $(FFLAGS) -c $(INCLUDES) -o $@ $<

%.o: %.F
	$(FC) $(FFLAGS) -c $(INCLUDES) -o $@ $<

all:  $(EXE) lib

lib: lib$(LIBNAME).a

$(EXE): $(OBJECTS) $(M_OBJ)
	$(NICE) $(FC) $(FFLAGS) -o $@ $^ $(L_NETCDF) $(L_MKL)
#	chmod 755 $(EXE)

lib$(LIBNAME).a: $(OBJECTS)
	ar rcs $@ $^

clean:
	rm -f *.o $(EXE) *.mod *.a

install:
	@echo "Installing here $(PREFIX)"
	@mkdir -p $(PREFIX)/{bin,include,lib}
	@cp $(EXE) $(PREFIX)/bin
	@cp *.mod $(PREFIX)/include
	@cp lib$(LIBNAME).a $(PREFIX)/lib

uninstall:
	@echo "Uninstalling from $(PREFIX)"
	@rm -rf $(PREFIX)

.PHONY: all lib clean install uninstall $(EXE)
