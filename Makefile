.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

## GFORTRAN OPTIONS for saga ##
# module load netCDF-Fortran/4.5.2-gompi-2020a
FC = gfortran
LIB = ${EBROOTNETCDFMINFORTRAN}/lib
INC = ${EBROOTNETCDFMINFORTRAN}/include
# On clusters netCDF-C lives in a separate module path (EBROOTNETCDF).
# On local conda builds it is in the same directory as the Fortran library.
LIBC = $(if ${EBROOTNETCDF},${EBROOTNETCDF}/lib,$(LIB))

objdir = obj
libname = libncio.a

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make gpdd       : compiles the gpdd program gpdd.x"
	@echo " make gpdd_monthly : compiles the gpdd program gpdd_monthly.x"
	@echo " make clean      : cleans object and executable files"
	@echo ""

# Command-line options at make call
debug ?= 0 

FLAGS  = -I$(objdir) -J$(objdir) -I$(INC)
LFLAGS = -L$(LIB) -L$(LIBC) -lnetcdff -lnetcdf -Wl,-rpath,$(LIB) -Wl,-rpath,$(LIBC)

DFLAGS = -O3 -ffree-line-length-none
ifeq ($(debug), 1)
    DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all -ffree-line-length-none
endif

## Individual libraries or modules ##
$(objdir)/ncio.o: ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/massbalance_module.o: massbalance_module.f90 $(objdir)/ncio.o 
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<


## Share library 
$(objdir)/ncio.so: ncio.f90 
	$(FC) -c -shared -fPIC $(DFLAGS) $(FLAGS) -o ncio.so $^

## Static library
lib: $(objdir)/$(libname)

$(objdir)/$(libname): $(objdir)/ncio.o
	ar -rv $@ $^
	ranlib $@

## Complete programs

gpdd: $(objdir)/ncio.o $(objdir)/massbalance_module.o
	$(FC) $(DFLAGS) $(FLAGS) -o gpdd.x $^ gpdd.f90 $(LFLAGS)
	@echo " "
	@echo "    gpdd.x is ready."
	@echo " "

gpdd_monthly: $(objdir)/ncio.o $(objdir)/massbalance_module.o
	$(FC) $(DFLAGS) $(FLAGS) -o gpdd_monthly.x $^ gpdd_monthly.f90 $(LFLAGS)
	@echo " "
	@echo "    gpdd_monthly.x is ready."
	@echo " "

gpdd_monthly_inout: $(objdir)/ncio.o $(objdir)/massbalance_module.o
	$(FC) $(DFLAGS) $(FLAGS) -o gpdd_monthly_inout.x $^ gpdd_monthly_inout.f90 $(LFLAGS)
	@echo " "
	@echo "    gpdd_monthly_inout.x is ready."
	@echo " "

gpdd_monthly_inout_mpi: $(objdir)/ncio.o $(objdir)/massbalance_module.o
	mpifort $(DFLAGS) $(FLAGS) -o gpdd_monthly_inout_mpi.x $^ gpdd_monthly_inout_mpi.f90 $(LFLAGS)
	@echo " "
	@echo "    gpdd_monthly_inout_mpi.x is ready."
	@echo " "

test_massbalance: $(objdir)/massbalance_module.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_massbalance.x $(objdir)/massbalance_module.o test_massbalance.f90
	@echo " "
	@echo "    test_massbalance.x is ready."
	@echo " "

clean:
	rm -f gpdd.x gpdd_monthly.x gpdd_monthly_inout_mpi.x test_massbalance.x $(objdir)/*.o $(objdir)/*.mod $(objdir)/$(libname)

