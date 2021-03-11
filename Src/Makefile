
#F90CFLAGS   = '-Wall -fbounds-check -fcheck-array-temporaries -fcheck=all -g'
F90CFLAGS   = '-Wall -Ofast'
MPIF90CFLAGS = '-Ofast'

LFLAGS   = 
DPLIBS   = 

F90COMPILER =  gfortran
MPIF90COMPILER = mpifort

DEPS = modules.o

HAUSER_OBJECTS =           \
modules.o                  \
YAHFC-Master.o             \
binding-energy.o           \
cleb-lib.o                 \
compound-cs.o              \
decay-chain.o              \
fission-model.o            \
fresco-potentials.o        \
fresco-setup.o             \
Gamma-gamma.o              \
Gauss-quadrature.o         \
Legendre-expand.o          \
GDR.o                      \
HF-denominator.o           \
HF-decay-MC.o              \
HF-primary-decay-setup.o   \
interp.o                   \
level-density.o            \
Moldauer-WF.o              \
optical-setup.o            \
Output-2.o                 \
print-spectra.o            \
print-library-data.o       \
parse-command.o            \
particle-properties.o      \
pre-equilibrium-1.o        \
PREEQ-Samp.o               \
random.o                   \
range_char.o               \
spectrum.o

# this is only made when building YAHFC, not YAHFC-MPI
MPI_OBJECTS = dummy_MPI_libs.o

modules.o : modules.f90
	$(COMPILER) -c $(CFLAGS) $^ -o $@



%.o : %.f90 $(DEPS)
	$(COMPILER) -c $(CFLAGS) $(<F) -o $@

# default: $(HAUSER_OBJECTS) $(MPI_OBJECTS)
# 	$(COMPILER) -o $(EX).x  $(HAUSER_OBJECTS) $(MPI_OBJECTS) $(LFLAGS) $(DPLIBS)

# YAHFC:
# 	cp mpif.h.template mpif.h
# 	$(MAKE) COMPILER=$(F90COMPILER) CFLAGS=$(F90CFLAGS) MPI_OBJECTS=dummy_MPI_libs.o EX=$@
# 	rm mpif.h
# 	cp $@.x $(YAHFC_BIN_DIR)/.

YAHFC:  
	cp mpif.h.template mpif.h
	$(MAKE) $(HAUSER_OBJECTS) $(MPI_OBJECTS) COMPILER=$(F90COMPILER) CFLAGS=$(F90CFLAGS) 
	$(F90COMPILER) -o $@.x $(HAUSER_OBJECTS) $(MPI_OBJECTS)
	rm mpif.h
	cp $@.x $(YAHFC_BIN_DIR)/.

YAHFC-MPI:  
	$(MAKE) $(HAUSER_OBJECTS) COMPILER=$(MPIF90COMPILER) CFLAGS=$(MPIF90CFLAGS) 
	$(MPIF90COMPILER) -o $@.x $(HAUSER_OBJECTS) 
	cp $@.x $(YAHFC_BIN_DIR)/.

clean :
	rm -f *.o *.mod *.x mpif.h

