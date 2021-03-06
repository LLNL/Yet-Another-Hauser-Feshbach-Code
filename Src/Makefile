
#F90CFLAGS   = -Wall -fbounds-check -fcheck-array-temporaries -fcheck=all -g
#CFLAGS   = -Wall -Ofast
F90CFLAGS   = "-Wall -Ofast"

LFLAGS   = 
DPLIBS   = 

F90COMPILER =  gfortran
MPI90COMPILER = mpifort

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

%.o : %.f90 $(DEPS)
	$(COMPILER) -c $(CFLAGS) $(<F)  


default: $(HAUSER_OBJECTS) $(MPI_OBJECTS)
	$(COMPILER) -o $(EX).x  $(HAUSER_OBJECTS) $(MPI_OBJECTS) $(LFLAGS) $(DPLIBS)

YAHFC:
	cp mpif.h.template mpif.h
	$(MAKE) COMPILER=$(F90COMPILER) CFLAGS=$(F90CFLAGS) MPI_OBJECTS=dummy_MPI_libs.o EX=$@
	rm mpif.h
	cp $@.x $(YAHFC_BIN_DIR)/.

YAHFC-MPI:
	cp mpif.h.template mpif.h
	$(MAKE) COMPILER=$(MPIF90COMPILER) CFLAGS=$(F90CFLAGS) MPI_OBJECTS= EX=$@
	rm mpif.h
	cp $@.x $(YAHFC_BIN_DIR)/.

#YAHFC : $(HAUSER_OBJECTS) $(DUMMY_MPI)
#	cp mpif.h.template mpif.h
#	$(F90COMPILER) -o $@.x $(HAUSER_OBJECTS) $(DUMMY_MPI) $(LFLAGS) $(DPLIBS)
#	rm mpif.h
#	cp $@.x $(YAHFC_BIN_DIR)/.

clean :
	rm *.o *.mod *.x mpif.h

