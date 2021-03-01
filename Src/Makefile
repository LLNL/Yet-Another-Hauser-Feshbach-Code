
#F77CFLAGS   = -Wall -fbounds-check -fcheck-array-temporaries -fcheck=all -g
#F90CFLAGS   = -Wall -fbounds-check -fcheck-array-temporaries -fcheck=all -g
F77CFLAGS   = -Wall -Ofast
F90CFLAGS   = -Wall -Ofast
#F77CFLAGS   = -Wall -O3
#F90CFLAGS   = -Wall -O3

LFLAGS   = 
DPLIBS   = 

F77COMPILER =  gfortran
F90COMPILER =  gfortran

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
	$(F90COMPILER) -c $(F90CFLAGS) $(<F)  

%.o : %.f $(DEPS)
	$(F77COMPILER) -c $(F77CFLAGS) $(<F)  

default: $(HAUSER_OBJECTS) $(MPI_OBJECTS)
	$(F90) -o $(EX).x $(HAUSER_OBJECTS) $(MPI_OBJECTS) $(LFLAGS) $(DPLIBS)

YAHFC:
	cp mpif.h.template mpif.h
	$(MAKE) F90=$(F90COMPILER) MPI_OBJECTS=dummy_MPI_libs.o EX=$@
	rm mpif.h
	cp $@.x $(YAHFC_BIN_DIR)/.

#YAHFC : $(HAUSER_OBJECTS) $(DUMMY_MPI)
#	cp mpif.h.template mpif.h
#	$(F90COMPILER) -o $@.x $(HAUSER_OBJECTS) $(DUMMY_MPI) $(LFLAGS) $(DPLIBS)
#	rm mpif.h
#	cp $@.x $(YAHFC_BIN_DIR)/.

clean :
	rm *.o *.mod *.x mpif.h

