
#F77CFLAGS   = -Wall -fbounds-check -fcheck-array-temporaries -fcheck=all -g
#F90CFLAGS   = -Wall -fbounds-check -fcheck-array-temporaries -fcheck=all -g
F77CFLAGS   = -Wall -Wno-maybe-uninitialized -Ofast
F90CFLAGS   = -Wall -Wno-maybe-uninitialized -Ofast
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
parse-command.o            \
particle-properties.o      \
pre-equilibrium-1.o        \
PREEQ-Samp.o               \
random-64.o                \
range_char.o               \
spectrum.o


%.o : %.f90 $(DEPS)
	$(F90COMPILER) -c $(F90CFLAGS) $(<F)  

%.o : %.f
	$(F77COMPILER) -c $(F77CFLAGS) $(<F)  

YAHFC : $(HAUSER_OBJECTS)
	$(F90COMPILER) -o $@.x $(HAUSER_OBJECTS) $(LFLAGS) $(DPLIBS)
	cp $@.x $(YAHFC_BIN_DIR)/.

clean :
	rm *.o *.mod 

