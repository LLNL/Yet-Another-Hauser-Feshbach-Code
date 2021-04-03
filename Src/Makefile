
#F90CFLAGS   = '-Wall -fbounds-check -fcheck-array-temporaries -fcheck=all -g'
F90CFLAGS   = '-Wall -Ofast'
#F90CFLAGS   = '-Wall -Wintrinsic-shadow -Wsurprising -Wextra -Ofast'
#F90CFLAGS   = '-Wall -g -fcheck=bounds -Ofast'
MPIF90CFLAGS = '-Ofast'

LFLAGS   = 
DPLIBS   = 

F90COMPILER =  gfortran
MPIF90COMPILER = mpifort

DEPS = modules.o

HAUSER_OBJECTS =           \
modules.o                  \
YAHFC-Master.o             \
Setup_Energy_bins.o        \
parse-command.o            \
get_binding_energy.o       \
cleb-lib.o                 \
compound-cs.o              \
decay-chain.o              \
fission-model.o            \
fix-incident-energies.o    \
fresco-potentials.o        \
Gamma_gamma.o              \
fit_nuke_Gamma_gamma.o     \
Gauss-quadrature.o         \
Legendre-expand.o          \
GDR.o                      \
EM_str_param.o             \
HF_denominator.o           \
pack_data.o                \
unpack_data.o              \
HF-decay-MC.o              \
HF_primary_decay_setup.o   \
sum_HFden.o                \
interp.o                   \
level-density.o            \
Moldauer-WF.o              \
optical_setup.o            \
state_index.o              \
fresco_make_tco.o          \
run_fresco.o               \
tco_interpolate.o          \
start_IO.o                 \
output_trans_coef.o        \
output_nucleus_data.o      \
check_directories.o        \
print_nuke_data.o          \
print_x_particle_spectra.o \
print_preeq_spectra.o      \
print_direct_spectra.o     \
print_cs_data.o            \
print_elastic.o            \
print_inelastic.o          \
print_channels.o           \
print_channel_gammas.o     \
print_pop_mult_data.o      \
print_primary_decay.o      \
particle-properties.o      \
pre-equilibrium-1.o        \
PREEQ-Samp.o               \
random_32.o                \
random_64.o                \
range_char.o               \
real8_equal.o              \
spectrum.o                 \
remove_states.o

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

