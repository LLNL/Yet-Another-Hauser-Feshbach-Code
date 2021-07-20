#
#---   Uncomment lines to perform desired compilation, and set compilers for your system
#
#---   pre-processor variable USE_MPI is used to compile section of the code
#---   with MPI calls. USE_MPI=0 means no MPI, USE_MPI=1 is with MPI
#
#----   Set compilers - default is gfortran for serial and mpifort for mpi 
#
F90COMPILER =  gfortran
MPIF90COMPILER = mpifort
#
#************************************************************************************************
#---   Standard set up with optimization -------------------------------------------------------*
#************************************************************************************************
F90CFLAGS   = '-Wall -Ofast -cpp -DUSE_MPI=0'
#
#----   Bounds checking, all checks and debugging with no optimization   ------------------------
#
#F90CFLAGS   = '-Wall -fbounds-check -fcheck-array-temporaries -fcheck=all -g -cpp -DUSE_MPI=0'
#
#----   Bounds checking, with optimization     --------------------------------------------------
#
#F90CFLAGS   = '-Wall -g -fcheck=bounds -Ofast -cpp -DUSE_MPI=0'
#
#************************************************************************************************
#----   Flags for compilation with MPI         -------------------------------------------------*
#************************************************************************************************
MPIF90CFLAGS = '-Ofast -cpp -DUSE_MPI=1'
#
#----   Library flags - currently none         --------------------------------------------------
#
LFLAGS   = 
DPLIBS   = 
#
#----   Dependencies - all modules defined in modules.f90, which must be compiled first and
#----   everything must be recompiled if it is edited
#
DEPS = modules.o
#
#----    Objects comprised of YAHFC code system
#
HAUSER_OBJECTS =                 \
modules.o                        \
YAHFC-Master.o                   \
fail_safe_check.o                \
Setup_Energy_bins.o              \
parse-command.o                  \
get_binding_energy.o             \
cleb-lib.o                       \
compound-cs.o                    \
decay-chain.o                    \
fission-model.o                  \
fix-incident-energies.o          \
fresco-potentials.o              \
Gamma_gamma.o                    \
fit_nuke_Gamma_gamma.o           \
Gauss-quadrature.o               \
eig_eis.o                        \
Legendre-expand.o                \
gamma_functions.o                \
EM_str_param.o                   \
HF_denominator.o                 \
pack_data.o                      \
unpack_data.o                    \
HF-decay-MC.o                    \
HF_primary_decay_setup.o         \
interp.o                         \
level-density.o                  \
Moldauer-WF.o                    \
optical_setup.o                  \
state_index.o                    \
fresco_make_tco.o                \
run_fresco.o                     \
tco_interpolate.o                \
start_IO.o                       \
output_trans_coef.o              \
output_nucleus_data.o            \
check_directories.o              \
print_nuke_data.o                \
print_x_particle_spectra.o       \
print_preeq_spectra.o            \
print_direct_spectra.o           \
print_cs_data.o                  \
print_elastic.o                  \
print_inelastic.o                \
print_channels.o                 \
print_channel_gammas.o           \
print_channel_gammas_all_ein.o   \
print_pop_mult_data.o            \
print_primary_decay.o            \
particle-properties.o            \
pre-equilibrium-exciton.o        \
PREEQ-Samp.o                     \
random_32.o                      \
random_64.o                      \
range_char.o                     \
real8_equal.o                    \
memory_used.o                    \
spectrum.o                       \
remove_states.o                  \
read_saved_parameters.o
#
#----   Rules to compile .o files 
#
%.o : %.f90 $(DEPS)
	$(COMPILER) -c $(CFLAGS) $(<F) -o $@
#
#----   Rules to compile serial verison, creates YAHFC.x and copies to YAHFC bin directory 
#
YAHFC:  
	$(MAKE) $(HAUSER_OBJECTS) COMPILER=$(F90COMPILER) CFLAGS=$(F90CFLAGS) 
	$(F90COMPILER) -o $@.x $(HAUSER_OBJECTS) 
	cp $@.x $(YAHFC_BIN_DIR)/.
#
#----   Rules to compile MPI verison, creates YAHFC.x and copies to YAHFC bin directory 
#
YAHFC-MPI:  
	$(MAKE) $(HAUSER_OBJECTS) COMPILER=$(MPIF90COMPILER) CFLAGS=$(MPIF90CFLAGS) 
	$(MPIF90COMPILER) -o $@.x $(HAUSER_OBJECTS) 
	cp $@.x $(YAHFC_BIN_DIR)/.
#
#----   Rules to make modules
#
modules.o : modules.f90
	$(COMPILER) -c $(CFLAGS) $^ -o $@
#
#----   Rules to clean up
#
clean :
	rm -f *.o *.mod *.x

