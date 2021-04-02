#
#----   Sample run for n + U238
#----   Commands may appear in any order and any case
#----   lines starting with # or ! are comments
#
time YAHFC.x << input
file U238                       !   output file
all_discrete_states y                    !   use all discrete states with decay path to g.s. or isomer
#
#-----   target can be specified with Z, A or Isotope
#
target U238                       !   Target:  Z_target  A_target
#target 92  238                   !   Target:  Z_target  A_target
#
#-----   projectile can be specified with Z, A or label g,n,p,d,t,h,a
#
projectile n                    !   projectile:  Z_proj  A_proj
max_particle p 0                !   max # of particle allowed for type - no protons
ran_seed 29383615               !  If you want to use the same seed every time.
#----   Option to use MC sampling of particle decay with HF weights.
#----   Default is equal probablity for n,p,d,t,h,a, and fission
#biased_sampling y
#
#----   Do DWBA calculation with FRESCO
#----   Coupled Channels states are in
#----   $YAHFC_DATA/Coupled-Channels.txt
#----   This can be manually overriden
#
do_dwba y
#
#----   Rescale elastic scattering cross section
#----   Here energy dependent scaling of 7%
#----   For actinides, Soukhoviskii is default
#
scale_elastic   1.0 0.07    0.5
#
#----   Size of continuous energy bins
#
delta_e        0.10
#
#----   Two ways to make energy list
#----   proj_e_file which can have emax < delta_e
#----   or proj_eminmax  emin, emax, de but emin >= delta_e
#
#proj_e_file   elist-e-0.1.input
proj_eminmax   0.1  20.1  0.1
#
lev_fit_d0 y                    !  fit level density a-param to D0
pair_model 1                    !  use pairing model #1
lev_option 2                    !  use level density model #2 (default for actinides)
lev_ecut 92  238   1.15         !  manually set ecut for U238
#
#----   Fission barrier parameters. Defaults are found in
#----   $YAHFC_DATA/Fission-Parameters.txt
#
#
#----   Barriers for U23
#
f_num_barrier       U239  2
#---  U239 barrier #1
f_barrier           U239  1   6.05   0.3    !  height, hbw
f_lev_delta         U239  1   0.82          !  pairing
f_lev_shell         U239  1   0.8           !  shell correction
f_lev_gamma         U239  1   0.0760028     !  Igantiuk damping (default)
f_lev_spin          U239  1   0.0138900     !  spin cutoff (defaut)
f_lev_ematch        U239  1   3.75          !  matching energy
f_barrier_symmetry  U239  1   3             !  type of level density above barrier
#---  U239 barrier #2
f_barrier           U239  2   5.5   0.3
f_lev_delta         U239  2   0.0
f_lev_shell         U239  2   1.8
f_lev_gamma         U239  2   0.0760028
f_lev_ematch        U239  2   4.3
f_barrier_symmetry  U239  2   2
#
#----   Barriers for U238
#
f_num_barrier       U238  2
f_barrier           U238  1   6.4  0.7
f_lev_delta         U238  1   1.0
f_lev_shell         U238  1   1.3
f_lev_gamma         U238  1   0.0760028
f_lev_spin          U238  1   0.0138900
f_lev_ematch        U238  1   3.65
f_barrier_symmetry  U238  1   3
f_barrier           U238  2   5.8   0.6
f_lev_delta         U238  2   0.4
f_lev_shell         U238  2   1.6
f_lev_gamma         U238  2   0.0760028
f_lev_spin          U238  2   0.0138900
f_lev_ematch        U238  2   3.2
f_barrier_symmetry  U238  2   2
#
#----   Barriers for U237
#
f_num_barrier       U237  2
f_barrier           U237  1   5.8   0.5
f_lev_delta         U237  1   0.9
f_lev_shell         U237  1   0.0
f_lev_gamma         U237  1   0.0760028
f_lev_spin          U237  1   0.0138900
f_lev_ematch        U237  1   3.75
f_barrier_symmetry  U237  1   3
f_barrier           U237  2   5.3   0.6
f_lev_delta         U237  2   0.5
f_lev_shell         U237  2   0.50
f_lev_gamma         U237  2   0.0760028
f_lev_spin          U237  2   0.0138900
f_lev_ematch        U237  2   3.4
f_barrier_symmetry  U237  2   2
#
#----   Barriers for U236
#
F_Num_Barrier       U236  2
F_Barrier           U236  1   6.41  0.35
F_lev_delta         U236  1   0.0
F_lev_shell         U236  1   1.1   
F_lev_gamma         U236  1   0.0760028 
F_lev_spin          U236  1   0.0138900
F_lev_ematch        U236  1   0.7 0.75
f_barrier_symmetry  U236  1   3
F_Barrier           U236  2   6.10  0.4
F_lev_delta         U236  2   0.25
F_lev_shell         U236  2   0.5
F_lev_gamma         U236  2   0.0760028 
F_lev_spin          U236  2   0.0138900
F_lev_ematch        U236  2   2.0
f_barrier_symmetry  U236  2   1
#
#----    Width Fluctuation model 0 = none, 1 = Moldauer
#
WF_model  0
#
#----   Prequilibrium model input and parameters
# 
preeq_model 1
preeq_M2_C1  1.00
preeq_M2_C2  1.00
preeq_M2_C3  1.20
#
#----   Track all gamma decays 0/n/f = no  1/y/t = yes 
#
track_gammas n
#
#----   Modify gamma strength function to fit to Gamma_gamma
#
fit_Gamma_gamma y
#
#----   Cross section only (y) or full library (n) (default = n)
#
xs_only y
#
#----   Option to print to unit=6 while running
#----   when fitting, might be best to set to n
#
verbose_output y
#
#----   Extra angles helps with angular distributions- default = 10
#
#num_theta_angles 100
#
#----   Number of Monte Carlo Samples
#
num_mc_samp  100000
#
#---   All finished
#
end
input
