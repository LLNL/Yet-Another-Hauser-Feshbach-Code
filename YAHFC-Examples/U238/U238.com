#
#----   Sample run for n + U238
#----   Commands may appear in any order and mixed case
#----   lines starting with # or ! are comments
#
time YAHFC.x << input
file u238                       !   output file
ran_seed 1920293
all_discrete_states y                    !   use all discrete states with decay path to g.s. or isomer
target U238                   !   Target:  Z_target  A_target
projectile n                  !   projectile:  Z_proj  A_proj
max_particle p 0                !   max # of particle allowed for type - no protons
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
#----   Note if you change this and have do_dwba t
#----   you need to rerun the DWBA calculation
#
delta_e        0.100
#
#----   Two ways to make energy list
#----   proj_e_file which can have emax < delta_e
#----   or proj_eminmax  emin, emax, de but emin >= delta_e
#----   Commands to specify incident energies for projectile
#----   If commented out, default energy grid from 
#----   0.001 - 20.0 MeV is used
#
#proj_e_file   elist-1.input
#proj_eminmax   0.1 6.1  0.2
#
global_pair_model 1                    !  use pairing model #1
global_lev_option 2                    !  use level density model #2 (default for actinides)
#
#----   Parameters defining nuclear properties. Defaults are stored in
#----   $YAHFC_DATA/Saved-Parameters.txt
#
#---  U239
lev_option     U239   2
e1_model       U239   1
lev_delta      U239   0.000
el_param       U239   1  1   11.00   2.62    299.0
el_param       U239   1  2   14.0    4.53    382.0
el_param       U239   1  3    5.0    5.0      -5.345
f_num_barrier       U239  2
f_use_tran_states   U239  n
f_barrier           U239  1   6.13  0.3 0.3    !  height, hbw
f_barrier_damp      U239  1  -1.0   80.0
f_lev_delta         U239  1   0.00          !  pairing
f_lev_shell         U239  1  -0.30          !  shell correction
f_lev_ematch        U239  1   1.78          !  matching energy
f_barrier_symmetry  U239  1   3             !  type of level density above barrier
f_barrier           U239  2   5.45   0.4
f_barrier_damp      U239  2  -1.0   80.0
f_lev_delta         U239  2   0.22
f_lev_shell         U239  2  -0.3
f_lev_ematch        U239  2   2.0
f_barrier_symmetry  U239  2   2
#---  U238
lev_option     U238   2
e1_model       U238   1
lev_delta      U238   0.778
lev_ecut       U238   1.15    1.15         !  manually set ecut for U238
el_param       U238   1  1   11.00   2.62    299.0   11.20  3.0   300.0  300.0
el_param       U238   1  2   14.0    4.53    382.0   14.40  4.8   330.0
el_param       U238   1  3    5.0    5.0      -4.7  18.8
f_num_barrier       U238  2
f_use_tran_states   U238  n
f_barrier           U238  1   6.15  0.65  6.10  0.7
f_barrier_damp      U238  1  -1.0   80.00
f_lev_delta         U238  1   0.55
f_lev_shell         U238  1  -1.0
f_lev_ematch        U238  1   1.9
f_barrier_symmetry  U238  1   3
f_barrier           U238  2   5.6  0.6 0.5  5.85   0.5
f_barrier_damp      U238  2  -1.0   80.00
f_lev_delta         U238  2   0.5
f_lev_shell         U238  2   0.8
f_lev_ematch        U238  2   1.7
f_barrier_symmetry  U238  2   2
#---  U237
lev_option     U237   2
e1_model       U237   1
el_param       U237   1  1   11.00   2.62    299.0
el_param       U237   1  2   14.0    4.53    382.0
el_param       U237   1  3    5.0    5.0      -6.1
lev_delta      U237   0.0         !  manually set ecut for U238
lev_ecut       U237   0.6         !  manually set ecut for U238
f_num_barrier       U237  2
f_barrier           U237  1   6.0  0.5    !  height, hbw
f_use_tran_states   U237  n
f_barrier_damp      U237  1   -1.0d0   80.00
f_lev_delta         U237  1   0.2          !  pairing
f_lev_shell         U237  1  -0.4         !  shell correction
f_lev_ematch        U237  1   2.0   3.75          !  matching energy
f_barrier_symmetry  U237  1   3             !  type of level density above barrier
f_barrier           U237  2   5.5   0.5
f_barrier_damp      U237  2   -1.0d0   80.00
f_lev_delta         U237  2   0.47
f_lev_shell         U237  2   1.0
f_lev_ematch        U237  2   2.0
f_barrier_symmetry  U237  2   2
#---  U236
lev_option     U236  2
e1_model       U236  1
lev_delta      U236  0.781
el_param       U236   1  1   11.00   2.62    299.0
el_param       U236   1  2   14.0    4.53    382.0
el_param       U236   1  3    5.0    5.0      -5.23
f_Num_Barrier       U236  2
f_Barrier           U236  1   5.55  0.7
f_use_tran_states   U236  n
f_barrier_damp      U236  1  -1.0   80.0
f_lev_delta         U236  1   0.52
f_lev_shell         U236  1   1.8
f_lev_ematch        U236  1   2.0
f_barrier_symmetry  U236  1   1
f_Barrier           U236  2   6.0  0.7 0.6
f_barrier_damp      U236  2  -1.0   80.0
f_lev_delta         U236  2   0.7
f_lev_shell         U236  2   0.5
f_lev_ematch        U236  2   2.15
f_barrier_symmetry  U236  2   2
#---  U235
lev_option     U235  2
e1_model       U235  1
lev_delta      U235   0.1
#lev_ecut       U235   0.6
el_param       U235   1  1   11.00   2.50    299.00 
el_param       U235   1  2   14.0    4.53    382.00
el_param       U235   1  3    5.0    5.0      -5.77
f_Barrier           U235  1   5.4  0.4 
f_use_tran_states   U235  n
f_barrier_damp      U235  1  -1.0   80.0
f_lev_delta         U235  1   0.1
f_lev_shell         U235  1   2.2
f_lev_ematch        U235  1   2.6
f_barrier_symmetry  U235  1   1
f_Barrier           U235  2   5.875   0.5
f_barrier_damp      U235  2  -1.0   80.0
f_lev_delta         U235  2   0.25
f_lev_shell         U235  2   2.0
f_lev_ematch        U235  2   2.8
f_barrier_symmetry  U235  2   2
#
#----    Width Fluctuation model 0 = none, 1 = Moldauer
#
wf_model  1
#
#----    Calculatue only cross section (y) default= n
#
xs_only n
#
#----   Prequilibirum mode data
# 
preeq_model 1
#
#----   Track all gamma decays 0/n/f = no  1/y/t = yes 
#
track_gammas n
#
#----   Number of Monte Carlo Samples
#
num_mc_samp  1000000
#
#---   All finished
#
end
input
