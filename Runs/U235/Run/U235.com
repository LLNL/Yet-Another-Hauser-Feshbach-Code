time YAHFC.x << input
file u235
all_discrete_states y
target U235
projectile n
max_particle p 0
do_dwba y
scale_elastic   1.0 0.07    0.5
Delta_e        0.100
#proj_e_file   elist-e-0.1.input
proj_eminmax   0.1  20.1  0.2 
lev_fit_d0 y
fit_gamma_gamma y
pair_model 1
lev_option 2
#
lev_delta   U235    0.0    1.700
lev_delta   U234    0.6    1.5627
#
F_Num_Barrier       U236  2
F_Barrier           U236  1   6.41  0.35
F_lev_delta         U236  1   0.0
F_lev_shell         U236  1   1.1
F_lev_gamma         U236  1   0.0760028
F_lev_spin          U236  1   0.0138900
F_lev_ematch        U236  1   0.7 
f_barrier_symmetry  U236  1   3
F_Barrier           U236  2   6.10  0.4
F_lev_delta         U236  2   0.25 
F_lev_shell         U236  2   0.5
F_lev_gamma         U236  2   0.0760028
F_lev_spin          U236  2   0.0138900
F_lev_ematch        U236  2   2.0
f_barrier_symmetry  U236  2   1
#
F_Num_Barrier       U235  2
F_Barrier           U235  1   5.86  0.6  6.45
F_lev_delta         U235  1   0.28 
F_lev_shell         U235  1   1.1
F_lev_gamma         U235  1   0.0760028
F_lev_spin          U235  1   0.0138900
F_lev_ematch        U235  1   2.5
f_barrier_symmetry  U235  1   3
F_Barrier           U235  2   5.5   0.4
F_lev_delta         U235  2   0.1
F_lev_shell         U235  2   1.0
F_lev_gamma         U235  2   0.0760028 
F_lev_spin          U235  2   0.0138900
F_lev_ematch        U235  2   2.9
f_barrier_symmetry  U235  2   1
#
F_Num_Barrier       U234  2
F_Barrier           U234  1   5.8   0.6
F_lev_delta         U234  1   0.15
F_lev_shell         U234  1   0.75
F_lev_gamma         U234  1   0.0760028 
F_lev_spin          U234  1   0.0138900
F_lev_ematch        U234  1   3.0
f_barrier_symmetry  U234  1   3
F_Barrier           U234  2   5.5   0.6 
F_lev_delta         U234  2   0.25
F_lev_shell         U234  2   0.75
F_lev_gamma         U234  2   0.0760028
F_lev_spin          U234  2   0.0138900
F_lev_ematch        U234  2   2.0
f_barrier_symmetry  U234  2   1
#
F_Num_Barrier       U233  2
F_Barrier           U233  1   6.1  0.5
F_lev_delta         U233  1   1.00 
F_lev_shell         U233  1   0.50
F_lev_gamma         U233  1   0.0760028 
F_lev_spin          U233  1   0.0138900
F_lev_ematch        U233  1   2.00
f_barrier_symmetry  U233  1   3
F_Barrier           U233  2   5.5   0.6
F_lev_delta         U233  2   0.5
F_lev_shell         U233  2   2.00
F_lev_gamma         U233  2   0.0760028 
F_lev_spin          U233  2   0.0138900
F_lev_ematch        U233  2   3.00
f_barrier_symmetry  U233  2   1
#
WF_model  0 
preeq_model 1
preeq_M2_C1  1.00
preeq_M2_C2  1.00
preeq_M2_C3  1.20
track_gammas n
num_mc_samp  100000
end
input
