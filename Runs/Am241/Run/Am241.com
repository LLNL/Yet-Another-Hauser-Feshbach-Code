time YAHFC.x << input
file am241
All_gammas 1
target 95 241
projectile 0 1
max_particle 2 0
do_dwba t
scale_elastic   1.0 0.07    0.5
Delta_e        0.100
#proj_e_file   elist-e-0.1.input
proj_eminmax   0.1  20.1  0.2 
fit_aparam 1
lev_shell   95 241  0.25
lev_shell   95 240  0.25
lev_shell   95 239  0.25
lev_shell   95 238  0.25
#
F_Num_Barrier 95 242  2
F_Barrier     95 242  1   6.47   0.5    6.45   0.55
F_lev_delta   95 242  1   0.07
F_lev_shell   95 242  1   2.12
F_lev_gamma   95 242  1   0.0760028
F_lev_spin    95 242  1   0.0138900
F_lev_ematch  95 242  1   3.4
f_barrier_symmetry  95 242  1   3
F_Barrier     95 242  2   5.9   0.6  5.6   0.8
F_lev_delta   95 242  2   0.0
F_lev_shell   95 242  2   3.05
F_lev_gamma   95 242  2   0.0760028
F_lev_spin    95 242  1   0.0138900
F_lev_ematch  95 242  2   3.15
f_barrier_symmetry  95 242  2   2
#
F_Num_Barrier 95 241  2
F_Barrier     95 241  1   6.59  0.6
F_lev_delta   95 241  1   0.0
F_lev_shell   95 241  1   2.5
F_lev_gamma   95 241  1   0.0760028
F_lev_spin    95 241  1   0.0138900
F_lev_ematch  95 241  1   3.8
f_barrier_symmetry  95 241  1   3
F_Barrier     95 241  2   6.2   0.7
F_lev_delta   95 241  2   0.
F_lev_shell   95 241  2   3.2
F_lev_gamma   95 241  2   0.0760028
F_lev_spin    95 241  2   0.0138900
F_lev_ematch  95 241  2   3.6
f_barrier_symmetry  95 241  2   2
#
F_Num_Barrier 95 240  2
F_Barrier     95 240  1   6.22  0.5 6.2   0.6   5.5    0.8
F_lev_delta   95 240  1   0.0
F_lev_shell   95 240  1   3.6  2.8
F_lev_gamma   95 240  1   0.0760028
F_lev_spin    95 240  1   0.0138900
F_lev_ematch  95 240  1   3.8  3.6
f_barrier_symmetry  95 240  1   3
F_Barrier     95 240  2   5.7  0.6  5.8   0.7   5.0    1.0
F_lev_delta   95 240  2   0.0
F_lev_shell   95 240  2   3.5
F_lev_gamma   95 240  2   0.0760028
F_lev_spin    95 240  2   0.0138900
F_lev_ematch  95 240  2   3.6
f_barrier_symmetry  95 240  2   2
#
F_Num_Barrier 95 239  2
F_Barrier     95 239  1   6.2  0.6   5.8    0.6
F_lev_delta   95 239  1   0.0
F_lev_shell   95 239  1   2.5
F_lev_gamma   95 239  1   0.0760028
F_lev_spin    95 239  1   0.0138900
F_lev_ematch  95 239  1   3.4
f_barrier_symmetry  95 239  1   3
F_Barrier     95 239  2   5.8  0.6  5.5    0.6
F_lev_delta   95 239  2   0.0
F_lev_shell   95 239  2   2.9
F_lev_gamma   95 239  2   0.0760028
F_lev_spin    95 239  2   0.0138900
F_lev_ematch  95 239  2   3.60
f_barrier_symmetry  95 239  2   2
#
WF_model  0 
preeq_M2_C1  0.7
preeq_M2_C2  0.5
preeq_M2_C3 1.00
preeq_pair_model 2 -1.50
Preeq_V_n 7.0
track_gammas 0
fit_Gamma_gamma 1
track_primary_gammas 0
num_mc_samp  100000
end
input
