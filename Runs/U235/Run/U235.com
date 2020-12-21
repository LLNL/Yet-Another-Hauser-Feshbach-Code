time YAHFC.x << input
file u235
All_gammas 1
target 92 235
projectile 0 1
max_particle 2 0
do_dwba t
scale_elastic   1.0 0.07    0.5
Delta_e        0.100
#proj_e_file   elist-e-0.1.input
proj_eminmax   0.1  20.1  0.2 
fit_aparam 1
pair_model 1
lev_option 2
#
lev_delta   92 235    0.0    1.700
lev_delta   92 234    0.6    1.5627
#
F_Num_Barrier       92 236  2
F_Barrier           92 236  1   6.41  0.35
F_lev_delta         92 236  1   0.0
F_lev_shell         92 236  1   1.1
F_lev_gamma         92 236  1   0.0760028
F_lev_spin          92 236  1   0.0138900
F_lev_ematch        92 236  1   0.7 
f_barrier_symmetry  92 236  1   3
F_Barrier           92 236  2   6.10  0.4
F_lev_delta         92 236  2   0.25 
F_lev_shell         92 236  2   0.5
F_lev_gamma         92 236  2   0.0760028
F_lev_spin          92 236  2   0.0138900
F_lev_ematch        92 236  2   2.0
f_barrier_symmetry  92 236  2   1
#
F_Num_Barrier       92 235  2
F_Barrier           92 235  1   5.86  0.6  6.45
F_lev_delta         92 235  1   0.28 
F_lev_shell         92 235  1   1.1
F_lev_gamma         92 235  1   0.0760028
F_lev_spin          92 235  1   0.0138900
F_lev_ematch        92 235  1   2.5
f_barrier_symmetry  92 235  1   3
F_Barrier           92 235  2   5.5   0.4
F_lev_delta         92 235  2   0.1
F_lev_shell         92 235  2   1.0
F_lev_gamma         92 235  2   0.0760028 
F_lev_spin          92 235  2   0.0138900
F_lev_ematch        92 235  2   2.9
f_barrier_symmetry  92 235  2   1
#
F_Num_Barrier       92 234  2
F_Barrier           92 234  1   5.8   0.6
F_lev_delta         92 234  1   0.15
F_lev_shell         92 234  1   0.75
F_lev_gamma         92 234  1   0.0760028 
F_lev_spin          92 234  1   0.0138900
F_lev_ematch        92 234  1   3.0
f_barrier_symmetry  92 234  1   3
F_Barrier           92 234  2   5.5   0.6 
F_lev_delta         92 234  2   0.25
F_lev_shell         92 234  2   0.75
F_lev_gamma         92 234  2   0.0760028
F_lev_spin          92 234  2   0.0138900
F_lev_ematch        92 234  2   2.0
f_barrier_symmetry  92 234  2   1
#
F_Num_Barrier       92 233  2
F_Barrier           92 233  1   6.1  0.5
F_lev_delta         92 233  1   1.00 
F_lev_shell         92 233  1   0.50
F_lev_gamma         92 233  1   0.0760028 
F_lev_spin          92 233  1   0.0138900
F_lev_ematch        92 233  1   2.00
f_barrier_symmetry  92 233  1   3
F_Barrier           92 233  2   5.5   0.6
F_lev_delta         92 233  2   0.5
F_lev_shell         92 233  2   2.00
F_lev_gamma         92 233  2   0.0760028 
F_lev_spin          92 233  2   0.0138900
F_lev_ematch        92 233  2   3.00
f_barrier_symmetry  92 233  2   1
#
WF_model  0 
preeq_model 1
preeq_M2_C1  1.00
preeq_M2_C2  1.00
preeq_M2_C3  1.20
output_mode 0
track_gammas 0
fit_Gamma_gamma 1
track_primary_gammas 0
num_mc_samp  1000000
end
input
