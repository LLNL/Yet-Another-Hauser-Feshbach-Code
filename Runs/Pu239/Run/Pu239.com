time YAHFC.x << input
file Pu239
All_gammas 1
target 94 239
projectile 0 1
max_particle 2 0
do_dwba t
scale_elastic   1.0 0.07    0.5
Delta_e        0.100
proj_e_file   elist-e-0.1.input
#proj_eminmax   0.1  20.1  0.2 
fit_aparam 1
lev_option 2
lev_delta    94 240   0.670
lev_delta    94 239   0.250 0.310
lev_delta    94 238   0.694
lev_delta    94 237   0.310
lev_delta    94 236   0.681
lev_shell    94 238   0.050
lev_shell    94 237   0.400
lev_shell    94 236   0.050
lev_ecut     94 240   0.750
lev_ecut     94 239   0.749
lev_ecut     94 239   0.800
lev_ecut     94 238   0.800
lev_ecut     94 237   0.330
lev_ecut     94 236   0.600
lev_parity_fac      94 239 -1.0 0.39  1.0
F_Barrier           94 240  1   6.5  0.45 
F_lev_delta         94 240  1   0.05 
F_lev_shell         94 240  1   1.7 
F_lev_ematch        94 240  1   0.85
F_Barrier           94 240  2   6.25  0.6  
F_lev_delta         94 240  2   0.45 
F_lev_shell         94 240  2   1.4 
F_lev_ematch        94 240  2   3.0
f_num_barrier       94 239  2
f_barrier           94 239  1   5.85  0.45
f_lev_delta         94 239  1   0.3 
f_lev_shell         94 239  1   1.2 
f_lev_gamma         94 239  1   0.0760028
f_lev_spin          94 239  1   0.0138900 
f_lev_ematch        94 239  1   2.4 2.5
f_barrier_symmetry  94 239  1   3
f_barrier           94 239  2   5.55    0.5
f_lev_delta         94 239  2   0.4 
f_lev_shell         94 239  2   1.0
f_lev_gamma         94 239  2   0.0760028
f_lev_spin          94 239  2   0.0138900 
f_lev_ematch        94 239  2   2.0
f_barrier_symmetry  94 239  2   2
WF_model  0 
preeq_model 1
preeq_M2_C1  1.00
output_mode 0
track_gammas 0
fit_Gamma_gamma 1
track_primary_gammas 0
num_mc_samp  1000000
end
input
