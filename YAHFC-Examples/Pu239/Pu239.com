time YAHFC.x << input
file Pu239
all_discrete_states y
target Pu239
projectile n
max_particle p 0
do_dwba y
track_gammas n
scale_elastic   1.0 0.07    0.5
Delta_e        0.100
#proj_e_file   elist-1.input
proj_eminmax  10.1  20.0  0.2 
#
#######     Pu240   #####################
lev_option   Pu240   2
lev_delta    Pu240   0.670
#lev_ecut     Pu240   0.750
e1_model     Pu240   1
el_param     Pu240   1 1  11.28 3.25  325.0
el_param     Pu240   1 2  13.73 4.25  384.0
el_param     Pu240   1 3    5.0 5.0    -7.075
f_num_barrier       Pu240  2
f_use_tran_states   Pu240  n
f_barrier           Pu240  1   6.5  0.4  
f_barrier_damp      Pu240  1  -1.0   80.0
f_lev_delta         Pu240  1   0.1
f_lev_shell         Pu240  1   1.7
f_lev_ematch        Pu240  1   0.8
f_barrier_symmetry  Pu240  1   3
f_barrier           Pu240  2   6.3  0.4
f_barrier_damp      Pu240  2  -1.0   80.0
f_lev_delta         Pu240  2   0.3
f_lev_shell         Pu240  2   1.2
f_lev_ematch        Pu240  2   3.0
f_barrier_symmetry  Pu240  2   2
#######     Pu239   #####################
lev_option   Pu239   2
lev_delta    Pu239   0.250 0.310
#lev_ecut     Pu239   0.749
e1_model     Pu239   1
el_param     Pu239   1 1  11.28 3.25  325.0
el_param     Pu239   1 2  13.73 4.25  384.0
el_param     Pu239   1 3   5.0 5.0    -7.02
lev_parity_fac      Pu239 -1.0 0.39  1.0
f_num_barrier       Pu239  2
f_use_tran_states   Pu239  n
f_barrier           Pu239  1   5.85  0.45
f_barrier_damp      Pu239  1  -1.0   80.0
f_lev_delta         Pu239  1   0.33 
f_lev_shell         Pu239  1   1.2 
f_lev_ematch        Pu239  1   2.4 2.5
f_barrier_symmetry  Pu239  1   3
f_barrier           Pu239  2   5.6    0.5
f_barrier_damp      Pu239  2  -1.0   80.0
f_lev_delta         Pu239  2   0.4 
f_lev_shell         Pu239  2   1.5
f_lev_ematch        Pu239  2   2.0
f_barrier_symmetry  Pu239  2   2
#######     Pu240   #####################
lev_option   Pu238   2
lev_delta    Pu238   0.694
lev_shell    Pu238   0.050
e1_model     Pu238   1
el_param     Pu238   1 1  11.28 3.25  325.0
el_param     Pu238   1 2  13.73 4.25  384.0
el_param     Pu238   1 3    5.0 5.0    -7.02
f_num_barrier       Pu238  2
f_use_tran_states   Pu238  n
f_barrier           Pu238  1   5.5   0.6
f_barrier_damp      Pu238  1  -1.0   80.0
f_lev_delta         Pu238  1   0.0
f_lev_shell         Pu238  1   0.50
f_lev_ematch        Pu238  1   2.30
f_barrier_symmetry  Pu238  1   1
f_barrier           Pu238  2   5.20   0.6
f_barrier_damp      Pu238  2  -1.0   80.0
f_lev_delta         Pu238  2   0.3
f_lev_shell         Pu238  2   0.0
f_lev_ematch        Pu238  2   1.00
f_barrier_symmetry  Pu238  2   2
#######     Pu237   #####################
lev_option   Pu237   2
lev_delta    Pu237   0.310
lev_shell    Pu237   0.400
lev_ecut     Pu237   0.330
e1_model     Pu237   1
el_param     Pu237   1 1  11.28 3.25  325.0
el_param     Pu237   1 2  13.73 4.25  384.0
el_param     Pu237   1 3    5.0 5.0    -7.02
f_num_barrier       Pu237  2
f_use_tran_states   Pu237  n
f_barrier           Pu237  1   5.2  0.6
f_barrier_damp      Pu237  1  -1.0   80.0
f_lev_delta         Pu237  1   0.3
f_lev_shell         Pu237  1   1.5
f_lev_ematch        Pu237  1   2.6
f_barrier_symmetry  Pu237  1   1
f_barrier           Pu237  2   5.0  0.6
f_barrier_damp      Pu237  2  -1.0   80.0
f_lev_delta         Pu237  2   0.3
f_lev_shell         Pu237  2   1.5
f_lev_ematch        Pu237  2   2.0
f_barrier_symmetry  Pu237  2   2
#######     Pu236   #####################
lev_option   Pu236   2
lev_delta    Pu236   0.681
lev_shell    Pu236   0.050
lev_ecut     Pu236   0.600
e1_model     Pu236   1
el_param     Pu236   1 1  11.28 3.25  325.0
el_param     Pu236   1 2  13.73 4.25  384.0
el_param     Pu236   1 3    5.0 5.0    -7.02
f_num_barrier       Pu236  2
f_use_tran_states   Pu236  n
f_barrier           Pu236  1   5.55    0.8
f_barrier_damp      Pu236  1  -1.0   80.0
f_lev_delta         Pu236  1   0.4
f_lev_shell         Pu236  1   0.0
f_lev_ematch        Pu236  1   2.6
f_barrier_symmetry  Pu236  1   1
f_barrier           Pu236  2   5.15    1.0
f_barrier_damp      Pu236  2  -1.0   80.0
f_lev_delta         Pu236  2   0.50
f_lev_shell         Pu236  2   0.00
f_lev_ematch        Pu236  2   2.0
f_barrier_symmetry  Pu236  2   2
#######     Pu235   #####################
lev_option   Pu235   2
lev_delta    Pu235   0.670
e1_model     Pu235   1
el_param     Pu235   1 1  11.28 3.25  325.0
el_param     Pu235   1 2  13.73 4.25  384.0
el_param     Pu235   1 3    5.0 5.0    -7.02
f_num_barrier       Pu235  2
f_use_tran_states   Pu235  n
f_barrier           Pu235  1   5.50    0.8
f_barrier_damp      Pu235  1  -1.0   80.0
f_lev_delta         Pu235  1   0.4
f_lev_shell         Pu235  1   0.0
f_lev_ematch        Pu235  1   2.6
f_barrier_symmetry  Pu235  1   1
f_barrier           Pu235  2   5.15    1.0
f_barrier_damp      Pu235  2  -1.0   80.0
f_lev_delta         Pu235  2   0.50
f_lev_shell         Pu235  2   0.00
f_lev_ematch        Pu235  2   2.0
f_barrier_symmetry  Pu235  2   2
#
xs_only y
WF_model  0 
preeq_model 1
num_mc_samp  100000
end
input
