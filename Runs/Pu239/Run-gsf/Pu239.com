time YAHFC.x << input
file Pu239
all_discrete_states y
target Pu239
projectile n
max_particle p 0
do_dwba y
track_gammas n
fit_Gamma_gamma y
scale_elastic   1.0 0.07    0.5
Delta_e        0.100
#proj_e_file   elist-1.input
proj_eminmax   0.1  20.1  0.2 
lev_fit_d0 y
lev_fit_aparam y
lev_option 2
lev_delta    Pu240   0.670
lev_delta    Pu239   0.250 0.310
lev_delta    Pu238   0.694
lev_delta    Pu237   0.310
lev_delta    Pu236   0.681
lev_shell    Pu238   0.050
lev_shell    Pu237   0.400
lev_shell    Pu236   0.050
lev_ecut     Pu240   0.750
lev_ecut     Pu239   0.749
lev_ecut     Pu239   0.800
lev_ecut     Pu238   0.800
lev_ecut     Pu237   0.330
lev_ecut     Pu236   0.600
lev_parity_fac      Pu239 -1.0 0.39  1.0
#
f_num_barrier       Pu240  2
f_barrier           Pu240  1   6.52  0.45  
f_lev_delta         Pu240  1   0.06
f_lev_shell         Pu240  1   1.7
f_lev_gamma         Pu240  1   0.0760028
f_lev_spin          Pu240  1   0.0138900
f_lev_ematch        Pu240  1   0.75
f_barrier_symmetry  Pu240  1   3
f_barrier           Pu240  2   6.15  0.6
f_lev_delta         Pu240  2   0.5
f_lev_shell         Pu240  2   1.5
f_lev_gamma         Pu240  2   0.0760028
f_lev_spin          Pu240  2   0.0138900 
f_lev_ematch        Pu240  2   3.0
f_barrier_symmetry  Pu240  2   2
#
f_num_barrier       Pu239  2
f_barrier           Pu239  1   5.85  0.45
f_lev_delta         Pu239  1   0.3 
f_lev_shell         Pu239  1   1.2 
f_lev_gamma         Pu239  1   0.0760028
f_lev_spin          Pu239  1   0.0138900 
f_lev_ematch        Pu239  1   2.4 2.5
f_barrier_symmetry  Pu239  1   3
f_barrier           Pu239  2   5.55    0.5
f_lev_delta         Pu239  2   0.4 
f_lev_shell         Pu239  2   1.3
f_lev_gamma         Pu239  2   0.0760028
f_lev_spin          Pu239  2   0.0138900 
f_lev_ematch        Pu239  2   2.0
f_barrier_symmetry  Pu239  2   2
#
f_num_barrier       Pu238  2
f_barrier           Pu238  1   5.5   0.6
f_lev_delta         Pu238  1   0.0
f_lev_shell         Pu238  1   0.50
f_lev_gamma         Pu238  1   0.0760028
f_lev_spin          Pu238  1   0.0138900
f_lev_ematch        Pu238  1   2.30
f_barrier_symmetry  Pu238  1   1
f_barrier           Pu238  2   5.20   0.6
f_lev_delta         Pu238  2   0.3
f_lev_shell         Pu238  2   0.0
f_lev_gamma         Pu238  2   0.0760028
f_lev_spin          Pu238  2   0.0138900
f_lev_ematch        Pu238  2   1.00
f_barrier_symmetry  Pu238  2   2
#
start  94  237   2
f_num_barrier       Pu237  2
f_barrier           Pu237  1   5.2  0.6
f_lev_delta         Pu237  1   0.3
f_lev_shell         Pu237  1   1.5
f_lev_gamma         Pu237  1   0.0760028
f_lev_spin          Pu237  1   0.0138900
f_lev_ematch        Pu237  1   2.6
f_barrier_symmetry  Pu237  1   1
f_barrier           Pu237  2   5.0  0.6
f_lev_delta         Pu237  2   0.3
f_lev_shell         Pu237  2   1.5
f_lev_gamma         Pu237  2   0.0760028
f_lev_spin          Pu237  2   0.0138900
f_lev_ematch        Pu237  2   2.0
f_barrier_symmetry  Pu237  2   2
#
f_num_barrier       Pu236  2
f_barrier           Pu236  1   5.55    0.8
f_lev_delta         Pu236  1   0.4
f_lev_shell         Pu236  1   0.0
f_lev_gamma         Pu236  1   0.0760028
f_lev_spin          Pu236  1   0.0138900
f_lev_ematch        Pu236  1   2.6
f_barrier_symmetry  Pu236  1   1
f_barrier           Pu236  2   5.15    1.0
f_lev_delta         Pu236  2   0.50
f_lev_shell         Pu236  2   0.00
f_lev_gamma         Pu236  2   0.0760028
f_lev_spin          Pu236  2   0.0138900
f_lev_ematch        Pu236  2   2.0
f_barrier_symmetry  Pu236  2   2
#
f_num_barrier       Pu235  2
f_barrier           Pu235  1   5.50    0.8
f_lev_delta         Pu235  1   0.4
f_lev_shell         Pu235  1   0.0
f_lev_gamma         Pu235  1   0.0760028
f_lev_spin          Pu235  1   0.0138900
f_lev_ematch        Pu235  1   2.6
f_barrier_symmetry  Pu235  1   1
f_barrier           Pu235  2   5.15    1.0
f_lev_delta         Pu235  2   0.50
f_lev_shell         Pu235  2   0.00
f_lev_gamma         Pu235  2   0.0760028
f_lev_spin          Pu235  2   0.0138900
f_lev_ematch        Pu235  2   2.0
f_barrier_symmetry  Pu235  2   2
#
#e1_param Pu240 3 7.0 3.0  2.0
read_el_gsf Pu240  1  f  Pu240-E1-1.dat
read_el_gsf Pu240  1  f  Pu240-E1-2.dat
WF_model  1 
preeq_model 1
preeq_M2_C1  1.00
preeq_M2_C2  1.00
preeq_M2_C3  1.20
num_mc_samp  100000
end
input
