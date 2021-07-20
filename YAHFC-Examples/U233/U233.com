time YAHFC.x << input
file U233
all_discrete_states y
target U233
projectile n
max_particle p 0
do_dwba y
scale_elastic   1.0 0.07    0.5
Delta_e        0.100
#proj_e_file   elist-1.input
proj_eminmax   5.1  14.1  0.2 
global_e1_model 1
global_lev_fit_d0 y
global_lev_option 2
# ----   U234
lev_option     U234   2
lev_delta      U234   0.781
e1_model       U234   1
el_param       U234   1  1   11.00   2.62    299.00 
el_param       U234   1  2   14.0    4.53    382.00
el_param       U234   1  3    5.0    5.0      -5.868
f_num_barrier       U234  2
f_use_tran_states   U234  n 
f_Barrier           U234  1   4.9  0.4
f_barrier_damp      U234  1  -1.0   80.0
F_lev_delta         U234  1   1.35
F_lev_shell         U234  1   2.2
F_lev_ematch        U234  1   2.7
f_barrier_symmetry  U234  1   1
F_Barrier           U234  2   5.25  0.4  
f_barrier_damp      U234  2  -1.0   80.0
F_lev_delta         U234  2   1.37  
F_lev_shell         U234  2   2.0
F_lev_ematch        U234  2   2.7
f_barrier_symmetry  U234  2   2
# ----   U233
lev_option     U233   2
lev_delta      U233   0.0
e1_model       U233   1
el_param       U233   1  1   11.00   2.50    299.00 
el_param       U233   1  2   14.0    4.53    382.00
el_param       U233   1  3    5.0    5.0      -5.868
f_num_barrier       U233  2
f_use_tran_states   U233  n 
f_barrier           U233  1   5.1   0.4
f_barrier_damp      U233  1  -1.0   80.0
f_lev_delta         U233  1   0.5 
f_lev_shell         U233  1   1.8 
f_lev_ematch        U233  1   3.3
f_barrier_symmetry  U233  1   1
f_barrier           U233  2   5.65    0.3
f_barrier_damp      U233  2  -1.0   80.0
f_lev_delta         U233  2   0.5
f_lev_shell         U233  2   1.8
f_lev_ematch        U233  2   3.3
f_barrier_symmetry  U233  2   2
# ----   U232
lev_option     U232  2
lev_shell      U232   -0.45
lev_delta      U232    0.784 
e1_model       U232   1
el_param       U232   1  1   11.00   2.62    299.00 
el_param       U232   1  2   14.0    4.53    382.00
el_param       U232   1  3    5.0    5.0      -5.868
f_num_barrier       U232  2
f_use_tran_states   U232  n 
f_barrier           U232  1   5.0  0.45
f_barrier_damp      U232  1  -1.0   80.0
f_lev_delta         U232  1   0.3 
f_lev_shell         U232  1   1.2 
f_lev_ematch        U232  1   3.0
f_barrier_symmetry  U232  1   1
f_barrier           U232  2   5.5    0.5
f_barrier_damp      U232  2  -1.0   80.0
f_lev_delta         U232  2   0.3 
f_lev_shell         U232  2   1.2 
f_lev_ematch        U232  2   3.0
f_barrier_symmetry  U232  2   2
# ----   U231
lev_option     U231   2
lev_delta      U231   0.0    
lev_shell      U231   0.0
lev_ecut       U231   0.90
e1_model       U231   1
el_param       U231   1  1   11.00   2.62    299.00 
el_param       U231   1  2   14.0    4.53    382.00
el_param       U231   1  3    5.0    5.0      -5.868
f_num_barrier       U231  2
f_use_tran_states   U231  n 
f_barrier           U231  1   5.85  0.45
f_lev_delta         U231  1   0.3 
f_lev_shell         U231  1   1.2 
f_lev_ematch        U231  1   2.4 2.5
f_barrier_symmetry  U231  1   1
f_barrier           U231  2   5.55    0.5
f_lev_delta         U231  2   0.4 
f_lev_shell         U231  2   1.0
f_lev_ematch        U231  2   2.0
f_barrier_symmetry  U231  2   2
#
WF_model  0
preeq_model 1
track_gammas n
xs_only y
num_mc_samp  200000
end
input
