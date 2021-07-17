time YAHFC.x << input
file U235
all_discrete_states y
target U235
projectile n
max_particle p 0
do_dwba y
scale_elastic   1.0 0.04    0.5
Delta_e        0.100
#proj_e_file   elist-1.input
proj_eminmax  5.1 14.1  0.2 
#---  U236
lev_option     U236  2
e1_model       U236  1
lev_delta      U236  0.781
el_param       U236   1  1   11.00   2.62    299.0
el_param       U236   1  2   14.0    4.53    382.0
el_param       U236   1  3    5.0    5.0      -5.23
f_num_Barrier       U236  2
f_use_tran_states   U236  n
f_Barrier           U236  1   5.6  0.7
f_barrier_damp      U236  1  -1.0   80.0
f_lev_delta         U236  1   0.4
f_lev_shell         U236  1   1.8
f_lev_ematch        U236  1   2.0
f_barrier_symmetry  U236  1   1
f_Barrier           U236  2   6.1  0.5 0.6
f_barrier_damp      U236  2  -1.0   80.0
f_lev_delta         U236  2   0.6
f_lev_shell         U236  2   0.5
f_lev_ematch        U236  2   2.3
f_barrier_symmetry  U236  2   2
#---  U235
lev_option     U235  2
e1_model       U235  1
lev_delta      U235   0.1
el_param       U235   1  1   11.00   2.50    299.00 
el_param       U235   1  2   14.0    4.53    382.00
el_param       U235   1  3    5.0    5.0      -5.525
f_num_Barrier       U235  2
f_Barrier           U235  1   5.6    0.4 
f_use_tran_states   U235  n
f_barrier_damp      U235  1  -1.0   80.0
f_lev_delta         U235  1   0.23
f_lev_shell         U235  1   2.2
f_lev_ematch        U235  1   2.75
f_barrier_symmetry  U235  1   1
f_Barrier           U235  2   5.8   0.4
f_barrier_damp      U235  2  -1.0   80.0
f_lev_delta         U235  2   0.35
f_lev_shell         U235  2   2.2
f_lev_ematch        U235  2   3.15
f_barrier_symmetry  U235  2   2
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
f_lev_delta         U234  1   1.35
f_lev_shell         U234  1   2.2
f_lev_ematch        U234  1   2.7
f_barrier_symmetry  U234  1   1
f_Barrier           U234  2   5.25  0.4  
f_barrier_damp      U234  2  -1.0   80.0
f_lev_delta         U234  2   1.37  
f_lev_shell         U234  2   2.0
f_lev_ematch        U234  2   2.7
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
lev_option U232  2
lev_shell      U232  -0.45
lev_delta      U232   0.784 
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
#
xs_only y
Wf_model  0
track_gammas n
num_mc_samp  100000
end
input
