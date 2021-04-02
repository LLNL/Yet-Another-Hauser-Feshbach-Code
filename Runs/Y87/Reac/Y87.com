time YAHFC.x << input
file Y87 
all_discrete_states y
target Y87 
projectile n
delta_e        0.05
proj_eminmax   0.1 10.0  0.2
#proj_e_file   elist.input 
lev_fit_d0 y
lev_option 1
E1_param  Y88 1  16.790  3.960    205.3
E1_param  Y88 2   5.000  5.000      0.5
E1_param  Y87 1  16.790  3.960    205.3
E1_param  Y87 2   5.000  5.000      0.5
WF_model 1
preeq_model 1 
track_gammas 0
num_mc_samp  100000
end
input
