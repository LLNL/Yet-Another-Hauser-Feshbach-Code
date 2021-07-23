ltime YAHFC.x << input
file Y87 
all_discrete_states y
target Y87 
projectile n
delta_e        0.05
proj_eminmax   0.1 10.0  0.2
#proj_e_file   elist.input 
global_lev_fit_d0 y
global_lev_option 1
el_param  Y88 1 1  16.790  3.960    205.3
el_param  Y88 1 2   5.000  5.000      0.5
el_param  Y87 1 1  16.790  3.960    205.3
el_param  Y87 1 2   5.000  5.000      0.5
wf_model 1
preeq_model 1 
track_gammas 0
num_mc_samp  1000000
end
input
