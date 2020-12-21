time YAHFC.x << input
file y87 
target 39 87 
projectile 0 1
delta_e        0.2
proj_eminmax   0.1 20.3  0.2
#proj_e_file   elist.input 
fit_aparam 0
lev_option 1
#lev_shell   39 88  -0.8005 0   0.05     -0.4872866
WF_model  1
E1_param  39 88 1  16.790  3.960    205.3
E1_param  39 88 2   5.000  5.000      0.5
E1_param  39 87 1  16.790  3.960    205.3
E1_param  39 87 2   5.000  5.000      0.5
preeq_model 1 
track_gammas 0
All_gammas 1 
Out_Gammas_vs_E 0  
output_mode 1
num_mc_samp  100000
end
input
