time YAHFC.x << input
file Y88-Model
All_gammas 1
target 39 88
projectile -1  0   input-Model.pop
Delta_e        0.05
fit_aparam 0
lev_shell   39 88  -0.8005 0   0.05     -0.4872866
lev_option 1
WF_model  1
E1_param  39 88 1  16.790  3.960    205.3
E1_param  39 88 2   5.000  5.000      0.5
E1_param  39 87 1  16.790  3.960    205.3
E1_param  39 87 2   5.000  5.000      0.5
preeq_model 0 
track_gammas 1
output_mode 2
num_mc_samp  1000000
end
input
