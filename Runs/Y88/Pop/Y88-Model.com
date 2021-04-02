time YAHFC.x << input
file Y88-Model
all_discrete_states 1
target Y88
projectile -1  0   input-Model.pop
Delta_e        0.05
lev_fit_d0 y
lev_shell   39 88  -0.8005 0   0.05     -0.4872866
lev_option 1
WF_model  1
e1_param  Y88 1  16.790  3.960    205.3
e1_param  Y88 2   5.000  5.000      0.5
e1_param  Y87 1  16.790  3.960    205.3
e1_param  Y87 2   5.000  5.000      0.5
preeq_model 0 
track_gammas y
#----   "xs_only y" without angular distributions 
#----   (default for pop calcs, use"xs_only n" to overide
#xs_only n
num_mc_samp  1000000
end
input
