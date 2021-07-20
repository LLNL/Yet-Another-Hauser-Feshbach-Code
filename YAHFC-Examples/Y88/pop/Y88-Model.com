time YAHFC.x << input
file Y88-Model
all_discrete_states 1
target Y88
#projectile -1  0   input-Model.pop
projectile -1  0   Test.pop
delta_e        0.05
global_lev_fit_d0 y
global_lev_option 1
lev_shell   Y88  -0.8005 0   0.05     -0.4872866
el_param    Y88 1 1  16.790  3.960    205.3
el_param    Y88 1 2   5.000  5.000      0.5
#
el_param    Y87 1 1  16.790  3.960    205.3
el_param    Y87 1 2   5.000  5.000      0.5
preeq_model 0 
track_gammas y
#----   "xs_only y" without angular distributions 
#----   (default for pop calcs, use"xs_only n" to overide
#xs_only n
num_mc_samp  100000
end
input
