time YAHFC.x << input
file zr90 
target 40 90 
projectile p
use_unequal_bins y
all_discrete_states y 
delta_e        0.25
proj_eminmax   2.0 50.1  0.5
wf_model 0
preeq_model 1 
num_mc_samp  100000
end
input
