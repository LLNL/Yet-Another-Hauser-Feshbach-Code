time YAHFC.x << input
file zr90 
target Zr90 
projectile n
all_discrete_states y  
delta_e        0.1
proj_eminmax   0.1  10.1  0.2
wf_model 0
preeq_model 1 
num_mc_samp  100000
end
input
