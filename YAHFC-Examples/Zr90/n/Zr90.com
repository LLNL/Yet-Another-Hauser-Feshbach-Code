time YAHFC.x << input
file Zr90 
target Zr90 
projectile n
all_discrete_states y  
delta_e        0.1
#----   Commands to specify incident energies for projectile
#----   If commented out, default energy grid from 
#----   0.0001 - 20.0 MeV is used
proj_eminmax   0.1  20.1  0.2
wf_model 1
preeq_model 1 
num_mc_samp  100000
end
input
