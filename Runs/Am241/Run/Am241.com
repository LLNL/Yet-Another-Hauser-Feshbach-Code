time YAHFC.x << input
file Am241
all_discrete_states y
target Am241
projectile n
max_particle p 0
do_dwba y
lev_fit_d0 y
fit_Gamma_gamma y
scale_elastic   1.0 0.07    0.5
Delta_e        0.100
#proj_e_file   elist.input
proj_eminmax   0.1  20.1  0.2 
#
#---   No fission parameters entered here
#---   Relying on default parameters found in
#---   $YAHFC_DATA/Fission-Parameters.txt
#
WF_model  0 
preeq_model 1
preeq_M2_C1  1.00
preeq_M2_C2  1.00
preeq_M2_C3  1.20
track_gammas 0
num_mc_samp  100000
end
input
