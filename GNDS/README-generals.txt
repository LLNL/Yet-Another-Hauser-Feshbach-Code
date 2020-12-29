Simple demonstrations and examples of extracting data from GNDS files in XML form
Ian Thompson
December 2019

Self-contained python scripts with LC fudgePath include last as default.

These scripts should be copied and extended by users, putting in desired loops over energies and angles.


gnds_reaction_list.py  <xml>
 - print list of all reactions


gnds_reaction_crosssection.py <xml> <reaction> <Ein>
 - cross section (barns) for <reaction> at incident energy <Ein> (MeV)


gnds_reaction_angulardistribution.py <xml> <reaction> <Ein> <product> 
 - angular distribution (barns/sr) for <reaction> at incident energy <Ein> (MeV) 
   at angles [10.,180.] (deg) in steps of 10 deg.


gnds_reaction_production.py: <xml> <reaction> <product> <Emax>
 - inclusive production cross-section,  for specific <reaction>, 
   for <product> at incident energies Ein from <Emax>/10 to <Emax>, 
   and exit energies in range [0,Ein],
   integrated over all angles.
     (options for muOut and phiOut ranges are not used).


gnds_inclusive_production.py <xml> <product> <Ein> <Eoutmin> <Eoutmax>
 - inclusive production cross-section,  summed over all reactions, 
   for <product> at <Ein> incident energy, 
   and exit energies in range [<Eoutmin>,<Eoutmax>],
   integrated over all angles.
     (options for muOut and phiOut ranges are not used).


<xml> : xml file name
<product> : exit product, such as  gamma or n or He4
<Ein>:  energy in MeV (lab), such as  10.0. 
<Eout*> similarly 
<reaction>: one of the reactions from reaction_list, such as '2n + Fe54'  (has blanks, so must use quotes)


CM and LAB values:
The results from ‘gnds_reaction_crosssection.py’ are for angle-integrated reactions, and are the same in CM and LAB.
 
The results from ‘gnds_reaction_angulardistribution.py’ are for 2-body reactions, and are CM angles and CM cross-sections.
 
The ‘production’ results from ‘gnds_inclusive_production.py’ and ‘gnds_reaction_production.py’ are in LAB angles and LAB cross-sections.


Examples:
using xml files at  C = /usr/workspace/ndg/translations/current   on LC.

gnds_reaction_list.py $C/endl2009.3/yi01/za008016.xml  

gnds_reaction_crosssection.py  $C/ENDF-VII.1/neutrons/n-026_Fe_056.endf.gnds.xml '2n + Fe55 + photon' 12.5

gnds_reaction_angulardistribution.py $C/ENDF-VII.1/neutrons/n-026_Fe_056.endf.gnds.xml '2n + Fe55 + photon' 9.5 n (FAILS)

gnds_reaction_production.py $C/ENDF-VIII/neutrons/n-008_O_016.endf.gnds.xml 'O17 + photon' gamma 10  (FAILS)

gnds_inclusive_production.py $C/ENDF-VIII/neutrons/n-008_O_016.endf.gnds.xml gamma 18 1.0 7.0 (GIVES 0.0)                                                                                                                                                                                                                                                   

Outputs:
gnds_reaction_list.py $C/endl2009.3/yi01/za008016.xml 
n + O16
n + (O16_e1 -> O16)
n + (O16_e2 -> O16)
n + (O16_e3 -> O16)
n + (O16_e4 -> O16)
n + (O16_e5 -> O16)
n + (O16_e6 -> O16)
n + (O16_e7 -> O16)
n + O16 + photon [continuum]
2n + O15
n + H1 + N15 + photon
n + He4 + C12 + photon
n + 4He4 + photon
2He4 + Be9 + photon
H1 + N16
H1 + (N16_e1 -> N16)
H1 + (N16_e2 -> N16)
H1 + (N16_e3 -> N16)
H2 + N15
H2 + (N15_e1 -> N15)
H2 + (N15_e2 -> N15)
H2 + (N15_e3 -> N15)
H2 + (N15_e4 -> N15)
H2 + (N15_e5 -> N15)
H2 + (N15_e6 -> N15)
H2 + (N15_e7 -> N15)
H2 + (N15_e8 -> N15)
H3 + N14
H3 + (N14_e1 -> N14)
H3 + (N14_e2 -> N14)
He4 + C13
He4 + (C13_e1 -> C13)
He4 + (C13_e2 -> C13)
He4 + (C13_e3 -> C13)
O17 + photon


gnds_reaction_list.py $C/ENDF-VII.1/neutrons/n-026_Fe_056.endf.gnds.xml
n + Fe56
n + (Fe56_e1 -> Fe56 + photon)
n + (Fe56_e2 -> Fe56 + photon)
n + (Fe56_e3 -> Fe56 + photon)
n + (Fe56_e4 -> Fe56 + photon)
n + (Fe56_e5 -> Fe56 + photon)
n + (Fe56_e6 -> Fe56 + photon)
n + (Fe56_e7 -> Fe56 + photon)
n + (Fe56_e8 -> Fe56 + photon)
n + (Fe56_e9 -> Fe56 + photon)
n + (Fe56_e10 -> Fe56 + photon)
n + (Fe56_e11 -> Fe56 + photon)
n + (Fe56_e12 -> Fe56 + photon)
n + (Fe56_e13 -> Fe56 + photon)
n + (Fe56_e14 -> Fe56 + photon)
n + (Fe56_e15 -> Fe56 + photon)
n + (Fe56_e16 -> Fe56 + photon)
n + (Fe56_e17 -> Fe56 + photon)
n + (Fe56_e18 -> Fe56 + photon)
n + (Fe56_e19 -> Fe56 + photon)
n + (Fe56_e20 -> Fe56 + photon)
n + (Fe56_e21 -> Fe56 + photon)
n + (Fe56_e22 -> Fe56 + photon)
n + (Fe56_e23 -> Fe56 + photon)
n + (Fe56_e24 -> Fe56 + photon)
n + (Fe56_e25 -> Fe56 + photon)
n + Fe56 + photon [continuum]
2n + Fe55 + photon
n + H1 + Mn55 + photon
Fe57 + photon
n + He4 + Cr52 + photon
H1 + Mn56 + photon
H2 + Mn55
H3 + Mn54
He3 + Cr54
He4 + Cr53 + photon
sumOfRemainingOutputChannels

gnds_reaction_list.py  $C/ENDF-VIII/neutrons/n-003_Li_006.endf.gnds.xml 
n + Li6
n + (Li6_e1 -> H2 + He4)
n + (Li6_e2 -> H2 + He4)
n + (Li6_e3 -> H2 + He4)
n + (Li6_e4 -> H2 + He4)
n + (Li6_e5 -> H2 + He4)
n + (Li6_e6 -> H2 + He4)
n + (Li6_e7 -> Li6 + photon)
n + (Li6_e8 -> H2 + He4)
n + (Li6_e9 -> H2 + He4)
n + (Li6_e10 -> H2 + He4)
n + (Li6_e11 -> H2 + He4)
n + (Li6_e12 -> H2 + He4)
n + (Li6_e13 -> H2 + He4)
n + (Li6_e14 -> H2 + He4)
n + (Li6_e15 -> H2 + He4)
n + (Li6_e16 -> H2 + He4)
n + (Li6_e17 -> H2 + He4)
n + (Li6_e18 -> H2 + He4)
n + (Li6_e19 -> H2 + He4)
n + (Li6_e20 -> H2 + He4)
n + (Li6_e21 -> H2 + He4)
n + (Li6_e22 -> H2 + He4)
n + (Li6_e23 -> H2 + He4)
n + (Li6_e24 -> H2 + He4)
n + (Li6_e25 -> H2 + He4)
n + (Li6_e26 -> H2 + He4)
n + (Li6_e27 -> H2 + He4)
n + (Li6_e28 -> H2 + He4)
n + (Li6_e29 -> H2 + He4)
n + (Li6_e30 -> H2 + He4)
n + (Li6_e31 -> H2 + He4)
Li7 + photon
2n + He4 + H1
H1 + He6
H3 + He4



gnds_reaction_crosssection.py  $C/ENDF-VII.1/neutrons/n-026_Fe_056.endf.gnds.xml '2n + Fe55 + photon' 12.5
partialProductionIntegral from evaluation file /usr/workspace/ndg/translations/current/ENDF-VII.1/neutrons/n-026_Fe_056.endf.gnds.xml with reaction '2n + Fe55 + photon' at energy 12.5
Incident E  12.5 has cross-section   0.12907500000000016


gnds_reaction_angulardistribution.py $C/ENDF-VII.1/neutrons/n-008_O_016.endf.gnds.xml 'He4 + C13' 6 He4 
partialProductionIntegral from evaluation file /usr/workspace/ndg/translations/current/ENDF-VII.1/neutrons/n-008_O_016.endf.gnds.xml with reaction 'He4 + C13' at energy 6.0
Incident E  6.0  to angle = 10.0  deg :: diff-xs  0.002437942691214376
Incident E  6.0  to angle = 20.0  deg :: diff-xs  0.001818711815961977
Incident E  6.0  to angle = 30.0  deg :: diff-xs  0.0013312144796630127
Incident E  6.0  to angle = 40.0  deg :: diff-xs  0.0011875929535997804
Incident E  6.0  to angle = 50.0  deg :: diff-xs  0.0012058128842558984
Incident E  6.0  to angle = 60.0  deg :: diff-xs  0.001113449480352453
Incident E  6.0  to angle = 70.0  deg :: diff-xs  0.0008776102787217605
Incident E  6.0  to angle = 80.0  deg :: diff-xs  0.0006687776812216692
Incident E  6.0  to angle = 90.0  deg :: diff-xs  0.0005783206895363768
Incident E  6.0  to angle = 100.0  deg :: diff-xs  0.0005136103477989951
Incident E  6.0  to angle = 110.0  deg :: diff-xs  0.0004069410596484066
Incident E  6.0  to angle = 120.0  deg :: diff-xs  0.0004112007149563378
Incident E  6.0  to angle = 130.0  deg :: diff-xs  0.0007587547183009655
Incident E  6.0  to angle = 140.0  deg :: diff-xs  0.0014318574409447288
Incident E  6.0  to angle = 150.0  deg :: diff-xs  0.0020898829285631824
Incident E  6.0  to angle = 160.0  deg :: diff-xs  0.0023966863672151783
Incident E  6.0  to angle = 170.0  deg :: diff-xs  0.0023669757970352927
Incident E  6.0  to angle = 180.0  deg :: diff-xs  0.002298197356822963

gnds_reaction_angulardistribution.py $C/ENDF-VIII/neutrons/n-003_Li_006.endf.gnds.xml 'H3 + He4' 1 H3
partialProductionIntegral from evaluation file /usr/workspace/ndg/translations/current/ENDF-VIII/neutrons/n-003_Li_006.endf.gnds.xml with reaction 'H3 + He4' at energy 1.0
Incident E  1.0  to angle = 10.0  deg :: diff-xs  0.032230059499488456
Incident E  1.0  to angle = 20.0  deg :: diff-xs  0.030527838162466403
Incident E  1.0  to angle = 30.0  deg :: diff-xs  0.028061050475821202
Incident E  1.0  to angle = 40.0  deg :: diff-xs  0.025258666802967993
Incident E  1.0  to angle = 50.0  deg :: diff-xs  0.022527471380679003
Incident E  1.0  to angle = 60.0  deg :: diff-xs  0.020163813464945366
Incident E  1.0  to angle = 70.0  deg :: diff-xs  0.01832648005513571
Incident E  1.0  to angle = 80.0  deg :: diff-xs  0.01706317108678031
Incident E  1.0  to angle = 90.0  deg :: diff-xs  0.016357467094670022
Incident E  1.0  to angle = 100.0  deg :: diff-xs  0.016161455614819884
Incident E  1.0  to angle = 110.0  deg :: diff-xs  0.016400080470037097
Incident E  1.0  to angle = 120.0  deg :: diff-xs  0.016960290006588698
Incident E  1.0  to angle = 130.0  deg :: diff-xs  0.017690948786623667
Incident E  1.0  to angle = 140.0  deg :: diff-xs  0.018429104597993646
Incident E  1.0  to angle = 150.0  deg :: diff-xs  0.019043591906457893
Incident E  1.0  to angle = 160.0  deg :: diff-xs  0.019468529002271764
Incident E  1.0  to angle = 170.0  deg :: diff-xs  0.019702817702139604
Incident E  1.0  to angle = 180.0  deg :: diff-xs  0.019775543690602054


gnds_reaction_production.py $C/ENDF-VIII/neutrons/n-008_O_016.endf.gnds.xml 'O17 + photon' photon 10
partialProductionIntegral from evaluation file /usr/workspace/ndg/translations/current/ENDF-VIII/neutrons/n-008_O_016.endf.gnds.xml with reaction 'O17 + photon' giving photon
Incident Ein= 1.0  to Eout= [0.0, 1.0]  MeV has partialProductionIntegral= 8.19e-05
Incident Ein= 2.0  to Eout= [0.0, 2.0]  MeV has partialProductionIntegral= 5.39514e-05
Incident Ein= 3.0  to Eout= [0.0, 3.0]  MeV has partialProductionIntegral= 9.542591216e-05
Incident Ein= 4.0  to Eout= [0.0, 4.0]  MeV has partialProductionIntegral= 0.00013178556196
Incident Ein= 5.0  to Eout= [0.0, 5.0]  MeV has partialProductionIntegral= 0.0001621989561
Incident Ein= 6.0  to Eout= [0.0, 6.0]  MeV has partialProductionIntegral= 0.00017047669770000001
Incident Ein= 7.0  to Eout= [0.0, 7.0]  MeV has partialProductionIntegral= 0.00014224384518
Incident Ein= 8.0  to Eout= [0.0, 8.0]  MeV has partialProductionIntegral= 0.0001221921549
Incident Ein= 9.0  to Eout= [0.0, 9.0]  MeV has partialProductionIntegral= 0.00010954268464
Incident Ein= 10.0  to Eout= [0.0, 10.0]  MeV has partialProductionIntegral= 0.00010015131223999999


gnds_reaction_production.py $C/ENDF-VII.1/neutrons/n-026_Fe_056.endf.gnds.xml '2n + Fe55 + photon' n 18
partialProductionIntegral from evaluation file /usr/workspace/ndg/translations/current/ENDF-VII.1/neutrons/n-026_Fe_056.endf.gnds.xml with reaction '2n + Fe55 + photon' giving n
Incident Ein= 1.8  to Eout= [0.0, 1.8]  MeV has n partialProductionIntegral= 0.0
Incident Ein= 3.6  to Eout= [0.0, 3.6]  MeV has n partialProductionIntegral= 0.0
Incident Ein= 5.4  to Eout= [0.0, 5.4]  MeV has n partialProductionIntegral= 0.0
Incident Ein= 7.2  to Eout= [0.0, 7.2]  MeV has n partialProductionIntegral= 0.0
Incident Ein= 9.0  to Eout= [0.0, 9.0]  MeV has n partialProductionIntegral= 0.0
Incident Ein= 10.8  to Eout= [0.0, 10.8]  MeV has n partialProductionIntegral= 0.0
Incident Ein= 12.6  to Eout= [0.0, 12.6]  MeV has n partialProductionIntegral= 0.266245849728
Traceback (most recent call last):



gnds_inclusive_production.py $C/ENDF-VIII/neutrons/n-008_O_016.endf.gnds.xml photon 18 6.0 7.0   
partialProductionIntegral from evaluation file /usr/workspace/ndg/translations/current/ENDF-VIII/neutrons/n-008_O_016.endf.gnds.xml giving photon
Args: photon 18.0 [6.0, 7.0]
Incident Ein= 18.0  to Eout= [6.0, 7.0]  MeV has partialProductionIntegral= 0.172639207509050

