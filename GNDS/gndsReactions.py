#! /usr/bin/env python3
from __future__ import print_function

# <<BEGIN-copyright>>
# <<END-copyright>>

# This script takes a gnd/XML file, reads it in and rewrites it to a file with the same name 
# (only in the currently directory) with the extension '.g2g' added. The intent of this script 
# is for testing the reading/writing of gnd/XML file. More than one file can be inputted, and each
# is read/written independently.

import sys
from fudge.gnds import reactionSuite

for gndFile in sys.argv[1:] :
    gnd = reactionSuite.readXML( gndFile )
    gnd.convertUnits( {'eV':'MeV'} )
    p = gnd.projectile; t = gnd.target; incoming = p + ' + ' +t

    for reac in gnd.reactions:
        Qi = None; Qf = None; emin = None
        Qi = reac.getQ('MeV', final=False)
        Qf = reac.getQ('MeV', final=True)
        try:
            emin = reac.getQ.domainMin('MeV')
            print(gndFile,'::   ',reac.label,',   Q(i,f) =',Qi,Qf,' domainMin=',emin,' (MT =',reac.ENDF_MT,')')
        except:
            print(gndFile,'::   ',reac.label,',   Q(i,f) =',Qi,Qf,' (MT =',reac.ENDF_MT,')')
        if True: # and reac.label != incoming and reac.label!='sumOfRemainingOutputChannels':
            if reac.ENDF_MT==2: continue
            xsc = reac.crossSection[ 'eval' ]
            xFile = gndFile + '-' + reac.label.replace(' ','')
            with open(xFile,'w') as fout:
                try: xsc.clip( rangeMin=1e-20 )
                except: pass
                fout.write( xsc.toPointwise_withLinearXYs(lowerEps=1e-8).toString() )
