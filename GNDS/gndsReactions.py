#! /usr/bin/env python3
from __future__ import print_function

# <<BEGIN-copyright>>
# <<END-copyright>>


import sys
from fudge.gnds import reactionSuite

for gndFile in sys.argv[1:] :
    gnd = reactionSuite.readXML( gndFile )
    gnd.convertUnits( {'eV':'MeV'} )
    p = gnd.projectile; t = gnd.target; incoming = p + ' + ' +t

    for reac in gnd.reactions:
        try:
            Qi = reac.getQ('MeV', final=False)
            Qf = reac.getQ('MeV', final=True)
            print(gndFile,'::   ',reac.label,',   Q(i,f) =',Qi,Qf,' (MT =',reac.ENDF_MT,')')
        except:
            print(gndFile,'::   ',reac.label,' (MT =',reac.ENDF_MT,')')
        xsc = reac.crossSection[ 'eval' ]
        xFile = gndFile + '-' + reac.label.replace(' ','')
        with open(xFile,'w') as fout:
            try: xsc.clip( rangeMin=1e-20 )
            except: pass
            fout.write( xsc.toPointwise_withLinearXYs(lowerEps=1e-8).toString() )
