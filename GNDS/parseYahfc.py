#!/usr/bin/env python3

import argparse,os,re,datetime,sys,math

from PoPs import alias as PoPsAliasModule
from PoPs.groups import misc as popsNamesModule
from PoPs import database as databasePoPsModule
from PoPs import IDs as idsPoPsModule
from PoPs.groups import isotope as isotopeModule
from PoPs.groups import chemicalElement as chemicalElementModule
from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import lepton as leptonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nuclide as nuclideModule
from PoPs.quantities import charge as chargeModule    
from PoPs.quantities import mass as massModule    
from PoPs.quantities import spin as spinModule    
from PoPs.quantities import parity as parityModule    
from PoPs.quantities import halflife as halflifeModule    
from PoPs.quantities import nuclearEnergyLevel as nuclearEnergyLevelModule
from PoPs.decays import decayData as decayDataModule
from PoPs.decays import probability as probabilityModule
from PoPs.decays import product as popsProductModule

from pqu import PQU as PQUModule

import xData.standards as standardsModule
import xData.XYs as XYsModule
import xData.multiD_XYs as multiD_XYsModule
from xData.Documentation import documentation as documentationModule
from xData.Documentation import computerCode  as computerCodeModule
from xData import date

from fudge.gnds import physicalQuantity as physicalQuantityModule
from fudge.gnds import reactionSuite as reactionSuiteModule
from fudge.gnds import styles as stylesModule
from fudge.gnds import channels as channelsModule
from fudge.gnds import product as productModule
from fudge.gnds.reactions import reaction as reactionModule
from fudge.gnds.reactionData import crossSection as crossSectionModule
from fudge.gnds.channelData import Q as QModule
from fudge.gnds.productData import multiplicity as multiplicityModule
from fudge.gnds.productData.distributions import angular as angularModule
from fudge.gnds.productData.distributions import LLNL_angularEnergy as angularEnergyModule
from fudge.gnds.productData.distributions import reference as referenceModule
from fudge.gnds.productData.distributions import unspecified as unspecifiedModule
from fudge.gnds.productData.distributions import branching3d as branching3dModule
from fudge.gnds.productData.distributions import energyAngular as energyAngularModule

from fudge.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc as endfMiscModule 
from fudge.lib import GNDSType as GNDSTypeModule

FUDGE_EPS           = 1e-8
energyUnit          = "MeV"
styleName           = 'eval'
evaluationLibrary   = 'LLNL2020'
evaluationVersion   = '1'
evaluationEmax      = '20.0'

defaultProjectile   = 'n'
defaultVerbose      = 0
defaultFalse        = 'False'
defaultNone         = 'False'
defaultEndl         = ''
defaultInputDeck = 'YAHFC-commands.txt'

# Physics for calculating Rutherford elastic cross-sections:
amu = 931.494013
ifsc = 137.03599976
hbc = 197.32705
pi = 3.1415826536

fmscal = 2.0  * amu / hbc**2
coulcn = hbc/ifsc
etacns = coulcn * fmscal**0.5 * 0.5

maxMuCutoff         = 0.96
thin_eps            = 1e-4


### record of script calls to convert old endl files to GNDS for splicing 
###     python /path/to/your/fudge/brownies/LLNL/ENDL2gnds.py 95241 -l endl2009.3 -p n -o endl9p3_yi01za095241
### and the call to convert the new yahfc files into GNDS
###     python parseYahfc.py 241Am --endl endl9p3_yi01za095241.xml 

parser = argparse.ArgumentParser("Translate YAHFC to GNDS")
parser.add_argument('target', type=str, help='target nuclide' )
parser.add_argument('-p', '--projectile', type=str, default = defaultProjectile, help='projectile Default = %s'%defaultProjectile )
parser.add_argument('-E', '--Emax', type=float, default = evaluationEmax, help='Max projectile energy. Default=%s' % evaluationEmax)
parser.add_argument( '--endl', type=str, default = None, help='existing endl library to splice into. Default = %s'%defaultEndl )
parser.add_argument('-v', '--verbose', action='count', default = defaultVerbose, help='  Default = %d'%defaultVerbose )
parser.add_argument('--saved',  action='store_true', default = None, help=' Use a pickle file to save the data in an intermediate step, once read.' )
parser.add_argument('--resave',  action='store_true', default = None, help=' Use a pickle file to save the data in an intermediate step, once read.' )
parser.add_argument( '--ENDL_I_1_3',  action='store_true', default = None, help = '''Converts P(E',mu|E) distributions to ENDL I = 1 and 3 data for outputting.''' )
parser.add_argument( '--inputDeck', type=str, default = defaultInputDeck, help='Yahfc input parameter deck. Default = %s'%defaultInputDeck )

args = parser.parse_args()

print('\nRun:',' '.join([sys.argv[0].split('/')[-1]]+sys.argv[1:]),'\n')
aliasNameDict = {}

CWD = os.getcwd()
YAHFC_DATA = os.getenv('YAHFC_DATA')
if args.target[-1] == '/': args.target = args.target[:-1]

def now():
    return date.Date( resolution = date.Resolution.time )

def collateByFirstElement( data ) :

        #   handles  0.0              > 0.0            < 0.0 values.
    prior = min( data[0][0] - 1.0, 0.9 * data[0][0], 1.1 * data[0][0] )
    collatedList = []
    subList = []

    for datum in data :
        domain = datum[0]
        theRest = datum[1:]
        if( prior != domain ) :
            collatedList.append( [ prior, subList ] )
            subList = []
            prior = domain
        subList.append( theRest )
    collatedList.append( [ prior, subList ] )
    collatedList.pop( 0 )

    return( collatedList )

def main():

    ### get data from Erich's files
    channels, decayDataStructure = getEvaluationData() 
    
    PoPs = fillPops(decayDataStructure)
    PoPs.addFile(os.path.join(YAHFC_DATA, 'lightNuclei_pops.xml'))
#    print('Constructed PoPs')
#    print('PopsData:',PoPs.keys())    
    
    RxnSte = sortChannels(channels,PoPs,decayDataStructure)
    
    
    parameterFile = open(os.path.join(CWD,args.inputDeck),'r')
    inputParameters = parameterFile.readlines()
    print('Read parameter file',os.path.join(CWD,args.inputDeck),'with',len(inputParameters),'lines')
    inputDataSpecs = computerCodeModule.InputDeck( "Yahfc input commands" , ('\n  %s\n' % now() )  + (''.join( inputParameters ))+'\n' ) 
    computerCodeModel = computerCodeModule.ComputerCode( label = 'Yahfc fit', name = 'YAHFC', version = '', date = now() )
    computerCodeModel.inputDecks.add( inputDataSpecs )  
    
    computerCodes = RxnSte.documentation.computerCodes
    computerCodes.add( computerCodeModel )
#     print(RxnSte.toXML())
    
    
    ### now load old endl and convert
    if args.endl != None:
        from brownies.legacy.toENDL import reactionSuite
        args.endl = os.path.realpath(args.endl)
        oldReactionSuite = GNDSTypeModule.read( args.endl )
        RxnSte = getMissingEvaluationParts(RxnSte,oldReactionSuite)
    
    targetS,targetZ,targetA,ZAtarget = getSymAZfromIsotope(args.target)
    yi = 'npdtha'.find(args.projectile)+1
    protare = 'yi%02dza%03d%03d.yahfc.xml'%(yi,targetZ,targetA)
    print('\nWriting protare  ',protare)
    
    RxnSte.saveToFile(os.path.join(CWD,protare) , xs_pdf_cdf1d_singleLine = True )
    
    #RxnSte.check()
    
    if args.endl != None:
        RxnSte.toENDL( os.path.join(CWD,'yi%02dza%03d%03d.yahfc-ENDL/asciii'%(yi,targetZ,targetA)) )
    
def sortChannels(channels,PopsData,decayData):
#     print("sortChannels. PopsData for",list(PopsData.keys()))
    crossSectionAxes = crossSectionModule.defaultAxes( energyUnit )
    multiplicityAxes = multiplicityModule.defaultAxes( energyUnit )
    QAxes = QModule.defaultAxes( energyUnit )

    ### get projectile and target
    targetS,targetZ,targetA,ZAtarget = getSymAZfromIsotope(args.target)
    args.target = popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( targetS, targetA )
    target = PopsData[args.target]
    projectile = PopsData[getPopsPartId(args.projectile)]
    proj = projectile
    if hasattr(proj, 'nucleus'): proj = proj.nucleus
    ZpZt = targetZ * proj.charge[0].value
#    print('ZpZt:', ZpZt, targetZ, projectile.charge[0].value)

    thisTemp = physicalQuantityModule.temperature( 0 , 'MeV/k' )
    E_ratio = target.getMass('amu') /(projectile.getMass('amu') + target.getMass('amu'))
    reduced_mass = projectile.getMass('amu')*E_ratio

    projectileDomain = stylesModule.projectileEnergyDomain( 1e-11, 20, 'MeV' )  # will be overwritten after parsing all reactions
    evaluation = "%s-%s" % ( evaluationLibrary, evaluationVersion )
    datenow = datetime.datetime.today()
    maxDate = '%04d-%02d-%02d'%(int(datenow.year),int(datenow.month),int(datenow.day))
    evaluatedStyle = stylesModule.evaluated( styleName, '', thisTemp, projectileDomain, evaluationLibrary, evaluationVersion, date = maxDate )
    reactionSuite = reactionSuiteModule.reactionSuite( projectile.id, target.id, evaluation, style = evaluatedStyle, PoPs = PopsData, interaction = 'nuclear', formatVersion = '2.0.LLNL_3' )
    
    NTdict = defineInitialNTs(args.projectile)
    
    reactionStoreDict = {}
    reactionStoreMultiples = {}
    for k in channels:
        kbase = k.replace('_g','').split('_')    
        particles = kbase[0]
        level = kbase[1][1:] if len(kbase)>1 else '0'
        digits = re.search('(\d+)',particles)                 # numbers
        outgoingParts = re.search('(\D+)',particles).group()  # names
        particleMultiplicity = {}
        if digits is None:  # no digits: sequential products
            NTbase = 0
            for p in ['g','a','h','t','d','p','n']:   # make NT depending  on numbers of particles (but not on gammas)
                multiplicity = particles.count(p)
                NTbase = NTbase*10 + multiplicity
                if multiplicity > 0: particleMultiplicity[p] = multiplicity
            NT = NTbase*1000 + int(level)
            
        else:     # digits, so products already grouped eg "5nd2p"
            particle = ''
            multiplicity = 0
            for d in reversed(list(particles)):
                if d.isalpha():
                    if multiplicity > 0: particleMultiplicity[particle] = multiplicity
                    particle = d
                    multiplicity = 1
                else:  # digit
                    multiplicity = int(d)
            particleMultiplicity[particle] = multiplicity
            NTbase = 0
            for p in ['g','a','h','t','d','p','n']:   # make NT depending  on numbers of particles (but not on gammas)
                NTbase = NTbase*10 + particleMultiplicity.get(p,0)
            NT = NTbase*1000 + int(level)                    
        
        if '_g' in k: NT = (NT // 1000) * 1000 + 900  # continuum products
        if 'f' in k: NT = 18 # fission
        print('\nChannel k = ',k,'gives level=',level,' NT=',NT,'with multiplicities',particleMultiplicity)
        
        chan = channels[k]
        ENDF_NT = NTdict.get(k, NT )   # use NT numbers from above if missing in initial dictionary
        if ENDF_NT==NT:  NTdict[k] = NT
        cp_elastic =   ZpZt != 0 and ENDF_NT == 2 
        residualPartner = None
        
        if args.verbose>1: print('\n',k,' ===> ', chan.get('residual','?'),'NT:',ENDF_NT)
        
        ### channel is 2-body, N-body, or fission
        outgoingParts = list(outgoingParts)
        channelName =   outgoingParts.copy()
        gammaMult = getMultiplicityFromData(chan,'g') 
        if gammaMult is not None and 'g' not in outgoingParts: outgoingParts.append('g')
#         print('channelName',channelName,'outgoingParts',outgoingParts)
        
        genre = channelsModule.Genre.NBody
        if( ( k[0] not in ['f','g'] ) and ( len(  outgoingParts ) == 1 ) ) : genre = channelsModule.Genre.twoBody
        reaction = reactionModule.reaction( genre, ENDF_NT )

        outputChannel = reaction.outputChannel
        if( k[0] == 'f' ) : outputChannel.fissionGenre = channelsModule.fissionGenre.total
        if k.startswith('%s_g'%args.projectile): outputChannel.process = channelsModule.processes.continuum 
 
        crossSection = getCrossectionFromData(chan,k)
        if crossSection == None: # and ENDF_NT != 2:
            print('   Skip channel giving',''.join(outgoingParts),' (NT',ENDF_NT,') with missing or zero crossSection.')
            continue
        else:
            print('   Include channel giving',''.join(outgoingParts),' (NT',ENDF_NT,')')

        XSmin,XSmax = crossSection.domainMin, crossSection.domainMax
        axesMult = multiplicityModule.defaultAxes( crossSection.domainUnit )
        if not cp_elastic: 
            reaction.crossSection.add( crossSection )
        
        Qval = getQFromData(chan) if ENDF_NT != 2 else 0.
        Qaxes = QModule.defaultAxes( crossSection.domainUnit )
        Qconst = QModule.constant1d( Qval, crossSection.domainMin, crossSection.domainMax, axes = Qaxes, label = styleName ) 
        outputChannel.Q.add( Qconst )  

        if k=='f':
            reaction.updateLabel()
            reactionSuite.reactions.add(reaction)
            reactionStoreDict[NTdict[k]] = reaction 
#             print('Added',reactionStoreDict[NTdict[k]].label,'for NT=',NTdict[k])
            continue

        productList = outgoingParts
#         print(k,'productList',productList,'particleMultiplicity=',particleMultiplicity)
        found = []
        for productName in productList:
            popsProdName = getPopsPartId(productName)
            product = productModule.product( popsProdName, label = popsProdName) #,  outputChannel = outputChannel )
            
            if ENDF_NT==2: 
                multiplicityVal = 1 if productName == args.projectile else 0
                print('   elastic product',productName,' * ',multiplicityVal)         
                
            elif particleMultiplicity.get(productName,None) is not None:
                multiplicityVal = particleMultiplicity[productName]   
            else: 
                multiplicityVal = getMultiplicityFromData( chan, productName )
            if multiplicityVal is None:
                print('ERROR: Missing multiplicity for',productName,'in NT=',ENDF_NT,'from',k)
                continue
#             if productNameprint('    -- product',productName,'=',popsProdName,' has m=',multiplicityVal,'in',ENDF_NT)
            
            if( isinstance( multiplicityVal, ( int, float ) ) ) :
                if args.verbose>2: print('multiplicity constant',multiplicityVal)
                multiplicity = multiplicityModule.constant1d( multiplicityVal, crossSection.domainMin, crossSection.domainMax, axes = axesMult, label = styleName )
            else:
                multiplicity = multiplicityModule.XYs1d( data = multiplicityVal, dataForm = 'xsandys',  axes = axesMult, label = styleName ) 
                multiplicity = multiplicity.domainSlice(domainMin = XSmin, domainMax = XSmax)
                multiplicity.label = styleName
            product.multiplicity.add( multiplicity )
            
            PmuE = getAngGivenEin(channelName,chan,productName,crossSection)
            if PmuE is None: 
                found = []
                continue
            
            if cp_elastic:  # convert ratio-to-Rutherford to difference-to-Rutherford  for GNDS. Rewrite crossSection and PmuE.
                from fudge.gnds.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as  CoulombPlusNuclearElasticModule
                from fudge.gnds.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import nuclearPlusInterference as nuclearPlusInterferenceModule
                import fudge.gnds.productData.distributions as distributionsModule

                crossSectionAxes = crossSectionModule.defaultAxes( energyUnit )
                angularAxes = distributionsModule.angular.defaultAxes( 'MeV' )
                effDist = distributionsModule.angular.XYs2d( axes = angularAxes )
                
                effXsc  = []
                muCutoff = -1.
                for Elab in PmuE.domainGrid:
                    Ecm = Elab*E_ratio
                    SRratios = PmuE.evaluate(Elab)
                    wvk = (fmscal * reduced_mass * Ecm)**0.5
                    eta = ZpZt * etacns * (reduced_mass/Ecm)**0.5
                    sigdiffs = []
                    for mu in SRratios.domainGrid:
                        if mu > maxMuCutoff: continue
                        SRratio = SRratios.evaluate(mu)
                        muCutoff = max(mu,muCutoff)
                        
                        Rutherford = (eta/(wvk*(1.-mu)))**2 * 0.01 *2*pi  # convert to barns
                        sigdiffs.append([mu,Rutherford * (SRratio-1.)])
                    
                    angdist = distributionsModule.angular.XYs1d( data = sigdiffs, dataForm = 'xys',  outerDomainValue = Elab, axes = angularAxes ) 
                    norm = angdist.integrate()
                    if norm == 0.0: norm = 1.0  # zero angdist == pure Rutherford
                    if args.verbose>1: print('Elab=',Elab,' cm=',Ecm, ' sig-Rutherford norm :',norm)
                    angdist = angdist/norm  
                    effDist.append( angdist.thin( thin_eps ) )
                    effXsc.append( float(norm) )

                effXsc = crossSectionModule.XYs1d( axes = crossSectionAxes, data=(PmuE.domainGrid,effXsc), dataForm="XsAndYs" )     
                            
                NCPI = nuclearPlusInterferenceModule.nuclearPlusInterference( muCutoff=muCutoff,
                    crossSection=nuclearPlusInterferenceModule.crossSection( effXsc ),
                    distribution=nuclearPlusInterferenceModule.distribution( effDist ) )

                CoulombElastic = CoulombPlusNuclearElasticModule.form( projectile.id, styleName, nuclearPlusInterference = NCPI, identicalParticles=False )
                reaction.doubleDifferentialCrossSection.add( CoulombElastic )
                
                reaction.crossSection.add( crossSectionModule.CoulombPlusNuclearElastic( link = reaction.doubleDifferentialCrossSection[styleName],
                    label = styleName, relative = True ) )

            if k==productName and args.verbose>2:
                print( 'Outputs from', k, 'for', productName, ' PmuE on grid:', PmuE.domainGrid )

            PEpMuE, xys3d, frame = getEoutGivenAngEin( channelName, chan, productName, crossSection )
            if PEpMuE is None:
                if not cp_elastic: 
                    form = angularModule.twoBodyForm(styleName, standardsModule.frames.centerOfMassToken, PmuE )  ### all twobody reactions use the form : angularModule.twoBodyForm()
                else:
                    form =  referenceModule.CoulombPlusNuclearElastic( link = reaction.doubleDifferentialCrossSection[styleName], label = styleName, relative = True )
                residualPartner = form
                product.distribution.add( form )
            else:
                if( args.ENDL_I_1_3 ) :
                    angularSubform = angularEnergyModule.LLNLAngularOfAngularEnergySubform( PmuE )
                    angularEnergySubform = angularEnergyModule.LLNLAngularEnergyOfAngularEnergySubform( PEpMuE )
                    form = angularEnergyModule.LLNLAngularEnergyForm( styleName,  standardsModule.frames.labToken, angularSubform, angularEnergySubform )
                else :
                    form = energyAngularModule.form( styleName, frame, xys3d )
                product.distribution.add( form )
                
            outputChannel.products.add(outputChannel.products.uniqueLabel(product))
            found.append(productName)
        
        if len(found)==0: 
            if args.verbose>0: print('None of products found:',productList)
            continue
        if args.verbose>2: print('Products found:',found,'of',productList)
        ### add the residual channels[k]['residual'] , distribution is unknown or recoil for twobody
#         print('channels:',channels[k])
        try:
            residual =  PopsData[channels[k]['residual']]
        except: 
            residualPartner = None
        product = productModule.product( getAliasOf(residual.id), label = getAliasOf(residual.id) )
        multiplicity = multiplicityModule.constant1d( 1, crossSection.domainMin, crossSection.domainMax, axes = axesMult, label = styleName )
        product.multiplicity.add( multiplicity )
        if residualPartner is not None and len(outgoingParts) == 1 : 
            form = angularModule.twoBodyForm(styleName, standardsModule.frames.centerOfMassToken, angularModule.recoil(residualPartner) ) 
        else:
            if len(outgoingParts) == 1: print('   Missing recoil residualPartner')
            form = unspecifiedModule.form(styleName)
        product.distribution.add( form )
        
        levelEnergy = residual.energy[0].value
        if levelEnergy > 0.0 and not isIsomer(residual,decayData) :
                
            decayChannel = channelsModule.outputChannel( channelsModule.Genre.NBody )
            Qres = QModule.constant1d( levelEnergy, crossSection.domainMin, crossSection.domainMax, axes = Qaxes, label = styleName ) 
            decayChannel.Q.add( Qres ) 
            
            decayResidualName = reactionSuite.PoPs[residual.id].isotope.symbol
            decayResidual = productModule.product( id = decayResidualName, label = decayResidualName )
            decayResidual.multiplicity.add( multiplicityModule.constant1d( 1, crossSection.domainMin, crossSection.domainMax, axes = axesMult, label = styleName ) )
            decayResidual.distribution.add( unspecifiedModule.form(styleName) )
            decayChannel.products.add( decayResidual )

            decayPhoton = productModule.product( id = idsPoPsModule.photon, label = idsPoPsModule.photon )
            decayPhoton.multiplicity.add( multiplicityModule.branching1d( styleName ) )
            decayPhoton.distribution.add( branching3dModule.form( styleName, standardsModule.frames.labToken ) )
            decayChannel.products.add( decayPhoton )

            product.addOutputChannel( decayChannel )

        outputChannel.products.add(outputChannel.products.uniqueLabel(product))

        reaction.updateLabel()
        if k not in NTdict.keys():
            print('Key',k,'not in',NTdict.keys())
        if NTdict[k] not in list(reactionStoreDict.keys()):
            reactionStoreMultiples[NTdict[k]] = 1
            reactionStoreDict[NTdict[k]] = reaction
        else:
            print('Merge with previous NT=',NTdict[k],'reaction')
            reactionStoreMultiples[NTdict[k]] += 1
            reactionStoreDict[NTdict[k]] = addToExistingReaction( reactionStoreDict[NTdict[k]], reaction, productList )   
    
    NTs = niceSortOfNTs( list( reactionStoreDict.keys( ) ) )
    print( 'Sorted keys:', NTs,'\n')

    print( '\nChecking normalization' ) 
    for NT in NTs :
        if args.verbose>-1: print( 'Store NT=', NT, 'for', reactionStoreDict[NT].label )
        if NT not in [1,18]:
            reactionStoreDict[NT] = checkNormalization( reactionStoreDict[NT], reactionStoreMultiples[NT] )
            reactionSuite.reactions.add( reactionStoreDict[NT] )

    if ZpZt == 0:
        allChannels = [reaction.label for reaction in reactionSuite.reactions ]
        if len(allChannels)==0:
            print("There are no (non-zero) reactions channels!!!")
            sys.exit()
        NT = 1
        label = 'total'
        totalreac = reactionSuite.buildCrossSectionSum( label, NT, 'eval', allChannels, 0.0)
        reactionSuite.sums.crossSections.add  ( totalreac )
        if args.verbose>-1: print( 'Store sum NT=', NT, 'for', label )
             
    if( 2 not in NTs ) : print( '\nNo elastic NT=2 found!' )   
#     reactionSuite.check()
    return reactionSuite

def getAliasOf(productStr):
    if productStr in list(aliasNameDict.keys()) :
        return aliasNameDict[productStr]
    else:
        return productStr  

def isIsomer(residual,decayData) :
    if '_e' not in residual.id:
        return False
    levelIndex = int(residual.id.split('_e')[1])
    isoName = residual.id.split('_e')[0]
    S,Z,A,ZA = getSymAZfromIsotope(isoName)
    isoName = '%d%s'%(A,S)
    if decayData[isoName][levelIndex]['isomer']=='True':
        return True
    else:
        return False

def addToExistingReaction( oldRxn, newRxn, productList ) :
        
    ### are Q vals the same?
    ### if the Qs are not the same then we shouldn't be merging these reactions 
    if oldRxn.getQ('MeV') != newRxn.getQ('MeV'):
        print( 'Attempting to merge %s and %s ' % ( oldRxn.label, newRxn.label ) )
        print('Q values : ', oldRxn.getQ('MeV'), newRxn.getQ('MeV'))
        raise Exception('these two channels have different Qs and should not be merged. Please assign different NT numbers ')
    
    ### add crossection
    oldXS = oldRxn.crossSection[styleName]
    newXS = newRxn.crossSection[styleName]
    summedXS = oldXS + newXS
    summedXS.label = oldRxn.crossSection[styleName].label
    oldRxn.crossSection.replace(summedXS)

    ### add distributions for each product, except residual
    for productName in productList:
        popsProdName = getPopsPartId(productName)
        oldProd = oldRxn.outputChannel.getProductWithName(popsProdName)
        newProd = newRxn.outputChannel.getProductWithName(popsProdName)
        oldI1 = oldProd.distribution[0].subforms[0].data
        newI1 = newProd.distribution[0].subforms[0].data
        
        for func in oldI1.functionals:
            Epnt = func.outerDomainValue
            oldXSVal = oldXS.evaluate(Epnt)
            newXSVal = newXS.evaluate(Epnt)
            if oldXSVal is not None and newXSVal is not None:
                func = func + newI1.evaluate(Epnt)*newXS.evaluate(Epnt)*newProd.multiplicity.evaluate(Epnt)
        
        oldI3 = oldProd.distribution[0].subforms[1].data
        newI3 = newProd.distribution[0].subforms[1].data
        for func2 in oldI3.functionals:
            for func in func2:
                muPnt = func.outerDomainValue
                Epnt = func2.outerDomainValue
                newXSVal = newXS.evaluate(Epnt)
                if newXSVal is not None:
                    newfunc = newI3.evaluate(Epnt).evaluate(muPnt)
                    func, newfunc = func.mutualify(FUDGE_EPS,FUDGE_EPS,True,newfunc,FUDGE_EPS,FUDGE_EPS,True)
                    func = func + newfunc*newXSVal*newProd.multiplicity.evaluate(Epnt)
                
    return oldRxn
    
def checkNormalization(rxn,mult):

    def removeValues(subform,energiesToRemove,msg):
        if len(energiesToRemove) > 0:
            for EinToRemove,AngToRemove in energiesToRemove:
                for i2,f2 in enumerate(subform) :
                    if abs(f2.outerDomainValue-EinToRemove) < FUDGE_EPS :
                        if AngToRemove=='all':  ## kill the whole Ein
                            old = subform.pop(i2)
                            break
                        else:
                            if subform.dimension == 2 :
                                for i1,f1 in enumerate(f2):
                                    if f1[0] == AngToRemove:
                                        old = f2.pop(i1)
                                        break
                            elif subform.dimension == 3 :
                                for i1,f1 in enumerate(f2):
                                    if f1.outerDomainValue == AngToRemove:
                                        old = f2.pop(i1)
                                        break
    
#     if mult>1:  ##TEMP
#         rxn.crossSection[styleName] = rxn.crossSection[styleName]/mult 

    for product in rxn.outputChannel.products:
        productStr = product.label
        
        energiesToRemove = []                
        for form in product.distribution[0].subforms:
            if( isinstance( form, ( angularModule.recoil, energyAngularModule.XYs3d ) ) ) : continue
            if hasattr(form,'dimension'): 
                subform = form
            else:
                subform = form.data  
            
            if subform.dimension == 2 :
                for f1 in subform:
                    currNorm = f1.integrate().getValue() 
                    if  currNorm > 0. :  
                        f1 = f1.normalize(dimension=1)
                    else:
                        energiesToRemove.append((f1.outerDomainValue,'all'))
                        
            elif subform.dimension == 3 :
                for i2,f2 in enumerate(subform) :
                    thisNorm = 0.
                    for i1,f1 in enumerate(f2):
                        currNorm = f1.integrate().getValue() 
                        thisNorm += currNorm
                        if currNorm > 0. :  
                            subform[i2][i1] = f1.normalize(dimension=1) 
                        else:
                            ### remove all zero P(Ep|E,mu) distriobutions for any combination of Ein, mu
                            ### remove entire Ein for any without endpoints of -1 and 1
                            energiesToRemove.append((f2.outerDomainValue,f1.outerDomainValue))  ### Ein,angle to remove
                            if abs(f1.outerDomainValue) == 1. : energiesToRemove.append((f2.outerDomainValue,'all'))

        energiesToRemove = list(set(energiesToRemove))
        if len(energiesToRemove) > 0 :
            for form in product.distribution[0].subforms:
                if isinstance(form,angularModule.recoil) : continue
                if hasattr(form,'dimension'): 
                    subform = form
                else:
                    subform = form.data  
                if subform.dimension == 2 :
                    removeValues(subform,energiesToRemove,productStr+' XYs2D Ein')
                elif subform.dimension == 3 :
                    removeValues(subform,energiesToRemove,productStr+' XYs3D Ein')

    return rxn
 
def getMissingEvaluationParts(RxnSte,oldReactionSuite):
    
    ### get low energy xs for several channels (capture, elastic, fission)
    ### will need to mutualify and blur the edges before joining
    for channelName in ['elastic','capture','fission','Am242_e2 + photon']:
        aliasChannelName = channelName
        for k in list(aliasNameDict.keys()):
            aliasChannelName = aliasChannelName.replace(k,aliasNameDict[k])
        newRxn = RxnSte.getReaction( aliasChannelName )
        oldRxn = oldReactionSuite.getReaction( channelName )
        if newRxn is None :
            print('cant find new version of this reaction : ',aliasChannelName)
            continue
        if oldRxn is None :
            print('cant find old version of this reaction : ',channelName)
            continue
        newXS = newRxn.crossSection[styleName].trim()
        oldXS = oldRxn.crossSection[styleName]
        newMin,newMax = newXS.domainMin,newXS.domainMax
        oldMin,oldMax = oldXS.domainMin,oldXS.domainMax
        oldXS = oldXS.domainSlice(oldMin,newMin)
        oldXS,newXS = oldXS.mutualify(FUDGE_EPS,FUDGE_EPS,True,newXS,FUDGE_EPS,FUDGE_EPS,True)        
        newXS = newXS+oldXS
        newXS.label = styleName
        newRxn.crossSection.replace(newXS)
    
    ### get fission parts : multiplicity, Q(E) (I=12), distributions 
    oldFission = oldReactionSuite.getReaction( 'fission' )  ### NT num for fission
    newFission = RxnSte.getReaction( 'fission' )            ### NT num for fission
    newFissXS = newFission.crossSection[styleName]          ### get the new xs as amended above with the old low energy data
    oldFission.crossSection.replace(newFissXS)              ### put the new xs in the old reaction
    RxnSte.reactions.remove(newFission.label)               ### remove the new XS-only reaction in the new reactionSuite
    RxnSte.reactions.add(oldFission)                        ### add the old full reaction in the new reactionSuite
    
    return RxnSte

def getAngGivenEin(channelName,channelData,productName,crossSection):
    
    xsMin,xsMax = crossSection.domainMin, crossSection.domainMax
    angularAxes = angularModule.defaultAxes( energyUnit )
    subform = angularModule.XYs2d( axes = angularAxes )

#     print('channelName,productName',channelName,productName)
    data = {}
    colHead = 'Prob(%s)'%productName
    datafilename = '%s_Ang_G_Ein'%''.join(channelName)
    myKey = findPandTinDictKeys(list(channelData.keys()),datafilename)
    if myKey == None: 
        print('NO ANG DIST : ', datafilename,'is not substring in any key', list(channelData.keys()))
        return None   #  No distributions
    
    D = channelData[myKey]

    if 'E_in' not in list(D.keys()) or colHead not in list(D.keys()):
        print(' No exit distributions for',channelName,'   SKip channel')
        return None   #  No distributions
#         print("From",datafilename,',',myKey,':','E_in',' not in', list(D.keys()))
#         print('D:\n',D)
    for i in range(len(D['E_in'])):
        if D['E_in'][i] not in data : data[D['E_in'][i]] = {}
        if D['cos(theta)'][i] not in data[D['E_in'][i]] : 
            data[D['E_in'][i]][D['cos(theta)'][i]] = D[colHead][i]
    

    Eins = sorted(data.keys())
    Mus = sorted(data[Eins[0]].keys())
    lastNonZeroMuP = [ [mu,.5] for mu in Mus] 
    for Ein in Eins :
        Mus = sorted(data[Ein].keys())
        muP = [ [mu,data[Ein][mu]] for mu in Mus]
        ### for zero cross section and empty distribution, use isotropic
        xsVal = crossSection.evaluate(Ein)
        if xsVal is not None and xsVal == 0.0 : 
            if sum( [ mu[1] for mu in muP] ) >= 1.e-8 : 
                lastNonZeroMuP = muP
            else: 
                muP = lastNonZeroMuP
       
        if Ein >= xsMin and Ein < xsMax :
            subform.append( angularModule.XYs1d( data = muP, outerDomainValue = Ein, axes = angularAxes ) )
        elif Ein >= xsMax :
            subform.append( angularModule.XYs1d( data = muP, outerDomainValue = xsMax, axes = angularAxes ) )
            break

    if D.get('Threshold',0.) > 0. :
        storedZero = [[mu,.5] for mu in [-1.,1.] ] 
        subform.insertAtValue(angularModule.XYs1d( data = storedZero, outerDomainValue = xsMin, axes = angularAxes ), outerDomainValue=xsMin  )

    return( subform )

def getEoutGivenAngEin( channelName, channelData, productName, crossSection ) :

    xsMin,xsMax = crossSection.domainMin, crossSection.domainMax
    data = {}
    colHead = 'Prob(%s)'%productName
    datafilename = '%s_Eout_G_Ein_Ang'%''.join(channelName)
    datafilename = '%s_Eout_Ang_G_Ein'%''.join(channelName)
    if findSubStrInDictKeys(list(channelData.keys()),'lastic_%s_Ang'%args.projectile) is not None : return( None, None, None )  ### nn' is twobody, no energy distribution necessary
    myKey = findPandTinDictKeys(list(channelData.keys()),datafilename)

    if myKey is None:
        print('Key',datafilename,'not found in',list(channelData.keys()))
        sys.exit()

    D = channelData[myKey]
    frame = standardsModule.frames.labToken
    if( D['frame'] == 'COM' ) : frame = standardsModule.frames.centerOfMassToken

    axes = energyAngularModule.defaultAxes( energyUnit )
    xys3d = energyAngularModule.XYs3d( axes = axes, interpolationQualifier = standardsModule.interpolation.unitBaseToken )
    data1 = []
    for index, E_in in enumerate( D['E_in'] ) : data1.append( [ E_in, D['E_out'][index], D['cos(theta)'][index], D[colHead][index] ] )

    E_in_list = collateByFirstElement( data1 )

    E_in_mins  = [ E_in for E_in, E_out__mu__probility in E_in_list if E_in <= xsMin ]
    if( len( E_in_mins ) == 0 ) :
        E_in_min = xsMin
        E_in_list.insert( 0, [ xsMin, [ [ 0.0, -1, 1e10 ], [ 0.0, 1.0, 1e10 ], [ 5e-11, -1, 1e10 ], [ 5e-11, 1.0, 1e10 ] ] ] )
    else :
        E_in_min = max( E_in_mins )

#     print('For',productName,':')
    for E_in, E_out__mu__probility in E_in_list :
        if( E_in < E_in_min ) : continue
        data2 = collateByFirstElement( E_out__mu__probility )
        xys2d = energyAngularModule.XYs2d( axes = axes, outerDomainValue = E_in, interpolation = standardsModule.interpolation.flatToken )
        skipLeadingZeros = True
        for E_out, mu__probility in data2 :
            xys1d = energyAngularModule.XYs1d( data = mu__probility, axes = axes, outerDomainValue = E_out, interpolation = standardsModule.interpolation.flatToken )
            if( skipLeadingZeros ) :
                if( float( xys1d.integrate( ) ) == 0.0 ) : continue
                skipLeadingZeros = False
#             if( not( skipLeadingZeros ) ) :
#                 if( float( xys1d.integrate( ) ) == 0.0 ) :
#                     if lastZero > 3:
#                         xys1d = energyAngularModule.XYs1d( data = mu__probility_prior, axes = axes, outerDomainValue = E_out, interpolation = standardsModule.interpolation.flatToken )
#                         xys2d.append( xys1d.thin( thin_eps ) )
#                         break
            xys2d.append( xys1d.thin( thin_eps ) )
            mu__probility_prior = mu__probility
        if( len( xys2d ) > 0 ) :
            integral = xys2d.integrate( )
            xys2d.normalize( insitu = True )
            integral2 = xys2d.integrate( )
            if( abs( integral - 1 ) > 1e-2 ) : 
                print( '   Normalization off by >1% in ch.', ''.join(channelName), 'giving',productName, 'at E_in=', E_in, ': %15.6f.' % integral ,'Fix:',' %15.6f.' % integral2)
#             else:
#                 print( '   Normalization ok  by <1% in ch.', ''.join(channelName), 'giving',productName, 'at E_in=', E_in, ': %15.6f.' % integral,'-> %15.6f.' % integral2)
            
            xys3d.append( xys2d )

        if( E_in >= xsMax ) : break

    for i in range(len(D['E_in'])):
        if D['E_in'][i] not in data : data[D['E_in'][i]] = {}
        if D['cos(theta)'][i] not in data[D['E_in'][i]] : 
            data[D['E_in'][i]][D['cos(theta)'][i]] = {}
        if D['E_out'][i] not in data[D['E_in'][i]][D['cos(theta)'][i]] : 
            data[D['E_in'][i]][D['cos(theta)'][i]][D['E_out'][i]] = D[colHead][i]

    angularEnergyAxes = angularEnergyModule.defaultAxes( "MeV", 'P(energy_out|energy_in,mu)' )
    subform = angularEnergyModule.XYs3d( axes = angularEnergyAxes, interpolationQualifier = standardsModule.interpolation.unitBaseToken )
    Eins = sorted(data.keys())
    for Ein in Eins :
        Mus = sorted(data[Ein].keys())
        w_xys = angularEnergyModule.XYs2d( outerDomainValue = Ein, axes = angularEnergyAxes )
        for mu in Mus  :
            Eps = sorted(data[Ein][mu].keys())
            EpP = [ [Ep,data[Ein][mu][Ep]] for Ep in Eps]
            w_xys.append( angularEnergyModule.XYs1d(  outerDomainValue = mu, data = EpP, axes = angularEnergyAxes ) )
        
        if( xsMin <= Ein < xsMax ) :
            subform.append( w_xys )
        elif Ein <= xsMin :
            storedZero = w_xys
        elif Ein >= xsMax :
            w_xys.outerDomainValue = xsMax
            subform.append( w_xys )
            break

    if D['Threshold'] > xsMin:
        storedZero = angularEnergyModule.XYs2d( outerDomainValue = D['Threshold'], axes = angularEnergyAxes )
        for mu in [-1.,1.] :  # Mus  :
            storedZero.append( angularEnergyModule.XYs1d(  outerDomainValue = mu, data = [[0,1.],[FUDGE_EPS,0.]], axes = angularEnergyAxes ) )
        subform.insertAtValue(storedZero,outerDomainValue=D['Threshold'])
    
    while subform[-1].outerDomainValue > xsMax:
        null = subform.pop()

    return( subform, xys3d, frame )


def getCrossectionFromData(channelData,channelName):
    
    crossSectionAxes = crossSectionModule.defaultAxes( energyUnit )
    for k in channelData.keys():
        if isXSfile(k):
            if 'xs' not in list(channelData[k].keys()) :  ### theres no cross section, some channels will have none, due to low sampling
                return None    ### xs was zero for this channel, and was therefore deleted from input dicts
            xsData =  [channelData[k]['E_in'],channelData[k]['xs']]
            if k.startswith('Elastic') and args.verbose>2: print('xsData:',xsData)
            crossSection = crossSectionModule.XYs1d( data = xsData, dataForm = 'xsandys', label = styleName, axes = crossSectionAxes ) ## interpolation = getXYInterpolation( I0 ), 
            crossSection = crossSection.trim()
            
            maxE = args.Emax
#             for ix,x in enumerate(crossSection):
#                 print('        xsec for',channelName,': ix,x=',ix,x[:2])

#           for ix,x in enumerate(crossSection):
#               maxE = x[0]
#               if ix > 0 and x[1] == 0.0:
#                   break
#           if maxE < args.Emax: print('        Maybe lower maxE from',args.Emax,'to',min(maxE,args.Emax))
#             maxE = min(maxE,args.Emax)
            if crossSection.domainMin > args.Emax : 
                print('        Skip channel',channelName,'starting at',crossSection.domainMin,'(threshold =', channelData[k]['Threshold'],') up to',maxE)
                return None
            
            if channelData[k]['Threshold'] > crossSection.domainMin: 
                crossSection = crossSection.domainSlice(domainMin = channelData[k]['Threshold'], domainMax = maxE )
                crossSection[0] = [channelData[k]['Threshold'],0.]
            else: 
                crossSection = crossSection.domainSlice(domainMin = crossSection[0][0], domainMax = maxE )
            crossSection.label = styleName
            return crossSection
     
def getQFromData(channelData):
    for k in list(channelData.keys()):
        if isXSfile(k):
            return channelData[k]['Q']

def isXSfile(filename):
    if ( filename.find('_cs')>=0 ) :  ### found a cross section file
        ### some cases to determine the proper file within the channel
        if ( ( filename.startswith('Channel') and filename.find('fission')==-1) or 
            filename.startswith('fission') or 
#             filename.startswith('Compound_Elastic')  or 
#             filename.startswith('Shape_Elastic')  or 
            filename.startswith('Elastic')  or 
            filename.startswith('Inelastic') ) : 
            
            return True

    return( False )
        

def findPandTinDictKeys(dictkeys,datafilename):
# modify so n_Ang_G_Ein is found in 'Channel_n2p_Ang_G_Ein_L000.dat'
    sep = datafilename.index('_') # expect 1
    product = datafilename[:sep]
    datatype = datafilename[sep+1:]
    myKeys = [ k for k in dictkeys if (datatype in k and product in k.split('_')[1]) ]
    if len(myKeys) <1 or  myKeys[0] not in dictkeys :
        return None
    else:
        return myKeys[0]

def findSubStrInDictKeys(dictkeys,datafilename):
    myKeys = [ k for k in dictkeys if (datafilename in k) ]
    if len(myKeys) <1 or  myKeys[0] not in dictkeys : 
        return None
    else:
        return myKeys[0]


def getMultiplicityFromData( channelData, productName ) :

    data = None
    k_found = None
    colHead = 'Mult(%s)' % productName
    for k in list(channelData.keys()):
        if isXSfile(k):
            if colHead in list(channelData[k].keys()) : 
                data = [channelData[k]['E_in'],channelData[k][colHead]]
                k_found = k

    if data is not None : 
        dataNoZeros = [a for a in data[1] if a != 0]
        if min(dataNoZeros)==max(dataNoZeros) : 
            data = max(dataNoZeros)
            if( int( data ) == data ) : data = 1
            if args.verbose>2: print('For',productName,'getMultiplicityFromData  data=',data, isinstance( data, ( int, float ) ))
    if args.verbose>1: 
        if data is None and productName!= 'g':
            print('   No multiplicities found for',productName,'from colHead:',colHead)
        else:
            print('   Multiplicities found for',productName,'from colHead:',colHead,'for channel',k_found)

    return data
    
def getSymAZfromIsotope(name):
    A = int(re.search('(\d+)',name).group())
    S = re.search('(\D+)',name).group()
    popsNamesModule.checkSymbol( S ) 
    Z = popsNamesModule.ZFromSymbol[S]
    ZA = 1000*Z+A
    return S,Z,A,ZA

def getPopsPartId(p):
#    print('getPopsPartId(',p,')')
    if p == 'n' : return idsPoPsModule.neutron
    if p == 'p' : return popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( 'H', 1 )
    if p == 'd' : return popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( 'H', 2 )
    if p == 't' : return popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( 'H', 3 )
    if p == 'h' : return popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( 'He', 3 )
    if p == 'a' : return popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( 'He', 4 )
    if p == 'g' : return idsPoPsModule.photon
    print('dont recognize your product : %s'%str(p))
    
def getGNDSnuclideName(p):
    print('getGNDSnuclideName(',p,')')
    GNDSnuclei = {'n':'n', 'p':'H1', 'd':'H2', 't':'H3', 'h':'He3', 'a':'He4', 'g':'photon'}
    gndsName = GNDSnuclei.get(p,None)
    if gndsName is not None: return gndsName
    print('dont recognize your product : %s'%str(p))
 


def fillPops(decayDataStructure):

    decdat = decayDataStructure
    PoPsData = databasePoPsModule.database('protare_internal', '1.0')
        

    photonPart = gaugeBosonModule.particle(idsPoPsModule.photon)
    photonPart.buildFromRawData( mass = [0,'amu'], spin = [1,'hbar'], parity = [1,''], charge = [0,'e'], halflife = ['stable','s'], label = 'default' )
    PoPsData.add(photonPart)

    neutronPart = baryonModule.particle(idsPoPsModule.neutron)
    neutronPart.buildFromRawData( mass = [1.0086649216939967,'amu'], spin = ['1/2','hbar'], parity = [1,''], charge = [0,'e'], halflife = [881.5,'s'], label = 'default' )
    PoPsData.add(neutronPart)

#     protonPart = nuclideModule.particle(getPopsPartId(idsPoPsModule.proton))
#     protonPart.buildFromRawData( mass = [1.007825029784053,'amu'], spin = ['1/2','hbar'], parity = [1,''], charge = [1,'e'], halflife = ['stable','s'], label = 'default' )
#     PoPsData.add(protonPart)

    for IsoName in list(decdat.keys()):
        chemS, chemZ, chemA, chemZA = getSymAZfromIsotope(IsoName)
        chemName = popsNamesModule.nameFromZ[chemZ]
        isoS = popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( chemS, chemA )
        
        for level in decdat[IsoName]:
            levelData = decdat[IsoName][level]
            index = level
            nuclide =  nuclideModule.particle(popsNamesModule.nuclideIDFromIsotopeSymbolAndIndex( isoS, index ))
            
            PoPsData.add(nuclide)
            if index == 0 :
                nuclide.buildFromRawData( mass = [levelData['mass'],'amu'], charge = [0,'e'] )
            else:
                nuclide.buildFromRawData( charge = [0,'e'] )
            if int(levelData['parity'])==0:
                levelData['parity'] = 1
                print("Parity 0 for",isoS,'level',index,'   Setting as ', levelData['parity'])
            nuclide.nucleus.buildFromRawData( spin = ['%d/2'%int(levelData['spin']*2),'hbar'], parity = [int(levelData['parity']),''], charge = [chemZ,'e'] )
            nuclide.nucleus.energy.add( nuclearEnergyLevelModule.double( 'default', levelData['energy'], 'MeV' ) ) 
            if 'isomerLevel' in levelData and levelData['isomerLevel'] > 0 :
                aliasName = PoPsAliasModule.metaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex( nuclide.id, int(levelData['isomerLevel']) )
                if nuclide.id not in list(aliasNameDict.keys()) : aliasNameDict[nuclide.id] = aliasName
                PoPsData.add( PoPsAliasModule.metaStable( aliasName, nuclide.id, int(levelData['isomerLevel']) ) )
                continue ### dont add decay data at the nuclide level for isomers
                
            for decayNum in range(levelData['numDecays']):
                decayData = decayDataModule.decayMode( str(decayNum), 'electroMagnetic' )
                decayData.probability.add( probabilityModule.double( 'default', levelData['decays'][decayNum]['branch'] ) )
                decPath = decayDataModule.decay( str(decayNum), '' )
                daughterState = popsNamesModule.nuclideIDFromIsotopeSymbolAndIndex( isoS, levelData['decays'][decayNum]['to'] ) 
                decPath.products.add( popsProductModule.product( daughterState, daughterState)  )
                decPath.products.add( popsProductModule.product( idsPoPsModule.photon, idsPoPsModule.photon) )
                decayData.decayPath.add( decPath )
                nuclide.decayData.decayModes.add( decayData )
                
    PoPsData.saveToFile('pops_test.xml')    
    
    return PoPsData


def getEvaluationData(): 
    
    protarePath = os.path.join(args.target,args.projectile)
    os.chdir(protarePath)
    
    if args.resave:
        args.saved = True
        if os.path.exists('reactionData.pkl')  : os.remove('reactionData.pkl')
        if os.path.exists('structureData.pkl') : os.remove('structureData.pkl')
        
    if args.saved :
        import pickle as pickle
        if os.path.exists('reactionData.pkl') and os.path.exists('structureData.pkl') :
            channels = pickle.load( open('reactionData.pkl', 'rb') )
            structureData = pickle.load( open('structureData.pkl', 'rb') )
            print('Loaded data from reactionData.pkl and structureData.pkl')
            return channels, structureData

    structureData = {}
    channels = {}
    ZAprojectile = getChannelZA(args.projectile)
    
    dirList = os.listdir('.') # list of directories
    for counter,drPath in enumerate(dirList):
        if os.path.isdir(drPath):
            if drPath == 'Spectra': continue   # for now
            channels[drPath] = {}
            targetS,targetZ,targetA,ZAtarget = getSymAZfromIsotope(args.target)
            args.target = popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( targetS, targetA )
            targetIsotopeS = popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( targetS, targetA )
            
            if drPath == 'f':
                channels[drPath]['residual'] = 'fission'
                channels[drPath]['products'] = ['n','g']
                ### read in fission file here then we will continue with the next channel (don't need to calculate all the ZA stuff for this channel)
                for chan in os.listdir(drPath) : 
                    if chan == 'Fission-Stats.dat': continue
                    fileData = parseDataFile(os.path.join(drPath,chan))
                    fileData['channel'] = drPath
                    
                    apparentThresh = min([x[0] for x in zip(fileData['E_in'],fileData['xs']) if x[1] > 0.0 ])
                    if fileData['E_in'][0] == apparentThresh : 
                        fileData['Threshold'] = 0.
                    else:
                        fileData['Threshold'] = apparentThresh
                        
                    fileData['Q'] = 180.
                    channels[drPath][chan] = fileData
                continue
            


            ZAproducts = getChannelZA(drPath)
            ZAresidual = ZAprojectile + ZAtarget - ZAproducts 
            residualS = popsNamesModule.symbolFromZ[ZAresidual//1000]
            residualA = ZAresidual%1000
            residualIsotopeS = popsNamesModule.isotopeSymbolFromChemicalElementIDAndA( residualS, residualA )
            products = list(set(list(drPath)))
                
            ### move residual and product ZA info to each individual file    
            chanList = os.listdir(drPath)
            for chan in chanList :
                if chan[0] == '.' or chan == 'gammas' : continue                     # Can happen when someone is looking at a file with an editor.
                ### FIXME  Erich needs to fix this in his file names.
                if 'cs_gammas' in chan or chan.startswith('Compound_Elastic') or chan.startswith('Shape_Elastic') or 'Leg' in chan:
                    if args.verbose>2: print('File',chan,'ignored')
                    continue
                else:
                    if args.verbose>1: print('\nRead',chan,chan.startswith('Compound_Elastic'), chan.startswith('Shape_Elastic'))
                chanStr = chan.replace('lastic_Ang','lastic_n_Ang_G_Ein').replace('lastic_cs','lastic_n_cs') 
                chanStr = chanStr.replace('tic_n_','tic_%s_'% args.projectile)
                if args.verbose and chan != chanStr: print('\nChange file name',chan,'to internal:',chanStr)
                level = 0
                if re.search('L(\d+).dat',chan) is not None: 
                    level = int(re.search('L(\d+).dat',chan).groups()[0])
                
                if level > 0 :          ### direct inelastic neutron out
                    chanEx = '%s_L%03d'%(drPath,level)
                    myiso = popsNamesModule.nuclideIDFromIsotopeSymbolAndIndex( residualIsotopeS, level )
                    if drPath==args.projectile and 'Channel_' in chanStr : 
                        chanEx = '%s_g_L%03d'%(args.projectile,level)            ### Inelastic continuum to a metastable      
                    myproducts = products + ['g']
                else:                   ### all other channels
                    chanEx = drPath
                    myproducts = products
                    myiso = residualIsotopeS
                    if drPath==args.projectile:         ### for neutron out lets differentiate elastic from inelastic continuum
                        if 'Channel_' in chanStr :      ### Inelastic continuum
                            chanEx = '%s_g' % args.projectile
                            myproducts = products + ['g']
                        elif 'Elastic' in chanStr:      ### Elastic
                            chanEx = args.projectile
                            myproducts = products 
#                             if args.verbose:  print("Elastic chanEx,products=",args.projectile,myproducts,chanEx in channels)
                   
                if chanEx not in channels: channels[chanEx] = {}
#                 print('Channels residual for',drPath,':',myproducts,'is',myiso,ZAresidual)
                channels[chanEx]['products'] = myproducts
                channels[chanEx]['residual'] = myiso
                channels[chanEx]['residualZA'] = ZAresidual
                channels[chanEx]['productZA'] = ZAproducts
                channels[chanEx][chanStr] = parseDataFile( os.path.join( drPath, chan ) )
                channels[chanEx][chanStr]['channel'] = chanEx
        else:
            if '-' in drPath and drPath.split('-')[1] == 'decay' :
                structureData[drPath.split('-')[0]] = parseDecayFile(drPath)
#                 print('structureData for',drPath.split('-')[0],'from',drPath)
            elif drPath == 'Reaction.dat' : 
                xsReaction = parseDataFile(drPath)
            elif drPath == 'Total_cs.dat' : 
                xsTotal = parseDataFile(drPath)
   
    for chanEx in list(channels.keys()) :
        drPath = chanEx.split('_')[0]
        basePath = drPath
        if chanEx.startswith('%s_'%args.projectile) or chanEx==args.projectile : basePath = '%s_g'%args.projectile
        for chanStr in list(channels[chanEx].keys()) :
            if 'fission' not in chanStr :
                fileData = channels[chanEx][chanStr]
                if isinstance(fileData,dict):
                    baseStr = findSubStrInDictKeys(list(channels[basePath].keys()),'Channel_%s_cs_'%drPath)
                    if baseStr is None:
                        print('********* basePath:',basePath,'baseStr:',baseStr)
                        continue
                    if args.verbose > 1: print('  Check masses for',chanEx,'from',basePath,baseStr)
                    fileData = checkMasses(fileData,channels[basePath][baseStr])
                    channels[chanEx][chanStr] = fileData
                    if args.verbose: print('Ch',chanEx,'->',chanStr)
  
    if args.saved :
        pickle.dump(channels, open('reactionData.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)              ## this protocol should be faster than the python native pickle
        pickle.dump(structureData, open('structureData.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        print('Saved data to reactionData.pkl and structureData.pkl')
#     print('getEvaluationData - structureData=',structureData)  
    return channels, structureData


def checkMasses(fileData,channelsData) : 
    
       
    fileData['numPartOut'] = channelsData['numPartOut']
    strings = ['massProjectile','massTarget','massResidual','massAmu']+['massPart%d'%(c+1) for c in range(channelsData['numPartOut'])]
    for S in strings: 
        fileData[S] = channelsData[S]
        if args.verbose > 2: print(' The mass of', S,'is',fileData[S])
    if args.verbose > 2: print('Mass residual=',fileData['massResidual'])
    if 'excitation' not in fileData: fileData['excitation'] = 0.0
    
    ### check Q value (is there one?)
    ### compute from Masses and excitation level
    ### Mtarget+Mprojectile - Mresidual + excitation - Memitted = Q
    mysum = 0.
    for c in range(fileData['numPartOut']) :
        mysum += fileData['massPart%d'%(c+1)]
    
    massDiffAmu = (fileData['massTarget']+fileData['massProjectile'] - fileData['massResidual'] - mysum )
    massDiffMev = massDiffAmu*fileData['massAmu']
    calcQ =  massDiffMev - fileData['excitation']
    if 'Q' in fileData: 
        parsedQ = fileData['Q']
        if abs(parsedQ) > 0 and  (calcQ - parsedQ)/parsedQ > 1e-2 : 
            print('Calc Q (%g) doesnt agree with quoted Q (%g) '%(calcQ,parsedQ),'so mismatch =',calcQ-parsedQ,' E*=',fileData['excitation'])
            print('    Using ',fileData['massTarget'],'+',fileData['massProjectile'],'-', fileData['massResidual'] ,'-',mysum,'amu +',fileData['excitation'],'MeV')
        if abs(parsedQ + calcQ) <= 1e-4  and abs(calcQ)>1e-3: 
            if args.verbose: print('   but quoted Q (%g) and calculated Q (%g) only differ by a minus sign'%(parsedQ,calcQ))
            fileData['Q'] = -parsedQ
    else:
        if 'excitation' in fileData: 
            fileData['Q'] = -fileData['excitation']
        else:
            print('replacing missing Q with calculated : %f '%calcQ)
            fileData['Q']=calcQ
    
    Threshold = -fileData['Q']*(1. + fileData['massProjectile']/fileData['massTarget'])
    fileData['Threshold'] = float('%12.8f'%max(Threshold,0.0))
    
    return   fileData

def getChannelZA(chanstr):

    prodZAs = {'n':1,'p':1001,'d':1002,'t':1003,'h':2003,'a':2004,'g':0}
    particle = ''
    multiplicity = 0
#     particleMultiplicity = []
    totZA = 0
    for d in reversed(list(chanstr)):
        if d.isalpha():
            if multiplicity > 0: # particleMultiplicity.append(particle, multiplicity
                totZA += multiplicity * prodZAs[particle]
            particle = d
            multiplicity = 1
        else:  # digit
            multiplicity = int(d)
    totZA += multiplicity * prodZAs[particle]            

    return totZA     
 
def parseDecayFile(filename) :
    with open(filename,'r') as F:
        allLines = F.readlines()
    fileMass = 0.
    isomerLevel = 0.
    levels = {}
    Lines = iter(allLines)
    for L in Lines:
        if not L.strip().startswith('#') and not L.strip().startswith('--') :
            lvl = L.split()
            levNum = int(lvl[0])
            if levNum not in levels:
                levels[levNum] = {}
                levels[levNum]['energy'] = float(lvl[1])
                levels[levNum]['mass'] = fileMass
                levels[levNum]['spin']   = float(lvl[2])
                levels[levNum]['parity'] = float(lvl[3])
                levels[levNum]['numDecays'] = int(lvl[4])
                if levels[levNum]['numDecays'] > 0 : levels[levNum]['decays'] = []
                levels[levNum]['isomer'] = str(lvl[7])
                if levels[levNum]['isomer'] == 'True' :
                    isomerLevel += 1
                    levels[levNum]['isomerLevel'] = isomerLevel
            elif L.find('--->') > 0 :
                levels[levNum]['decays'].append({'from':int(lvl[0]), 'to':int(lvl[2]), 'branch':float(lvl[3]) , 'gammaProb' : float(lvl[4]) ,  'icProb': float(lvl[5])}) 
        else:
            if L.strip('#').startswith('Mass') and L.strip().endswith('amu'):
                fileMass = float(L.split()[2])

    return levels
    
def parseDataFile(filename) :
    channel,filestr = os.path.split(filename)
    with open(filename,'r') as F:
        Lines = F.readlines()
    if args.verbose>0: print('Input file',filename,'for',channel)

    headerLines = []
    allPlotDict = {}
    ### get the header lines
    for L in Lines:
        if L.startswith('#'):
            if re.match('# +E_in =',L): 
                continue
            headerLines.append(L.strip('#').strip())
        else:
            break ## header is over

    ### parse the header lines
    outPartCnt = 0
    for L in headerLines:
        Ls = L.split()
        if 'level' in L: 
            try : 
                allPlotDict['level'] = int(Ls[-1])
            except : 
                pass
                #allPlotDict['level'] = Ls[-1]
        if 'Q-value' in L: allPlotDict['Q'] = float(Ls[-2])
        if 'Final state' in L and 'Ex =' in L: 
            allPlotDict['excitation'] = float(Ls[-2])
            allPlotDict['energyUnit'] = Ls[-1]
        if 'Mass amu' in L: 
            allPlotDict['massAmu'] = float('%20.12f'%float(Ls[-2]))
            allPlotDict['massUnit'] = 'amu'
        if 'Mass of target' in L: 
            allPlotDict['massTarget'] = float('%20.12f'%float(Ls[-2]))
        if 'Mass of projectile' in L: 
            allPlotDict['massProjectile'] = float('%20.12f'%float(Ls[-2]))
        if 'Mass of residual' in L: 
            allPlotDict['massResidual'] = float('%20.12f'%float(Ls[-2]))
        if 'Mass particle' in L: 
            allPlotDict['massPart%d'%int(Ls[3])] = float('%20.12f'%float(Ls[-2]))
            outPartCnt = int(Ls[3])
        if 'Frame' in L: allPlotDict['frame'] = Ls[-1]
        if 'E_in' in L: 
            ### some column header cleanups
            ### FIXME  these should all be fixed in the format outputs
            L = L.replace('Eout','E_out') 
            L = L.replace('<EL/Rutherford>','xs')
            L = L.replace('Prob','Prob(%s)' % args.projectile).replace('Prob(%s)('% args.projectile,'Prob(') 
            L = L.replace('EL/Rutherford','Prob(%s)' % args.projectile).replace('Prob(%s)('% args.projectile,'Prob(') 
            L = L.replace('( ','(')
            L = L.replace('xs(b)','xs')
            L = L.replace('xs( b)','xs')
            allPlotDict['columns'] = L.split()
            for k in allPlotDict['columns']:
                allPlotDict[k] = []
    
    allPlotDict['numPartOut'] = outPartCnt
    if args.verbose > 2: 
        try:
            print('p,T,R,e1=',allPlotDict['massProjectile'],allPlotDict['massTarget'],allPlotDict['massResidual'], allPlotDict['massPart1'])
        except:
            pass
            
    ### get data lines
    for L in Lines:
        if not L.startswith('#'):
            nums = []
            for nc in L.split():
                nn = nc if nc[-1] not in ['+','-'] else nc[:-1]
                nums.append(float(nn.replace('E','e')))
            for index,num in enumerate(nums):
                allPlotDict[allPlotDict['columns'][index]].append(num) 

    if 'columns' in allPlotDict.keys():
        for index,label in enumerate(allPlotDict['columns']):
            if  len(allPlotDict[allPlotDict['columns'][index]]) > 0 and max(allPlotDict[allPlotDict['columns'][index]]) == 0.0 : 
                del allPlotDict[allPlotDict['columns'][index]]
    
    ### FIXME : Erich should add these for the parser's sake
    ### add in (n,n) multiplicity if it wasnt stated 
    if channel == args.projectile and 'Mult(%s)'%args.projectile not in allPlotDict : 
        if 'Mult(%s)'%args.projectile not in allPlotDict : 
            allPlotDict['Mult(%s)'%args.projectile] = [1.]

    multis = {}
    for product in channel:  # string
        if product != 'g':
            mult = 'Mult(%s)'%product
            if mult not in allPlotDict : multis[product] = 1 # if product not in multis.keys() else multis[product]+1
    for product in multis.keys():  # string
        allPlotDict['Mult(%s)'%product] = [multis[product] ]
            
    return allPlotDict

def defineInitialNTs(projectile):

    NTdict = {projectile: 2,   'f' : 18}  # Elastic & fission channels

    for i in range(1,500):  # do both gs and excited states
    # discrete levels
        NTdict['n_L%03d'%i] = 1000+i
        NTdict['p_L%03d'%i] = 10000+i 
        NTdict['d_L%03d'%i] = 100000+i
        NTdict['t_L%03d'%i] = 1000000+i
        NTdict['h_L%03d'%i] = 10000000+i
        NTdict['a_L%03d'%i] = 100000000+i
        NTdict['g_L%03d'%i] = 1000000000+i   

    for j in range(1,100):  
    # continuum (unresolved) spectra leaving excited states
        i = j+900
        NTdict['n_g_L%03d'%i] = 1000+i
        NTdict['p_g_L%03d'%i] = 10000+i 
        NTdict['d_g_L%03d'%i] = 100000+i
        NTdict['t_g_L%03d'%i] = 1000000+i
        NTdict['h_g_L%03d'%i] = 10000000+i
        NTdict['a_g_L%03d'%i] = 100000000+i  
    
    return(NTdict)

def niceSortOfNTs( NTs, verbose = 0, logFile = sys.stderr ) :
    import copy

    def removeGetIfPresent( NT, NTs ) :

        if( NT not in NTs ) : return( [] )
        NTs.remove( NT )
        return( [ NT ] )

    NTs = sorted( list( copy.deepcopy( NTs ) ) )
    
    newNTs = []
    newNTs += removeGetIfPresent(   1, NTs )
    newNTs += removeGetIfPresent(   2, NTs )
    newNTs += removeGetIfPresent(   4, NTs )

    newNTs += removeGetIfPresent(  18, NTs )   # (z,f)
    newNTs += removeGetIfPresent(  19, NTs )   # (z,0nf)
    newNTs += removeGetIfPresent(  20, NTs )   # (z,1nf)
    newNTs += removeGetIfPresent(  21, NTs )   # (z,2nf)
    newNTs += removeGetIfPresent(  38, NTs )   # (z,3nf)

    for NT in sorted(NTs):  # try simplest numerical method that does not have to make all possible integers!!   # also use copy of initial NTs list
        if NT >= 1000:   # deal with all new NT numbers, leaving ENDF_MTs
            newNTs += removeGetIfPresent( NT, NTs )        

    for NT in sorted(NTs):  # all the other NT values, ie. those added as NT values.
        if NT not in [5,10, 27, 101, 151] and not (500<=NT<573)  and not (201<=NT<600)  and not (850<=NT<875) : # skip those dealt with just below
            newNTs += removeGetIfPresent( NT, NTs )

    NT5 = removeGetIfPresent( 5, NTs )                 # (z,everything else)

    NTAtomics = []
    for NT in range( 500, 573 ) : NTAtomics += removeGetIfPresent( NT, NTs )

    skippingNTs = []
    for NT in [ 10, 27, 101, 151 ] : skippingNTs += removeGetIfPresent( NT, NTs )
    for NT in range( 201, 600 ) : skippingNTs += removeGetIfPresent( NT, NTs )

    if( ( verbose > 0 ) and ( len( skippingNTs ) > 0 ) ) : logFile.write( 'Skipping NTs = %s\n' % skippingNTs )

    newNTs += NTs + NT5 + NTAtomics
    return( newNTs )

if __name__ == "__main__":
  main()
