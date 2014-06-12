import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso

process = cms.Process("p")

MC_flag = False

HLTProcessName = "HLT"
OUTPUT_FILE_NAME = "TnPTree_Data.root"


ELECTRON_ET_CUT_MIN = 10.0
ELECTRON_COLL = "gsfElectrons"
ELECTRON_CUTS = "(abs(superCluster.eta)<2.5) && pt > 10.0"
#ELECTRON_CUTS = "(abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"



#########################################################
#NEUTRAL TRIGGER
#########################################################
HLTEle17Mass50Path1 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v3"
HLTEle17Mass50Path2 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v4"
HLTEle17Mass50Path3 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v5"
HLTEle17Mass50Path4 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v6"

HLTEle8Mass50Path1 = HLTEle17Mass50Path1
HLTEle8Mass50Path2 = HLTEle17Mass50Path2
HLTEle8Mass50Path3 = HLTEle17Mass50Path3
HLTEle8Mass50Path4 = HLTEle17Mass50Path4

nuetHLTEle17Mass50Filter = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter"
#nuetHLTEle17Mass50Filter = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8DoubleEtFilter" main check kita c ki is naal koi entry nai bachdi
nuetHLTEle8Mass50Filter = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter"

HLTEle20Mass50Path1 = "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v3"
HLTEle20Mass50Path2 = "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v4"
HLTEle20Mass50Path3 = "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v5"
HLTEle20Mass50Path4 = "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v6"
HLTEle20Mass50Path5 = "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v7"

HLTSC4Mass50Path1 = HLTEle20Mass50Path1
HLTSC4Mass50Path2 = HLTEle20Mass50Path2
HLTSC4Mass50Path3 = HLTEle20Mass50Path3
HLTSC4Mass50Path4 = HLTEle20Mass50Path4
HLTSC4Mass50Path5 = HLTEle20Mass50Path5

nuetHLTEle20Mass50Filter = "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter"
nuetHLTSC4Mass50Filter = "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter"

hltTagsPassingEle8Mass50HLTGsf = cms.VInputTag(
     cms.InputTag(HLTEle8Mass50Path1, nuetHLTEle8Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle8Mass50Path2, nuetHLTEle8Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle8Mass50Path3, nuetHLTEle8Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle8Mass50Path4, nuetHLTEle8Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path1, nuetHLTSC4Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path2, nuetHLTSC4Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path3, nuetHLTSC4Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path4, nuetHLTSC4Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path5, nuetHLTSC4Mass50Filter, HLTProcessName),

     )

hltTagsPassingEle17Mass50HLT =cms.VInputTag(
     cms.InputTag(HLTEle17Mass50Path1, nuetHLTEle17Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle17Mass50Path2, nuetHLTEle17Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle17Mass50Path3, nuetHLTEle17Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle17Mass50Path4, nuetHLTEle17Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle20Mass50Path1, nuetHLTEle20Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle20Mass50Path2, nuetHLTEle20Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle20Mass50Path3, nuetHLTEle20Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle20Mass50Path4, nuetHLTEle20Mass50Filter, HLTProcessName),
     cms.InputTag(HLTEle20Mass50Path5, nuetHLTEle20Mass50Filter, HLTProcessName),

     )


#######################################################################
#TRIGGER TO BE TESTED                                                 #
#######################################################################
HLTEle17Path1 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15"
HLTEle17Path2 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16"
HLTEle17Path3 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17"
HLTEle17Path4 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18"
HLTEle17Path5 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19"

hltEle17Filter = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter"
hltEle8Filter = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter"#hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ"

hltTagsPassingEle17HLT= cms.VInputTag(
    cms.InputTag(HLTEle17Path1,hltEle17Filter, HLTProcessName),
    cms.InputTag(HLTEle17Path2,hltEle17Filter, HLTProcessName),
    cms.InputTag(HLTEle17Path3,hltEle17Filter, HLTProcessName),
    cms.InputTag(HLTEle17Path4,hltEle17Filter, HLTProcessName),
    cms.InputTag(HLTEle17Path5,hltEle17Filter, HLTProcessName),
     )
hltTagsPassingEle8HLT =cms.VInputTag(
    cms.InputTag(HLTEle17Path1,hltEle8Filter, HLTProcessName),
    cms.InputTag(HLTEle17Path2,hltEle8Filter, HLTProcessName),
    cms.InputTag(HLTEle17Path3,hltEle8Filter, HLTProcessName),
    cms.InputTag(HLTEle17Path4,hltEle8Filter, HLTProcessName),
    cms.InputTag(HLTEle17Path5,hltEle8Filter, HLTProcessName),
     )



########################

##    ___            _           _      
##   |_ _|_ __   ___| |_   _  __| | ___ 
##    | || '_ \ / __| | | | |/ _` |/ _ \
##    | || | | | (__| | |_| | (_| |  __/
##   |___|_| |_|\___|_|\__,_|\__,_|\___|
##
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'FT_53_V21_AN3::All'
##process.GlobalTag.globaltag = GLOBAL_TAG
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#      $inputFileNames
#'dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/data/Run2012C/DoubleElectron/AOD/22Jan2013-v1/20002/78704EA6-AD68-E211-9272-0030486792AC.root'
'dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/data/Run2012B/DoubleElectron/AOD/22Jan2013-v1/20000/025B5E30-3569-E211-B1E3-003048678ED2.root'
#'file:/hdfs/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/FEF61011-7E90-E211-A84B-003048FFD75C.root'
#'/hdfs/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/20000/FEAF0DFD-9667-E211-9A4E-0025905938A8.root'
                                                  ),
eventsToProcess = cms.untracked.VEventRange(
#'193557:8215621',
'194108:484540317')
)

#if not MC_flag:
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
myLumis = LumiList.LumiList(filename = '/afs/hep.wisc.edu/home/kaur/TnP/CMSSW_5_3_9/src/PhysicsTools/TagAndProbe/test/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt').getCMSSWString().split(',')
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    
process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")


#process.load('RecoJets.JetProducers.kt4PFJets_cfi')
#process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
#process.JetsForRho = cms.Sequence( process.kt6PFJets )

## ==== Fast Filters ====
process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
    filter = cms.bool(True),
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)
process.fastFilter = cms.Sequence(process.goodVertexFilter + process.noScraping)


##    ____      __ _____ _           _                   
##   / ___|___ / _| ____| | ___  ___| |_ _ __ ___  _ __  
##  | |  _/ __| |_|  _| | |/ _ \/ __| __| '__/ _ \| '_ \ 
##  | |_| \__ \  _| |___| |  __/ (__| |_| | | (_) | | | |
##   \____|___/_| |_____|_|\___|\___|\__|_|  \___/|_| |_|
##  
#  GsfElectron ################ 

process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')


process.goodElectrons = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( ELECTRON_COLL ),
    cut = cms.string( ELECTRON_CUTS )    
)
process.goodElectronsTAG = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( ELECTRON_COLL ),
    cut = cms.string( "(abs(superCluster.eta)<=2.5) && !(1.4442<=abs(superCluster.eta)<=1.566) && pt >= 25.0" )
)
process.goodElectronsTAG20 = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( ELECTRON_COLL ),
    cut = cms.string( "(abs(superCluster.eta)<2.5) && pt > 20.0" )
)


process.PassingEle8Mass50HLTGsf = cms.EDProducer("trgMatchedPatElectronProducer",    
    InputProducer = cms.InputTag("goodElectrons" ),
    isTriggerFilter = cms.untracked.bool(True),
#    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle8Mass50HLTGsf,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

           

##    _____ _           _                     ___    _ 
##   | ____| | ___  ___| |_ _ __ ___  _ __   |_ _|__| |
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \   | |/ _` |
##   | |___| |  __/ (__| |_| | | (_) | | | |  | | (_| |
##   |_____|_|\___|\___|\__|_|  \___/|_| |_| |___\__,_|
##   
# Electron ID  ######

process.PassingWPMedium = cms.EDProducer("ElectronCutBasedCandidateProducer",
   ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
   IsoDepElectron = cms.VInputTag(cms.InputTag('elPFIsoDepositChargedPFIso'),
                                  cms.InputTag('elPFIsoDepositGammaPFIso'),
                                  cms.InputTag('elPFIsoDepositNeutralPFIso')),
   IsoValElectronPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
   conversionsInputTag     = cms.InputTag("allConversions"),
   beamSpotInputTag        = cms.InputTag("offlineBeamSpot"),
   rhoIsoInputTag          = cms.InputTag("kt6PFJets", "rho"),
   primaryVertexInputTag   = cms.InputTag("goodVertexFilter"),#offlinePrimaryVertices"),
   nameIDBOOL = cms.string('medium'),
)
#///meri bari

process.PassingWPTight = process.PassingWPMedium.clone()
process.PassingWPTight.nameIDBOOL = cms.string('tight')

process.PassingWPLoose = process.PassingWPMedium.clone()
process.PassingWPLoose.nameIDBOOL = cms.string('loose')

process.PassingWPVeto = process.PassingWPMedium.clone()
process.PassingWPVeto.nameIDBOOL = cms.string('veto')

process.PassingWPMediumTAG20 = process.PassingWPMedium.clone()
process.PassingWPMediumTAG20.ReferenceElectronCollection = cms.untracked.InputTag("goodElectronsTAG20")

process.PassingWPMediumTAG = process.PassingWPMedium.clone()
process.PassingWPMediumTAG.ReferenceElectronCollection = cms.untracked.InputTag("goodElectronsTAG")

process.PassingWPTightTAG = process.PassingWPTight.clone()
process.PassingWPTightTAG.ReferenceElectronCollection = cms.untracked.InputTag("goodElectronsTAG")

process.PassingWPTightTAG20 = process.PassingWPTight.clone()
process.PassingWPTightTAG20.ReferenceElectronCollection = cms.untracked.InputTag("goodElectronsTAG20")

#///ithe tak 
##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
##   
# Trigger  ##################
process.PassingEle17Mass50HLTTightTAG = cms.EDProducer("trgMatchedGsfElectronProducer",  
    InputProducer = cms.InputTag("PassingWPTightTAG" ),
    isTriggerFilter = cms.untracked.bool(True),
    #matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle17Mass50HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)
)


process.ele_sequence = cms.Sequence(
    process.pfParticleSelectionSequence +
    process.eleIsoSequence +
    process.goodElectrons +
    process.goodElectronsTAG +
    process.PassingEle8Mass50HLTGsf +
    process.PassingWPMedium +
    process.PassingWPTightTAG +
    process.PassingWPLoose +
    process.PassingWPVeto +
    process.PassingWPTight+
    process.PassingEle17Mass50HLTTightTAG

    )



##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
#  Tag & probe selection ######

process.tagTightGsf = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingEle17Mass50HLTTightTAG@+ goodElectrons@-"), # charge conjugate states are implied
    checkCharge = cms.bool(True),
    cut   = cms.string("50 < mass < 150"),
)

process.allTagsAndProbes = cms.Sequence(
    process.tagTightGsf
)
############################################################################
##    _____           _       _ ____            _            _   _  ____  ##
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |/ ___| ##
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | |  _  ##
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | |_| | ##
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_|\____| ##
##              |___/                                                     ##
##                                                                        ##
############################################################################
##    ____                      _     _           
##   |  _ \ ___ _   _ ___  __ _| |__ | | ___  ___ 
##   | |_) / _ \ | | / __|/ _` | '_ \| |/ _ \/ __|
##   |  _ <  __/ |_| \__ \ (_| | |_) | |  __/\__ \
##   |_| \_\___|\__,_|___/\__,_|_.__/|_|\___||___/
##
## I define some common variables for re-use later.
## This will save us repeating the same code for each efficiency category
ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt  = cms.string("pt"),
    mass  = cms.string("mass"),
)   

ProbeVariablesToStore = cms.PSet(
    probe_gsfEle_eta = cms.string("eta"),
    probe_gsfEle_abseta = cms.string("abs(eta)"),
    probe_gsfEle_pt  = cms.string("pt"),
    probe_gsfEle_et  = cms.string("et"),
    probe_gsfEle_e  = cms.string("energy"),
    probe_gsfEle_q  = cms.string("charge"),
    probe_gsfEle_trackiso = cms.string("dr03TkSumPt"),
    probe_gsfEle_reltrackiso = cms.string("dr03TkSumPt/pt"),
## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),  
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_abseta    = cms.string("abs(superCluster.eta)"),

#id based
    probe_gsfEle_dEtaIn = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    probe_gsfEle_dPhiIn  = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    probe_gsfEle_sigmaIEtaIEta = cms.string("sigmaIetaIeta"),
    probe_gsfEle_hoe           = cms.string("hadronicOverEm"),
    probe_gsfEle_ooemoop       = cms.string("(1.0/ecalEnergy - eSuperClusterOverP/ecalEnergy)"),
    probe_gsfEle_mHits = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits")
)


TagVariablesToStore = cms.PSet(
    gsfEle_eta = cms.string("eta"),
    gsfEle_abseta = cms.string("abs(eta)"),
    gsfEle_pt  = cms.string("pt"),
    gsfEle_et  = cms.string("et"),
    gsfEle_e  = cms.string("energy"),
    gsfEle_q  = cms.string("charge"),
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    sc_eta    = cms.string("superCluster.eta"),
    sc_abseta    = cms.string("abs(superCluster.eta)"),
)

CommonStuffForGsfElectronProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStore),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool (True),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
          mass60to120 = cms.string("60 < mass < 120")
    ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags     =  cms.PSet(
          passingGsf = cms.InputTag("goodElectrons"),
          isWPMedium = cms.InputTag("PassingWPMedium"),
    ),    
)

CommonStuffForSuperClusterProbe = CommonStuffForGsfElectronProbe.clone()
CommonStuffForSuperClusterProbe.variables = cms.PSet(
    probe_eta = cms.string("eta"),
    probe_abseta = cms.string("abs(eta)"),
    probe_pt  = cms.string("pt"),
    probe_et  = cms.string("et"),
    probe_e  = cms.string("energy"),
#    probe_trackiso = cms.string("dr03TkSumPt"),
 #   probe_reltrackiso = cms.string("dr03TkSumPt/pt"),
#    probe_HoverE         = cms.string("hcalOverEcal"),
#    probe_sigmaIetaIeta = cms.string("sigmaIetaIeta"),
    )


if MC_flag:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(True),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_flag),
        checkMotherInUnbiasEff = cms.bool(MC_flag),
        mcVariables = cms.PSet(
        probe_eta = cms.string("eta"),
        probe_abseta = cms.string("abs(eta)"),
        probe_pt  = cms.string("pt"),
        probe_et  = cms.string("et"),
        probe_e  = cms.string("energy"),
        probe_mass  = cms.string("mass"),
        ),
        mcFlags     =  cms.PSet(
        probe_flag = cms.string("pt>0")
        ),      
        )
else:
     mcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )

##   ____      __       __    ___                 ___    _ 
##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
##                                           |/            
##  gsf electron --> isolation, electron id  etc.
process.GsfElectronToIdMedium = cms.EDAnalyzer("TagProbeFitTreeProducer",
#    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("tagTightGsf"),
    arbitration   = cms.string("None"),
    flags = cms.PSet(
#        probe_isWP95 = cms.InputTag("WP95MatchedSuperClusterCandsClean"),
        probe_isWPMedium = cms.InputTag("PassingWPMedium"),
    ),
#    probeMatches  = cms.InputTag("McMatchGsf"),
#    allProbes     = cms.InputTag("PassingEle8Mass50HLTGsf")
    allProbes     = cms.InputTag("goodElectrons")
)

process.ElecPt = cms.EDProducer("ProbeVarMap",
    probes = cms.InputTag("gsfElectrons"),
   IsoDepElectron = cms.VInputTag(cms.InputTag('elPFIsoDepositChargedPFIso'),
                                  cms.InputTag('elPFIsoDepositGammaPFIso'),
                                  cms.InputTag('elPFIsoDepositNeutralPFIso')),
   IsoValElectronPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
   rhoIsoInputTag          = cms.InputTag("kt6PFJets", "rho"),

)
process.GsfElectronToIdMedium.variables.probe_d0vtx = cms.InputTag("ElecPt","getd0vtx")
process.GsfElectronToIdMedium.variables.probe_dzvtx = cms.InputTag("ElecPt","getdzvtx")
process.GsfElectronToIdMedium.variables.probe_iso = cms.InputTag("ElecPt","getiso")
process.GsfElectronToIdMedium.variables.probe_iso_nh = cms.InputTag("ElecPt","getisonh")
process.GsfElectronToIdMedium.variables.probe_iso_em = cms.InputTag("ElecPt","getisoem")
process.GsfElectronToIdMedium.variables.probe_iso_ch = cms.InputTag("ElecPt","getisoch")
process.GsfElectronToIdMedium.variables.probe_vtxFitConversion=cms.InputTag("ElecPt","getvtxFitConversion")

process.GsfElectronToIdTight = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdTight.flags = cms.PSet(probe_isWPTight = cms.InputTag("PassingWPTight"))

process.GsfElectronToIdLoose = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdLoose.flags = cms.PSet(probe_isWPLoose = cms.InputTag("PassingWPLoose"))

process.GsfElectronToIdVeto = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdVeto.flags = cms.PSet(probe_isWPVeto = cms.InputTag("PassingWPVeto"))

process.tree_sequence = cms.Sequence(
    process.GsfElectronToIdMedium 
#    process.GsfElectronToIdTight +
#    process.GsfElectronToIdLoose +
#    process.GsfElectronToIdVeto 
    )    

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##
process.out = cms.OutputModule("PoolOutputModule", 
       fileName = cms.untracked.string("PFEE.root"),
       SelectEvents = cms.untracked.PSet( 
       SelectEvents = cms.vstring("p")
       )
    )
process.outpath = cms.EndPath(process.out)
process.outpath.remove(process.out)

process.p = cms.Path(
        process.fastFilter+
        process.ele_sequence + 
        process.allTagsAndProbes +
           process.ElecPt+
    process.tree_sequence
        )
    
process.TFileService = cms.Service(
    "TFileService", fileName = cms.string( OUTPUT_FILE_NAME )
#    "TFileService", fileName = cms.string( '$outputFileName') 
)
