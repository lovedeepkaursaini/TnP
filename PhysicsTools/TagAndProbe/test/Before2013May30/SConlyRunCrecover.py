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

SUPERCLUSTER_COLL_EB = "correctedHybridSuperClusters"
SUPERCLUSTER_COLL_EE = "correctedMulti5x5SuperClustersWithPreshower"
SUPERCLUSTER_CUTS = "abs(eta)<2.5 && et>" + str(ELECTRON_ET_CUT_MIN)


#########################################################
#NEUTRAL TRIGGER
#########################################################

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
     cms.InputTag(HLTSC4Mass50Path1, nuetHLTSC4Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path2, nuetHLTSC4Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path3, nuetHLTSC4Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path4, nuetHLTSC4Mass50Filter, HLTProcessName),
     cms.InputTag(HLTSC4Mass50Path5, nuetHLTSC4Mass50Filter, HLTProcessName),
     )

hltTagsPassingEle17Mass50HLT =cms.VInputTag(
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
process.GlobalTag.globaltag = 'GR_R_53_V14::All'
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
      $inputFileNames
#'file:/hdfs/store/data/Run2012A/DoubleElectron/AOD/13Jul2012-v1/00000/FE1CBE0B-03DA-E111-A2AD-00266CFAE20C.root'
                                                  ),
#eventsToProcess = cms.untracked.VEventRange(
#'193557:8215621',
#'193557:35284992')
)

#if not MC_flag:
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
myLumis = LumiList.LumiList(filename = '/afs/hep.wisc.edu/home/anil79/TnP/CMSSW_5_3_5/src/PhysicsTools/TagAndProbe/test/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt').getCMSSWString().split(',')
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
##   ____                         ____ _           _            
##  / ___| _   _ _ __   ___ _ __ / ___| |_   _ ___| |_ ___ _ __ 
##  \___ \| | | | '_ \ / _ \ '__| |   | | | | / __| __/ _ \ '__|
##   ___) | |_| | |_) |  __/ |  | |___| | |_| \__ \ ||  __/ |   
##  |____/ \__,_| .__/ \___|_|   \____|_|\__,_|___/\__\___|_|   
##  

#  SuperClusters  ################
process.superClusters = cms.EDProducer("SuperClusterMerger",
   src = cms.VInputTag(cms.InputTag( SUPERCLUSTER_COLL_EB ,"", "RECO"),
                       cms.InputTag( SUPERCLUSTER_COLL_EE ,"", "RECO") )  
)

process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
   src = cms.InputTag("superClusters"),
   particleType = cms.int32(11),
)

#   Get the above SC's Candidates and place a cut on their Et and eta
process.goodSuperClusters = cms.EDFilter("CandViewSelector",
      src = cms.InputTag("superClusterCands"),
      cut = cms.string( SUPERCLUSTER_CUTS ),
      filter = cms.bool(True)
)                                         
                                         

#### remove real jets (with high hadronic energy fraction) from SC collection
##### this improves the purity of the probe sample without affecting efficiency

process.JetsToRemoveFromSuperCluster = cms.EDFilter("CaloJetSelector",   
    src = cms.InputTag("ak5CaloJets"),
    cut = cms.string('pt>5 && energyFractionHadronic > 0.15')
)
process.goodSuperClustersClean = cms.EDProducer("CandViewCleaner",
    srcObject = cms.InputTag("goodSuperClusters"),
    module_label = cms.string(''),
    srcObjectsToRemove = cms.VInputTag(cms.InputTag("JetsToRemoveFromSuperCluster")),
    deltaRMin = cms.double(0.1)
)

process.sc_sequence = cms.Sequence(
    process.superClusters +
    process.superClusterCands +
    process.goodSuperClusters +
    process.JetsToRemoveFromSuperCluster +
    process.goodSuperClustersClean 
    )



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
    cut = cms.string( "(abs(superCluster.eta)<2.5) && pt > 25.0" )
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

process.GsfMatchedSuperClusterCandsNoTrig = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
   deltaR =  cms.untracked.double(0.3)
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
   primaryVertexInputTag   = cms.InputTag("offlinePrimaryVertices"),
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

process.PassingEle17Mass50HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
#    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle17Mass50HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)
#meri bari

process.PassingEle17Mass50HLTMediumTAG = process.PassingEle17Mass50HLT.clone()
process.PassingEle17Mass50HLTMediumTAG.InputProducer = cms.InputTag("PassingWPMediumTAG" )

process.PassingEle17Mass50HLTMediumTAG20 = process.PassingEle17Mass50HLT.clone()
process.PassingEle17Mass50HLTMediumTAG.InputProducer = cms.InputTag("PassingWPMediumTAG20" )

process.PassingEle17Mass50HLTTight= process.PassingEle17Mass50HLT.clone()
process.PassingEle17Mass50HLTTight.InputProducer = cms.InputTag("PassingWPTight" )

process.PassingEle17Mass50HLTTightTAG20= process.PassingEle17Mass50HLT.clone()
process.PassingEle17Mass50HLTTightTAG20.InputProducer = cms.InputTag("PassingWPTightTAG20" )

process.PassingEle17Mass50HLTTightTAG= process.PassingEle17Mass50HLT.clone()
process.PassingEle17Mass50HLTTightTAG.InputProducer = cms.InputTag("PassingWPTightTAG" )

process.PassingEle17Mass50HLTLoose= process.PassingEle17Mass50HLT.clone()
process.PassingEle17Mass50HLTLoose.InputProducer = cms.InputTag("PassingWPLoose" )

process.PassingEle17Mass50HLTVeto= process.PassingEle17Mass50HLT.clone()
process.PassingEle17Mass50HLTVeto.InputProducer = cms.InputTag("PassingWPVeto" )
#ithe tak



process.PassingEle17HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(True),
    hltTags =hltTagsPassingEle17HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingEle8HLT = process.PassingEle17HLT.clone()
process.PassingEle8HLT.hltTags = hltTagsPassingEle8HLT

process.PassingNotEle17HLT =  process.PassingEle17HLT.clone()
process.PassingNotEle17HLT.antiSelect = cms.untracked.bool(True)

process.PassingEle8NotEle17HLT = process.PassingEle8HLT.clone()
process.PassingEle8NotEle17HLT.InputProducer = cms.InputTag("PassingNotEle17HLT" )


##    _____             ____        __ _       _ _   _             
##   |_   _|_ _  __ _  |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __  
##     | |/ _` |/ _` | | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \ 
##     | | (_| | (_| | | |_| |  __/  _| | | | | | |_| | (_) | | | |
##     |_|\__,_|\__, | |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|
##              |___/
## 


process.TagMedium = process.PassingEle17Mass50HLT.clone()
process.TagTight = process.PassingEle17Mass50HLTTight.clone()
process.TagLoose = process.PassingEle17Mass50HLTLoose.clone()
process.TagVeto = process.PassingEle17Mass50HLTVeto.clone()
process.TagMedium.InputProducer = cms.InputTag( "PassingWPMedium" )
process.TagTight.InputProducer = cms.InputTag( "PassingWPTight" )
process.TagLoose.InputProducer = cms.InputTag( "PassingWPLoose" )
process.TagVeto.InputProducer = cms.InputTag( "PassingWPVeto" )

process.ele_sequence = cms.Sequence(
    process.pfParticleSelectionSequence + 
    process.eleIsoSequence + 
    process.goodElectrons +
    process.goodElectronsTAG +
    process.goodElectronsTAG20 +
    process.PassingEle8Mass50HLTGsf +
    process.GsfMatchedSuperClusterCandsNoTrig +
#    process.kt6PFJets *
    process.PassingWPMedium +
    process.PassingWPMediumTAG +
    process.PassingWPMediumTAG20 +
    process.PassingWPTight +
    process.PassingWPTightTAG +
    process.PassingWPTightTAG20 +
    process.PassingWPLoose +
    process.PassingWPVeto +
    process.PassingEle17Mass50HLT +
    process.PassingEle17Mass50HLTMediumTAG +
    process.PassingEle17Mass50HLTMediumTAG20 +
    process.PassingEle17Mass50HLTTight +
    process.PassingEle17Mass50HLTTightTAG +
    process.PassingEle17Mass50HLTTightTAG20 +
    process.PassingEle17Mass50HLTLoose +
    process.PassingEle17Mass50HLTVeto +
    process.PassingEle17HLT +
    process.PassingEle8HLT +
    process.PassingNotEle17HLT +
    process.PassingEle8NotEle17HLT +
    process.TagMedium +
    process.TagTight +
    process.TagLoose +
    process.TagVeto 

    )


##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
#  Tag & probe selection ######
process.tagTightSC = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingEle17Mass50HLTTightTAG goodSuperClustersClean"), # charge conjugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("40 < mass < 10000"),
)
process.tagTight20SC = process.tagTightSC.clone()
process.tagTight20SC.decay = cms.string("PassingEle17Mass50HLTTightTAG20 goodSuperClustersClean")
process.tagMediumSC = process.tagTightSC.clone()
process.tagMediumSC.decay = cms.string("PassingEle17Mass50HLTMediumTAG goodSuperClustersClean")
process.tagMedium20SC = process.tagTightSC.clone()
process.tagMedium20SC.decay = cms.string("PassingEle17Mass50HLTMediumTAG20 goodSuperClustersClean")

process.tagGsf = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingEle17Mass50HLT@+ goodElectrons@-"), # charge conjugate states are implied
    checkCharge = cms.bool(True),                           
    cut   = cms.string("40 < mass < 10000"),
)

process.Zsignal = process.tagGsf.clone()
process.Zsignal.decay = cms.string("PassingEle17HLT@+ PassingEle8HLT@-")

process.tagTightGsf = process.tagGsf.clone()
process.tagTightGsf.decay = cms.string("PassingEle17Mass50HLTTightTAG@+ goodElectrons@-")
process.tagTight20Gsf = process.tagGsf.clone()
process.tagTight20Gsf.decay = cms.string("PassingEle17Mass50HLTTightTAG20@+ goodElectrons@-")
process.tagMediumGsf = process.tagGsf.clone()
process.tagMediumGsf.decay = cms.string("PassingEle17Mass50HLTMediumTAG@+ goodElectrons@-")
process.tagMedium20Gsf = process.tagGsf.clone()
process.tagMedium20Gsf.decay = cms.string("PassingEle17Mass50HLTMediumTAG20@+ goodElectrons@-")

process.tagWPMedium = process.tagGsf.clone()
process.tagWPMedium.decay = cms.string("PassingEle17Mass50HLTTightTAG@+ PassingWPMedium@-")
process.tagWPLoose = process.tagGsf.clone()
process.tagWPLoose.decay = cms.string("PassingEle17Mass50HLTTightTAG@+ PassingWPLoose@-")

process.allTagsAndProbes = cms.Sequence(
    process.tagTightSC +
    process.tagTight20SC +
    process.tagMediumSC +
    process.tagMedium20SC +
    process.Zsignal +
    process.tagGsf +
    process.tagTightGsf+
    process.tagMediumGsf+
    process.tagTight20Gsf+
    process.tagMedium20Gsf+
    process.tagWPMedium+
    process.tagWPLoose 
)

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        
process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("TagTight"),
    distMin = cms.double(0.2),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)
process.McMatchGsf = process.McMatchTag.clone()
process.McMatchGsf.src = cms.InputTag("goodElectrons")
process.McMatchWPMedium = process.McMatchTag.clone()
process.McMatchWPMedium.src = cms.InputTag("PassingWPMedium")
process.McMatchWPLoose = process.McMatchTag.clone()
process.McMatchWPLoose.src = cms.InputTag("PassingWPLoose")

process.McMatchSC = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("goodSuperClustersClean"),
    distMin = cms.double(0.2),
    matched = cms.InputTag("genParticles")
)
    
process.mc_sequence = cms.Sequence(
   process.McMatchTag +
   process.McMatchGsf + 
   process.McMatchWPMedium +
   process.McMatchWPLoose +
   process.McMatchSC
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
    probe_gsfEle_HoverE         = cms.string("hcalOverEcal"),    
    probe_gsfEle_sigmaIetaIeta = cms.string("sigmaIetaIeta"),
  ## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_abseta    = cms.string("abs(superCluster.eta)"),
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
          passingHLT = cms.InputTag("PassingEle17Mass50HLT"),     
          isWPTight = cms.InputTag("PassingWPTight"),
          passingHLTTight = cms.InputTag("PassingEle17Mass50HLTTight"),
          isWPLoose = cms.InputTag("PassingWPLoose"),
          passingHLTLoose = cms.InputTag("PassingEle17Mass50HLTLoose"),
          isWPVeto = cms.InputTag("PassingWPVeto"),
          passingHLTVeto = cms.InputTag("PassingEle17Mass50HLTVeto"),
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

##    ____   ____       __     ____      __ 
##   / ___| / ___|      \ \   / ___|___ / _|
##   \___ \| |      _____\ \ | |  _/ __| |_ 
##    ___) | |___  |_____/ / | |_| \__ \  _|
##   |____/ \____|      /_/   \____|___/_|  
##
## super cluster --> gsf electron
process.SuperClusterToGsfElectron = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,#CommonStuffForGsfElectronProbe,
    # choice of tag and probe pairs, and arbitration                 
    tagProbePairs = cms.InputTag("tagTightSC"),
    arbitration   = cms.string("None"),                      
    flags = cms.PSet(
        probe_passingGsf = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig"),        
    ),
    probeMatches  = cms.InputTag("McMatchSC"),
    allProbes     = cms.InputTag("goodSuperClustersClean")
)
process.SuperClusterToGsfElectron20 = process.SuperClusterToGsfElectron.clone()
process.SuperClusterToGsfElectron20.tagProbePairs = cms.InputTag("tagTight20SC") 

process.SuperClusterToGsfElectronM20 = process.SuperClusterToGsfElectron.clone()
process.SuperClusterToGsfElectronM20.tagProbePairs = cms.InputTag("tagMedium20SC")
 
process.SuperClusterToGsfElectronM25 = process.SuperClusterToGsfElectron.clone()
process.SuperClusterToGsfElectronM25.tagProbePairs = cms.InputTag("tagMediumSC")

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
    probeMatches  = cms.InputTag("McMatchGsf"),
#    allProbes     = cms.InputTag("PassingEle8Mass50HLTGsf")
    allProbes     = cms.InputTag("goodElectrons")
)
process.GsfElectronToIdMedium20 = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdMedium20.tagProbePairs = cms.InputTag("tagTight20Gsf")
process.GsfElectronToIdMediumM20 = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdMediumM20.tagProbePairs = cms.InputTag("tagMedium20Gsf")
process.GsfElectronToIdMediumM25 = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdMediumM25.tagProbePairs = cms.InputTag("tagMediumGsf")

process.GsfElectronToIdTight = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdTight.flags = cms.PSet(probe_isWPTight = cms.InputTag("PassingWPTight"))
process.GsfElectronToIdTight20 = process.GsfElectronToIdTight.clone()
process.GsfElectronToIdTight20.tagProbePairs = cms.InputTag("tagTight20Gsf")
process.GsfElectronToIdTightM20 = process.GsfElectronToIdTight.clone()
process.GsfElectronToIdTightM20.tagProbePairs = cms.InputTag("tagMedium20Gsf")
process.GsfElectronToIdTightM25 = process.GsfElectronToIdTight.clone()
process.GsfElectronToIdTightM25.tagProbePairs = cms.InputTag("tagMediumGsf")

process.GsfElectronToIdLoose = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdLoose.flags = cms.PSet(probe_isWPTight = cms.InputTag("PassingWPLoose"))
process.GsfElectronToIdLoose20 = process.GsfElectronToIdLoose.clone()
process.GsfElectronToIdLoose20.tagProbePairs = cms.InputTag("tagTight20Gsf")
process.GsfElectronToIdLooseM20 = process.GsfElectronToIdLoose.clone()
process.GsfElectronToIdLooseM20.tagProbePairs = cms.InputTag("tagMedium20Gsf")
process.GsfElectronToIdLooseM25 = process.GsfElectronToIdLoose.clone()
process.GsfElectronToIdLooseM25.tagProbePairs = cms.InputTag("tagMediumGsf")

process.GsfElectronToIdVeto = process.GsfElectronToIdMedium.clone()
process.GsfElectronToIdVeto.flags = cms.PSet(probe_isWPTight = cms.InputTag("PassingWPVeto"))
process.GsfElectronToIdVeto20 = process.GsfElectronToIdVeto.clone()
process.GsfElectronToIdVeto20.tagProbePairs = cms.InputTag("tagTight20Gsf")
process.GsfElectronToIdVetoM20 = process.GsfElectronToIdVeto.clone()
process.GsfElectronToIdVetoM20.tagProbePairs = cms.InputTag("tagMedium20Gsf")
process.GsfElectronToIdVetoM25 = process.GsfElectronToIdVeto.clone()
process.GsfElectronToIdVetoM25.tagProbePairs = cms.InputTag("tagMediumGsf")

process.Zsignalcounter = cms.EDAnalyzer("TagProbeFitTreeProducer",
#    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("Zsignal"),
    arbitration   = cms.string("OnePair"),
    flags = cms.PSet(
#        probe_isWP95 = cms.InputTag("WP95MatchedSuperClusterCandsClean"),
        probe_isWPMedium = cms.InputTag("PassingWPMedium"),
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
#    allProbes     = cms.InputTag("PassingEle8Mass50HLTGsf")
    allProbes     = cms.InputTag("PassingWPMedium")
)

##    ___    _       __    _   _ _   _____ 
##   |_ _|__| |      \ \  | | | | | |_   _|
##    | |/ _` |  _____\ \ | |_| | |   | |  
##    | | (_| | |_____/ / |  _  | |___| |  
##   |___\__,_|      /_/  |_| |_|_____|_|
##
##  offline selection --> HLT. First specify which quantities to store in the TP tree. 
##  WPMedium --> HLTEle17
if MC_flag:
    HLTmcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_flag),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_flag),
        checkMotherInUnbiasEff = cms.bool(MC_flag),
        mcVariables = cms.PSet(
          probe_eta = cms.string("eta"),
          probe_abseta = cms.string("abs(eta)"),
          probe_phi  = cms.string("phi"),
          probe_et  = cms.string("et"),
          probe_pt  = cms.string("pt"),
          probe_charge = cms.string("charge"),
        ),
        mcFlags     =  cms.PSet(
          probe_flag = cms.string("pt>0")
        ),
        )
else:
     HLTmcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )

process.WPMediumToHLTEle17 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    HLTmcTruthCommonStuff,
    variables = cms.PSet(
      probe_gsfEle_eta = cms.string("eta"),
      probe_gsfEle_abseta = cms.string("abs(eta)"),
      probe_gsfEle_et  = cms.string("et"),
      probe_gsfEle_pt  = cms.string("pt"),
      probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),
      probe_sc_eta    = cms.string("superCluster.eta"),
      probe_sc_abseta    = cms.string("abs(superCluster.eta)"),
    ),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool (True),
    tagProbePairs = cms.InputTag("tagWPMedium"),
    arbitration   = cms.string("None"),
    flags = cms.PSet(
        probe_passingHLT = cms.InputTag("PassingEle17HLT")
    ),
    probeMatches  = cms.InputTag("McMatchWPMedium"),
    allProbes     = cms.InputTag("PassingWPMedium")
)

process.WPLooseToHLTEle17 = process.WPMediumToHLTEle17.clone()
process.WPLooseToHLTEle17.tagProbePairs = cms.InputTag("tagWPLoose")
 
##  WPMedium --> HLTEle8NotEle17
process.WPMediumToHLTEle8NotEle17 = process.WPMediumToHLTEle17.clone()
process.WPMediumToHLTEle8NotEle17.flags =cms.PSet(
        probe_passingHLT = cms.InputTag("PassingEle8NotEle17HLT"),
        probe_passingHLTEle8 = cms.InputTag("PassingEle8HLT")
    )

process.WPLooseToHLTEle8NotEle17 = process.WPMediumToHLTEle8NotEle17.clone()
process.WPLooseToHLTEle8NotEle17.tagProbePairs = cms.InputTag("tagWPLoose")


process.tree_sequence = cms.Sequence(
    process.SuperClusterToGsfElectron +
    process.SuperClusterToGsfElectronM20 
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
        process.sc_sequence +
        process.ele_sequence + 
        process.allTagsAndProbes +
        process.tree_sequence
        )
    
process.TFileService = cms.Service(
#    "TFileService", fileName = cms.string( OUTPUT_FILE_NAME )
    "TFileService", fileName = cms.string( '$outputFileName') 
)
