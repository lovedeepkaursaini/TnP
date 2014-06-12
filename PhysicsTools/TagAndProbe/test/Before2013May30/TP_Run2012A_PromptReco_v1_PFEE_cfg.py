##  *****************************************************************************
##  * Project: CMS detector at the CERN
##  *
##  * Package: PhysicsTools/TagAndProbe
##  *
##  *
##  * Authors:
##  *
##  *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
##  *
##  * Description:
##  *   - Produces tag & probe TTree for further analysis and computing efficiency
##  *
##  * History:
##  *   
##  * 
##  *****************************************************************************/


import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso

process = cms.Process("p")

MC_flag = False
GLOBAL_TAG = 'GR_R_44_V13::All'
if MC_flag:
    GLOBAL_TAG = 'START44_V9B::All'
    

HLTEle32Path1 = "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC170_v5" 
HLTEle32Path2 = "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4" 
HLTEle32Path3 = "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5" 

HLTEle17Path1 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"
HLTEle17Path2 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16"
HLTEle17Path3 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17"

HLTSC17Path1 = "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v5" 
HLTSC17Path2 = "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4" 
HLTSC17Path3 = "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5" 


HLTProcessName = "HLT"
if MC_flag:
    HLTElePath = "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"
    HLTProcessName = "REDIGI311X"

OUTPUT_FILE_NAME = "try.root"


ELECTRON_ET_CUT_MIN = 20.0
ELECTRON_COLL = "gsfElectrons"
ELECTRON_CUTS = "ecalDrivenSeed==1 && (abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
####

PHOTON_COLL = "photons"
PHOTON_CUTS = "hadronicOverEm<0.15 && (abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566) && ((isEB && sigmaIetaIeta<0.01) || (isEE && sigmaIetaIeta<0.03)) && (superCluster.energy*sin(superCluster.position.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
####
SUPERCLUSTER_COLL_EB = "correctedHybridSuperClusters"
SUPERCLUSTER_COLL_EE = "correctedMulti5x5SuperClustersWithPreshower"
SUPERCLUSTER_CUTS = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>" + str(ELECTRON_ET_CUT_MIN)


JET_COLL = "ak5PFJets"
JET_CUTS = "abs(eta)<2.6 && chargedHadronEnergyFraction>0 && electronEnergyFraction<0.1 && nConstituents>1 && neutralHadronEnergyFraction<0.99 && neutralEmEnergyFraction<0.99" 
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
process.GlobalTag.globaltag = GLOBAL_TAG
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
'file:/hdfs/store/data/Run2011A/DoubleElectron/AOD/03Oct2011-v1/0000/00108EDF-28EF-E011-8414-0026189438F9.root'
                                                  )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    
process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")


process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
process.JetsForRho = cms.Sequence( process.kt6PFJets )

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

#  Photons!!! ################ 
process.goodPhotons = cms.EDFilter(
    "PhotonSelector",
    src = cms.InputTag( PHOTON_COLL ),
    cut = cms.string(PHOTON_CUTS)
    )


process.PassingSC17HLTSC = cms.EDProducer("trgMatchedCandidateProducer",    
    InputProducer = cms.InputTag("goodSuperClustersClean" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),                                   
    hltTags = cms.VInputTag(
     cms.InputTag(HLTSC17Path1, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DoubleEtFilter", HLTProcessName),
     cms.InputTag(HLTSC17Path2, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     cms.InputTag(HLTSC17Path3, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)



process.sc_sequence = cms.Sequence(
    process.superClusters +
    process.superClusterCands +
    process.goodSuperClusters +
    process.JetsToRemoveFromSuperCluster +
    process.goodSuperClustersClean +
    process.PassingSC17HLTSC +
    process.goodPhotons 
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

process.PassingSC17HLTGsf = cms.EDProducer("trgMatchedPatElectronProducer",    
    InputProducer = cms.InputTag("goodElectrons" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = cms.VInputTag(
     cms.InputTag(HLTSC17Path1, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     cms.InputTag(HLTSC17Path2, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     cms.InputTag(HLTSC17Path3, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)


process.GsfMatchedSuperClusterCands = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("PassingSC17HLTSC"),
   ReferenceElectronCollection = cms.untracked.InputTag("PassingSC17HLTGsf"),
   deltaR =  cms.untracked.double(0.3)
)


process.GsfMatchedSuperClusterCandsNoTrig = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
   deltaR =  cms.untracked.double(0.3)
)

process.GsfMatchedPhotonCands = process.GsfMatchedSuperClusterCands.clone()
process.GsfMatchedPhotonCands.src = cms.InputTag("goodPhotons")

            

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
)

##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
##   
# Trigger  ##################

process.PassingEle32HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = cms.VInputTag(
     cms.InputTag(HLTEle32Path1, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter", HLTProcessName),
     cms.InputTag(HLTEle32Path2, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter", HLTProcessName),
     cms.InputTag(HLTEle32Path3, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingEle17HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(True),
    hltTags = cms.VInputTag(
    cms.InputTag(HLTEle17Path1,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter", HLTProcessName),
    cms.InputTag(HLTEle17Path2,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter", HLTProcessName),
    cms.InputTag(HLTEle17Path3,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingEle8HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(False),
    hltTags = cms.VInputTag(
    cms.InputTag(HLTEle17Path1,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ", HLTProcessName),
    cms.InputTag(HLTEle17Path2,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ", HLTProcessName),
    cms.InputTag(HLTEle17Path3,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingNotEle17HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(True),
    antiSelect = cms.untracked.bool(True),
    hltTags = cms.VInputTag(
    cms.InputTag(HLTEle17Path1,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter", HLTProcessName),
    cms.InputTag(HLTEle17Path2,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter", HLTProcessName),
    cms.InputTag(HLTEle17Path3,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingEle8NotEle17HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingNotEle17HLT" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(True),
    hltTags = cms.VInputTag(
    cms.InputTag(HLTEle17Path1,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ", HLTProcessName),
    cms.InputTag(HLTEle17Path2,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ", HLTProcessName),
    cms.InputTag(HLTEle17Path3,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingSC17HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = cms.VInputTag(
     cms.InputTag(HLTSC17Path1, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     cms.InputTag(HLTSC17Path2, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     cms.InputTag(HLTSC17Path3, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingSC17HLTWPMedium = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = cms.VInputTag(
     cms.InputTag(HLTSC17Path1, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     cms.InputTag(HLTSC17Path2, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     cms.InputTag(HLTSC17Path3, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", HLTProcessName),
     ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)



##    _____      _                        _  __     __             
##   | ____|_  _| |_ ___ _ __ _ __   __ _| | \ \   / /_ _ _ __ ___ 
##   |  _| \ \/ / __/ _ \ '__| '_ \ / _` | |  \ \ / / _` | '__/ __|
##   | |___ >  <| ||  __/ |  | | | | (_| | |   \ V / (_| | |  \__ \
##   |_____/_/\_\\__\___|_|  |_| |_|\__,_|_|    \_/ \__,_|_|  |___/
##   
## Here we show how to use a module to compute an external variable
## process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
## ak5PFResidual.useCondDB = False

process.superClusterDRToNearestJet = cms.EDProducer("DeltaRNearestJetComputer",
    probes = cms.InputTag("goodSuperClusters"),
       # ^^--- NOTA BENE: if probes are defined by ref, as in this case, 
       #       this must be the full collection, not the subset by refs.
    objects = cms.InputTag(JET_COLL),
    objectSelection = cms.string(JET_CUTS + " && pt > 20.0"),
)
process.JetMultiplicityInSCEvents = cms.EDProducer("CandMultiplicityCounter",
    probes = cms.InputTag("goodSuperClusters"),
    objects = cms.InputTag(JET_COLL),
    objectSelection = cms.string(JET_CUTS + " && pt > 20.0"),
)

process.PhotonDRToNearestJet = process.superClusterDRToNearestJet.clone()
process.PhotonDRToNearestJet.probes =cms.InputTag("goodPhotons")
process.JetMultiplicityInPhotonEvents = process.JetMultiplicityInSCEvents.clone()
process.JetMultiplicityInPhotonEvents.probes = cms.InputTag("goodPhotons")

process.GsfDRToNearestJet = process.superClusterDRToNearestJet.clone()
process.GsfDRToNearestJet.probes = cms.InputTag( ELECTRON_COLL )
process.JetMultiplicityInGsfEvents = process.JetMultiplicityInSCEvents.clone()
process.JetMultiplicityInGsfEvents.probes = cms.InputTag( ELECTRON_COLL )

process.ext_ToNearestJet_sequence = cms.Sequence(
    #process.ak5PFResidual + 
    process.superClusterDRToNearestJet +
    process.JetMultiplicityInSCEvents +
    process.PhotonDRToNearestJet +
    process.JetMultiplicityInPhotonEvents +    
    process.GsfDRToNearestJet +
    process.JetMultiplicityInGsfEvents
    )


##    _____             ____        __ _       _ _   _             
##   |_   _|_ _  __ _  |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __  
##     | |/ _` |/ _` | | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \ 
##     | | (_| | (_| | | |_| |  __/  _| | | | | | |_| | (_) | | | |
##     |_|\__,_|\__, | |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|
##              |___/
## 


process.Tag = process.PassingEle32HLT.clone()
process.Tag.InputProducer = cms.InputTag( "PassingWPMedium" )
process.TagMatchedSuperClusterCandsClean = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("Tag"),
   deltaR =  cms.untracked.double(0.3)
)




process.TagMatchedPhotonCands = process.TagMatchedSuperClusterCandsClean.clone()
process.TagMatchedPhotonCands.src     = cms.InputTag("goodPhotons")

process.ele_sequence = cms.Sequence(
    process.pfParticleSelectionSequence + 
    process.eleIsoSequence + 
    process.goodElectrons +
    process.PassingSC17HLTGsf +
    process.GsfMatchedSuperClusterCands +
    process.GsfMatchedSuperClusterCandsNoTrig +
    process.GsfMatchedPhotonCands +
    process.PassingWPMedium +
    process.PassingEle32HLT +
    process.PassingEle17HLT +
    process.PassingEle8HLT +
    process.PassingNotEle17HLT +
    process.PassingEle8NotEle17HLT +
    process.PassingSC17HLT +
    process.PassingSC17HLTWPMedium +
    process.Tag +
    process.TagMatchedSuperClusterCandsClean +
    process.TagMatchedPhotonCands 
    )


##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
#  Tag & probe selection ######
process.tagSC = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingEle32HLT goodSuperClustersClean"), # charge conjugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("40 < mass < 10000"),
)

process.tagPhoton = process.tagSC.clone()
process.tagPhoton.decay = cms.string("PassingEle32HLT goodPhotons")
process.GsfGsf = process.tagSC.clone()
process.GsfGsf.decay = cms.string("goodElectrons goodElectrons")
process.tagGsf = process.tagSC.clone()
process.tagGsf.decay = cms.string("PassingEle32HLT goodElectrons")
process.tagWPMedium = process.tagSC.clone()
process.tagWPMedium.decay = cms.string("PassingEle32HLT PassingWPMedium")
process.Zsignal = process.tagSC.clone()
process.Zsignal.decay = cms.string("PassingEle17HLT PassingEle8HLT")
process.elecMet = process.tagSC.clone()
process.elecMet.decay = cms.string("pfMet PassingWPMedium")
process.elecMet.cut = cms.string("mt > 0")

process.CSVarsTagGsf = cms.EDProducer("ColinsSoperVariablesComputer",
    parentBoson = cms.InputTag("tagGsf")
)
process.CSVarsGsfGsf = process.CSVarsTagGsf.clone()
process.CSVarsGsfGsf.parentBoson = cms.InputTag("GsfGsf")



process.allTagsAndProbes = cms.Sequence(
    process.Zsignal +
    process.tagSC +
    process.tagPhoton +
    process.tagGsf +
    process.GsfGsf +
    process.tagWPMedium +
    process.elecMet + 
    process.CSVarsTagGsf +
    process.CSVarsGsfGsf
)

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        
process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("Tag"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)
process.McMatchSC = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("goodSuperClustersClean"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)
process.McMatchPhoton = process.McMatchSC.clone()
process.McMatchPhoton.src = cms.InputTag("goodPhotons")
process.McMatchGsf = process.McMatchTag.clone()
process.McMatchGsf.src = cms.InputTag("goodElectrons")
process.McMatchWPMedium = process.McMatchTag.clone()
process.McMatchWPMedium.src = cms.InputTag("PassingWPMedium")
    
process.mc_sequence = cms.Sequence(
   process.McMatchTag +
   process.McMatchSC +
   process.McMatchPhoton +
   process.McMatchGsf + 
   process.McMatchWPMedium 
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
    phi  = cms.string("phi"),
    et  = cms.string("et"),
    e  = cms.string("energy"),
    p  = cms.string("p"),
    px  = cms.string("px"),
    py  = cms.string("py"),
    pz  = cms.string("pz"),
    theta  = cms.string("theta"),    
    vx     = cms.string("vx"),
    vy     = cms.string("vy"),
    vz     = cms.string("vz"),
    rapidity  = cms.string("rapidity"),
    mass  = cms.string("mass"),
    mt  = cms.string("mt"),    
)   

ProbeVariablesToStore = cms.PSet(
    probe_gsfEle_eta = cms.string("eta"),
    probe_gsfEle_abseta = cms.string("abs(eta)"),
    probe_gsfEle_pt  = cms.string("pt"),
    probe_gsfEle_phi  = cms.string("phi"),
    probe_gsfEle_et  = cms.string("et"),
    probe_gsfEle_e  = cms.string("energy"),
    probe_gsfEle_p  = cms.string("p"),
    probe_gsfEle_px  = cms.string("px"),
    probe_gsfEle_py  = cms.string("py"),
    probe_gsfEle_pz  = cms.string("pz"),
    probe_gsfEle_theta  = cms.string("theta"),    
    probe_gsfEle_charge = cms.string("charge"),
    probe_gsfEle_rapidity  = cms.string("rapidity"),
    probe_gsfEle_missingHits = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits"),
    probe_gsfEle_convDist = cms.string("convDist"),
    probe_gsfEle_convDcot = cms.string("convDcot"),
    probe_gsfEle_convRadius = cms.string("convRadius"),        
    probe_gsfEle_hasValidHitInFirstPixelBarrel = cms.string("gsfTrack.hitPattern.hasValidHitInFirstPixelBarrel"),
    ## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    probe_sc_x      = cms.string("superCluster.x"),
    probe_sc_y      = cms.string("superCluster.y"),
    probe_sc_z      = cms.string("superCluster.z"),
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_abseta    = cms.string("abs(superCluster.eta)"),
    probe_sc_theta  = cms.string("superClusterPosition.theta"),   
    probe_sc_phi    = cms.string("superCluster.phi"),
    probe_sc_size   = cms.string("superCluster.size"), # number of hits
    ## track quantities
    probe_track_p      = cms.string("gsfTrack.p"),
    probe_track_pt     = cms.string("gsfTrack.pt"),    
    probe_track_px     = cms.string("gsfTrack.px"),
    probe_track_py     = cms.string("gsfTrack.py"),
    probe_track_pz     = cms.string("gsfTrack.pz"),
    probe_track_eta    = cms.string("gsfTrack.eta"),
    probe_track_theta  = cms.string("gsfTrack.theta"),   
    probe_track_phi    = cms.string("gsfTrack.phi"),
    probe_track_vx     = cms.string("gsfTrack.vx"),
    probe_track_vy     = cms.string("gsfTrack.vy"),
    probe_track_vz     = cms.string("gsfTrack.vz"),    
    probe_track_dxy    = cms.string("gsfTrack.dxy"),
    probe_track_d0     = cms.string("gsfTrack.d0"),
    probe_track_dsz    = cms.string("gsfTrack.dsz"),
    probe_track_charge = cms.string("gsfTrack.charge"),
    probe_track_qoverp = cms.string("gsfTrack.qoverp"),
    probe_track_normalizedChi2 = cms.string("gsfTrack.normalizedChi2"),
    ## isolation 
    probe_gsfEle_trackiso = cms.string("dr03TkSumPt"),
    probe_gsfEle_ecaliso  = cms.string("dr03EcalRecHitSumEt"),
    probe_gsfEle_hcaliso  = cms.string("dr03HcalTowerSumEt"),
    ## classification, location, etc.    
    probe_gsfEle_classification = cms.string("classification"),
    probe_gsfEle_numberOfBrems  = cms.string("numberOfBrems"),     
    probe_gsfEle_bremFraction   = cms.string("fbrem"),
    probe_gsfEle_mva            = cms.string("mva"),        
    probe_gsfEle_deltaEta       = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    probe_gsfEle_deltaPhi       = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    probe_gsfEle_deltaPhiOut    = cms.string("deltaPhiSeedClusterTrackAtCalo"),
    probe_gsfEle_deltaEtaOut    = cms.string("deltaEtaSeedClusterTrackAtCalo"),
    probe_gsfEle_isEB           = cms.string("isEB"),
    probe_gsfEle_isEE           = cms.string("isEE"),
    probe_gsfEle_isGap          = cms.string("isGap"),
    ## Hcal energy over Ecal Energy
    probe_gsfEle_HoverE         = cms.string("hcalOverEcal"),    
    probe_gsfEle_EoverP         = cms.string("eSuperClusterOverP"),
    probe_gsfEle_eSeedClusterOverP = cms.string("eSeedClusterOverP"),    
    ## Cluster shape information
    probe_gsfEle_sigmaEtaEta  = cms.string("sigmaEtaEta"),
    probe_gsfEle_sigmaIetaIeta = cms.string("sigmaIetaIeta"),
    probe_gsfEle_e1x5               = cms.string("e1x5"),
    probe_gsfEle_e2x5Max            = cms.string("e2x5Max"),
    probe_gsfEle_e5x5               = cms.string("e5x5"),
    ## is ECAL driven ? is Track driven ?
    probe_gsfEle_ecalDrivenSeed     = cms.string("ecalDrivenSeed"),
    probe_gsfEle_trackerDrivenSeed  = cms.string("trackerDrivenSeed")
)


TagVariablesToStore = cms.PSet(
    gsfEle_eta = cms.string("eta"),
    gsfEle_abseta = cms.string("abs(eta)"),
    gsfEle_pt  = cms.string("pt"),
    gsfEle_phi  = cms.string("phi"),
    gsfEle_et  = cms.string("et"),
    gsfEle_e  = cms.string("energy"),
    gsfEle_p  = cms.string("p"),
    gsfEle_px  = cms.string("px"),
    gsfEle_py  = cms.string("py"),
    gsfEle_pz  = cms.string("pz"),
    gsfEle_theta  = cms.string("theta"),    
    gsfEle_charge = cms.string("charge"),
    gsfEle_rapidity  = cms.string("rapidity"),
    gsfEle_missingHits = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits"),
    gsfEle_convDist = cms.string("convDist"),
    gsfEle_convDcot = cms.string("convDcot"),
    gsfEle_convRadius = cms.string("convRadius"),     
    gsfEle_hasValidHitInFirstPixelBarrel = cms.string("gsfTrack.hitPattern.hasValidHitInFirstPixelBarrel"),
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    sc_x      = cms.string("superCluster.x"),
    sc_y      = cms.string("superCluster.y"),
    sc_z      = cms.string("superCluster.z"),
    sc_eta    = cms.string("superCluster.eta"),
    sc_abseta    = cms.string("abs(superCluster.eta)"),
    sc_theta  = cms.string("superClusterPosition.theta"),      
    sc_phi    = cms.string("superCluster.phi"),
    sc_size   = cms.string("superCluster.size"), # number of hits
    ## track quantities
    track_p      = cms.string("gsfTrack.p"),
    track_pt     = cms.string("gsfTrack.pt"),    
    track_px     = cms.string("gsfTrack.px"),
    track_py     = cms.string("gsfTrack.py"),
    track_pz     = cms.string("gsfTrack.pz"),
    track_eta    = cms.string("gsfTrack.eta"),
    track_theta  = cms.string("gsfTrack.theta"),   
    track_phi    = cms.string("gsfTrack.phi"),
    track_vx     = cms.string("gsfTrack.vx"),
    track_vy     = cms.string("gsfTrack.vy"),
    track_vz     = cms.string("gsfTrack.vz"),    
    track_dxy    = cms.string("gsfTrack.dxy"),
    track_d0     = cms.string("gsfTrack.d0"),
    track_dsz    = cms.string("gsfTrack.dsz"),
    track_charge = cms.string("gsfTrack.charge"),
    track_qoverp = cms.string("gsfTrack.qoverp"),
    track_normalizedChi2 = cms.string("gsfTrack.normalizedChi2"),    
    ## isolation 
    gsfEle_trackiso = cms.string("dr03TkSumPt"),
    gsfEle_ecaliso  = cms.string("dr03EcalRecHitSumEt"),
    gsfEle_hcaliso  = cms.string("dr03HcalTowerSumEt"),
    ## classification, location, etc.    
    gsfEle_classification = cms.string("classification"),
    gsfEle_numberOfBrems  = cms.string("numberOfBrems"),     
    gsfEle_bremFraction   = cms.string("fbrem"),
    gsfEle_mva            = cms.string("mva"),        
    gsfEle_deltaEta       = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    gsfEle_deltaPhi       = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    gsfEle_deltaPhiOut    = cms.string("deltaPhiSeedClusterTrackAtCalo"),
    gsfEle_deltaEtaOut    = cms.string("deltaEtaSeedClusterTrackAtCalo"),
    gsfEle_isEB           = cms.string("isEB"),
    gsfEle_isEE           = cms.string("isEE"),
    gsfEle_isGap          = cms.string("isGap"),
    ## Hcal energy over Ecal Energy
    gsfEle_HoverE         = cms.string("hcalOverEcal"),    
    gsfEle_EoverP         = cms.string("eSuperClusterOverP"),
    gsfEle_eSeedClusterOverP = cms.string("eSeedClusterOverP"),  
    ## Cluster shape information
    gsfEle_sigmaEtaEta  = cms.string("sigmaEtaEta"),
    gsfEle_sigmaIetaIeta = cms.string("sigmaIetaIeta"),
    gsfEle_e1x5               = cms.string("e1x5"),
    gsfEle_e2x5Max            = cms.string("e2x5Max"),
    gsfEle_e5x5               = cms.string("e5x5"),
    ## is ECAL driven ? is Track driven ?
    gsfEle_ecalDrivenSeed     = cms.string("ecalDrivenSeed"),
    gsfEle_trackerDrivenSeed  = cms.string("trackerDrivenSeed")
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
          passingHLT = cms.InputTag("PassingEle32HLT")     
    ),    
)

CommonStuffForSuperClusterProbe = CommonStuffForGsfElectronProbe.clone()
CommonStuffForSuperClusterProbe.variables = cms.PSet(
    probe_eta = cms.string("eta"),
    probe_abseta = cms.string("abs(eta)"),
    probe_pt  = cms.string("pt"),
    probe_phi  = cms.string("phi"),
    probe_et  = cms.string("et"),
    probe_e  = cms.string("energy"),
    probe_p  = cms.string("p"),
    probe_px  = cms.string("px"),
    probe_py  = cms.string("py"),
    probe_pz  = cms.string("pz"),
    probe_theta  = cms.string("theta"),
    )


if MC_flag:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_flag),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_flag),
        checkMotherInUnbiasEff = cms.bool(MC_flag),
        mcVariables = cms.PSet(
        probe_eta = cms.string("eta"),
         probe_abseta = cms.string("abs(eta)"),
       probe_pt  = cms.string("pt"),
        probe_phi  = cms.string("phi"),
        probe_et  = cms.string("et"),
        probe_e  = cms.string("energy"),
        probe_p  = cms.string("p"),
        probe_px  = cms.string("px"),
        probe_py  = cms.string("py"),
        probe_pz  = cms.string("pz"),
        probe_theta  = cms.string("theta"),    
        probe_vx     = cms.string("vx"),
        probe_vy     = cms.string("vy"),
        probe_vz     = cms.string("vz"),   
        probe_charge = cms.string("charge"),
        probe_rapidity  = cms.string("rapidity"),    
        probe_mass  = cms.string("mass"),
        probe_mt  = cms.string("mt"),    
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
    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration                 
    tagProbePairs = cms.InputTag("tagSC"),
    arbitration   = cms.string("Random2"),                      
    flags = cms.PSet(
        probe_passingGsf = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig"),        
        probe_passingHLT = cms.InputTag("TagMatchedSuperClusterCandsClean")
    ),
    probeMatches  = cms.InputTag("McMatchSC"),
    allProbes     = cms.InputTag("goodSuperClustersClean")
)
process.SuperClusterToGsfElectron.variables.probe_dRjet = cms.InputTag("superClusterDRToNearestJet")
process.SuperClusterToGsfElectron.variables.probe_nJets = cms.InputTag("JetMultiplicityInSCEvents")
#process.SuperClusterToGsfElectron.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
#process.SuperClusterToGsfElectron.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")


process.NonIsoTagSuperClusterToGsfElectron = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration                 
    tagProbePairs = cms.InputTag("nonisotagSC"),
    arbitration   = cms.string("Random2"),                      
    flags = cms.PSet(
        probe_passingGsf = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig"),        
        probe_passingHLT = cms.InputTag("TagMatchedSuperClusterCandsClean")
    ),
    probeMatches  = cms.InputTag("McMatchSC"),
    allProbes     = cms.InputTag("goodSuperClustersClean")
)
process.SuperClusterToGsfElectron.variables.probe_dRjet = cms.InputTag("superClusterDRToNearestJet")
process.SuperClusterToGsfElectron.variables.probe_nJets = cms.InputTag("JetMultiplicityInSCEvents")
#process.SuperClusterToGsfElectron.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
#process.SuperClusterToGsfElectron.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")



## good photon --> gsf electron
process.PhotonToGsfElectron = process.SuperClusterToGsfElectron.clone()
process.PhotonToGsfElectron.tagProbePairs = cms.InputTag("tagPhoton")
process.PhotonToGsfElectron.flags = cms.PSet(
    probe_passingGsf = cms.InputTag("GsfMatchedPhotonCands"),
    probe_passingHLT = cms.InputTag("TagMatchedPhotonCands"),
    )
process.PhotonToGsfElectron.probeMatches  = cms.InputTag("McMatchPhoton")
process.PhotonToGsfElectron.allProbes     = cms.InputTag("goodPhotons")
process.PhotonToGsfElectron.variables.probe_dRjet = cms.InputTag("PhotonDRToNearestJet")
process.PhotonToGsfElectron.variables.probe_nJets = cms.InputTag("JetMultiplicityInPhotonEvents")
process.PhotonToGsfElectron.variables.probe_trackiso = cms.string("trkSumPtHollowConeDR03")
process.PhotonToGsfElectron.variables.probe_ecaliso = cms.string("ecalRecHitSumEtConeDR03")
process.PhotonToGsfElectron.variables.probe_hcaliso = cms.string("hcalTowerSumEtConeDR03")
process.PhotonToGsfElectron.variables.probe_HoverE  = cms.string("hadronicOverEm")
process.PhotonToGsfElectron.variables.probe_sigmaIetaIeta = cms.string("sigmaIetaIeta")
process.PhotonToGsfElectron.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
process.PhotonToGsfElectron.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")


##   ____      __       __    ___                 ___    _ 
##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
##                                           |/            
##  gsf electron --> isolation, electron id  etc.
process.GsfElectronToId = cms.EDAnalyzer("TagProbeFitTreeProducer",
#    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("tagGsf"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
#        probe_isWP95 = cms.InputTag("WP95MatchedSuperClusterCandsClean"),
        probe_isWPMedium = cms.InputTag("PassingWPMedium"),
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
#    allProbes     = cms.InputTag("PassingSC17HLTGsf")
    allProbes     = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig")
)

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
#    allProbes     = cms.InputTag("PassingSC17HLTGsf")
    allProbes     = cms.InputTag("PassingWPMedium")
)


#process.GsfElectronToId.variables.probe_dRjet = cms.InputTag("GsfDRToNearestJet")
#process.GsfElectronToId.variables.probe_nJets = cms.InputTag("JetMultiplicityInGsfEvents")
#process.GsfElectronToId.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
#process.GsfElectronToId.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")
#process.GsfElectronToId.pairVariables.costheta = cms.InputTag("CSVarsTagGsf","costheta")
#process.GsfElectronToId.pairVariables.sin2theta = cms.InputTag("CSVarsTagGsf","sin2theta")
#process.GsfElectronToId.pairVariables.tanphi = cms.InputTag("CSVarsTagGsf","tanphi")


process.GsfElectronPlusGsfElectron = process.GsfElectronToId.clone()
process.GsfElectronPlusGsfElectron.tagProbePairs = cms.InputTag("GsfGsf")
process.GsfElectronPlusGsfElectron.tagMatches = cms.InputTag("McMatchGsf")
process.GsfElectronPlusGsfElectron.pairVariables.costheta = cms.InputTag("CSVarsGsfGsf","costheta")
process.GsfElectronPlusGsfElectron.pairVariables.sin2theta = cms.InputTag("CSVarsGsfGsf","sin2theta")
process.GsfElectronPlusGsfElectron.pairVariables.tanphi = cms.InputTag("CSVarsGsfGsf","tanphi")


process.GsfElectronPlusMet = process.GsfElectronToId.clone()
process.GsfElectronPlusMet.tagProbePairs = cms.InputTag("elecMet")
process.GsfElectronPlusMet.tagVariables = cms.PSet()
process.GsfElectronPlusMet.pairVariables =  cms.PSet(ZVariablesToStore)
process.GsfElectronPlusMet.pairFlags =  cms.PSet( isMTabove40 = cms.string("mt > 40") )
process.GsfElectronPlusMet.isMC = cms.bool(False)


##    ___    _       __    _   _ _   _____ 
##   |_ _|__| |      \ \  | | | | | |_   _|
##    | |/ _` |  _____\ \ | |_| | |   | |  
##    | | (_| | |_____/ / |  _  | |___| |  
##   |___\__,_|      /_/  |_| |_|_____|_|
##
##  offline selection --> HLT. First specify which quantities to store in the TP tree. 
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

##  WPMedium --> HLTEle17
process.WPMediumToHLTEle17 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    HLTmcTruthCommonStuff,                                
    variables = cms.PSet(
      probe_gsfEle_eta = cms.string("eta"),
      probe_gsfEle_abseta = cms.string("abs(eta)"),
      probe_gsfEle_phi  = cms.string("phi"),
      probe_gsfEle_et  = cms.string("et"),
      probe_gsfEle_charge = cms.string("charge"),
      probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
      probe_sc_eta    = cms.string("superCluster.eta"), 
      probe_sc_abseta    = cms.string("abs(superCluster.eta)"), 
      probe_sc_phi    = cms.string("superCluster.phi"),
      probe_gsfEle_isEB           = cms.string("isEB"),
      probe_gsfEle_isEE           = cms.string("isEE"),
      probe_gsfEle_isGap          = cms.string("isGap"),
    ),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool (True),                                                        
    tagProbePairs = cms.InputTag("tagWPMedium"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet( 
        probe_passingHLT = cms.InputTag("PassingEle17HLT")        
    ),
    probeMatches  = cms.InputTag("McMatchWPMedium"),
    allProbes     = cms.InputTag("PassingWPMedium")
)

##  WPMedium --> HLTEle8NotEle17
process.WPMediumToHLTEle8NotEle17 = process.WPMediumToHLTEle17.clone()
process.WPMediumToHLTEle8NotEle17.tagProbePairs = cms.InputTag("tagWPMedium")
process.WPMediumToHLTEle8NotEle17.probeMatches  = cms.InputTag("McMatchWPMedium")
process.WPMediumToHLTEle8NotEle17.allProbes     = cms.InputTag("PassingWPMedium")
process.WPMediumToHLTEle8NotEle17.flags =cms.PSet( 
        probe_passingHLT = cms.InputTag("PassingEle8NotEle17HLT")        
    ) 

process.tree_sequence = cms.Sequence(
    process.SuperClusterToGsfElectron +
    process.GsfElectronToId +
    process.Zsignalcounter +
    process.WPMediumToHLTEle17 +
    process.WPMediumToHLTEle8NotEle17 
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

if MC_flag:
    process.p = cms.Path(
        process.JetsForRho + process.sc_sequence + process.eIDSequence + process.ele_sequence + 
        process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.mc_sequence + 
        process.tree_sequence
        )
else:
    process.p = cms.Path(
        process.sc_sequence +
        process.ele_sequence + 
        process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.tree_sequence
        )
    
process.TFileService = cms.Service(
    "TFileService", fileName = cms.string( OUTPUT_FILE_NAME )
    )
