import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V11::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#      $inputFileNames
'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0001/6C495815-BCD3-E111-92A0-0025B3E05CC0.root',
'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0001/B8D6916C-CAD3-E111-9DDE-001E67398791.root',
'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/9E37831D-2CD4-E111-9B38-0025B3E05BC4.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/FE087384-FAD2-E111-992C-003048670B64.root'
                                                  ),
eventsToProcess=cms.untracked.VEventRange('1:28065077','1:28065427','1:28033364','1:54509262')
#eventsToProcess=cms.untracked.VEventRange('1:28065077','1:28065427','1:28033364','1:54509262','1:28033351','1:28033341','1:28033347')
#eventsToProcess = cms.untracked.VEventRange(
#'193557:8215621',
#'193557:35284992')
   
)
#
# Event output
#

process.load("Configuration.EventContent.EventContent_cff")
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("myTesthisto.root")
)


#
# particle flow isolation
#

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

process.demo = cms.EDAnalyzer('ElectronTests',
    electronsInputTag       = cms.InputTag("gsfElectrons"),
    conversionsInputTag     = cms.InputTag("allConversions"),
    beamSpotInputTag        = cms.InputTag("offlineBeamSpot"),
    rhoIsoInputTag          = cms.InputTag("kt6PFJets", "rho"),
    primaryVertexInputTag   = cms.InputTag("offlinePrimaryVertices"),
    isoValInputTags         = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                            cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                            cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),

    printDebug              = cms.bool(True)
)

process.p = cms.Path(process.pfiso * process.demo)
