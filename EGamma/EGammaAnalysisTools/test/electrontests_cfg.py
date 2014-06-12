import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V11::All'
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#      $inputFileNames
'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0001/6C495815-BCD3-E111-92A0-0025B3E05CC0.root',
#%'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0001/B8D6916C-CAD3-E111-9DDE-001E67398791.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/9E37831D-2CD4-E111-9B38-0025B3E05BC4.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/FE087384-FAD2-E111-992C-003048670B64.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0001/94F2ED39-09D3-E111-B06D-00304867400E.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0001/627CBC27-FED3-E111-940B-001E673989DF.root',
#why i fail 
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/AE6FBBD1-3FD2-E111-BA56-001E67396E32.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/1C4AC243-F7D1-E111-8DBB-001E67397D73.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/7644D67E-6AD4-E111-957E-003048D476E2.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/DC117077-27D4-E111-A15F-002481E941F4.root',
#'file:/hdfs/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/F641985A-EDD1-E111-8E80-001E673986A6.root'
                                                  ),
#eventsToProcess=cms.untracked.VEventRange(
#'1:1:1870391',
#'1:2128449'
#,'1:3518872','1:1251813','1:1770393',

#'1:28065077','1:28065427',
#'1:67178490',
#'1:54509262',
#'1:28033364',
#'1:61112890',
#'1:64966492',
#'1:839961',
#'1:68422360',
#'1:19697187',
#'1:23019670',
#'1:7943104',
#'1:45466058',
#'1:70174878',
#'1:8927048',
#'1:36563707'
#)
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
#    "TFileService", fileName = cms.string( '$outputFileName' )
    "TFileService",    fileName = cms.string("myTesthisto.root")
)


#
# particle flow isolation
#

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

process.demo = cms.EDAnalyzer('EGammaCutBasedEleIdAnalyzer',
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

process.p = cms.Path(process.fastFilter * process.pfiso * process.demo)
