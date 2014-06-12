import FWCore.ParameterSet.Config as cms

process = cms.Process("p")

##  
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#      $inputFileNames
'dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/data/Run2012C/DoubleElectron/AOD/22Jan2013-v1/20002/78704EA6-AD68-E211-9272-0030486792AC.root'
#'file:/hdfs/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/FEF61011-7E90-E211-A84B-003048FFD75C.root'
#'/hdfs/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/20000/FEAF0DFD-9667-E211-9A4E-0025905938A8.root'
                                                  ),
#eventsToProcess = cms.untracked.VEventRange(
#'193557:8215621',
#'193557:35284992')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    


process.out = cms.OutputModule("PoolOutputModule", 
       fileName = cms.untracked.string("/scratch/kaur/Run2012CDoubleElectron22Jan2013-v12000278704EA6-AD68-E211-9272-0030486792AC.root"),
       )
    
process.outpath = cms.EndPath(process.out)
#process.outpath.remove(process.out)
