
import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

fileName = "2011mrgdTnPTree_DYToEEFall11POWHEG"
stage = "GsfElectronToIdLoose"
pred = "probe_isWPLoose"

import os
#os.system("rm ZeeGenLevel.root")
#os.system("ln -s Templates/"+stage+".root ./ZeeGenLevel.root")

from PhysicsTools.TagAndProbe.FitConf.tnpEffiStudies_cfi import *
from PhysicsTools.TagAndProbe.FitConf.tnpFitModels_cfi import *

truth = 1
treeName = "treeForEffi"
if (truth):
	treeName = "fitter_tree"
else:
	treeName = "mcUnbias_tree"


##one more function is handy
def SetTreeAnalyzer( anlzr, fileName, pred, ptBin, etaBin,ptBinN,etaBinN):
	binVar = ptEtaBins.clone()
	binVar.BinnedVariables.probe_sc_et = ptBin#cms.vdouble(20,30)
	binVar.BinnedVariables.probe_sc_abseta =etaBin#cms.vdouble(0,1.44)
	sigModel = surrogate.clone()
	SetSignalShape(sigModel,str(ptBin[0]),str(ptBin[1]),str(etaBin[0]),str(etaBin[1]))
	ConfigureTPAnalyzer(anlzr, fileName, pred, binVar, treeName)
	ActivateAndConfigureFit(anlzr, sigModel, True, 40)
	print anlzr.PDFs.fittingFunction
								    
###################################################################################
# All the processes below should now be generated according
# to new scheme. And also they should be generated using a  logic loop...
####################################################################################

process.TnP_20To30_0To0p8 = CommonDummyTree.clone()
SetTreeAnalyzer(process.TnP_20To30_0To0p8,fileName,pred,ptBin=cms.vdouble(20,30),etaBin=cms.vdouble(0,0.8),ptBinN='"20To30"',etaBinN='"0To0.8"')



##define the efficiency step
##create the subprocessces
analyzers = ["TnP_20To30_0To0p8"]
for anlzr in analyzers:
	subProcess = stage+"_"+anlzr
	setattr(process,subProcess,getattr(process,anlzr))
	getattr(process,subProcess).InputDirectoryName=stage
	getattr(process,subProcess).OutputFileName="/scratch/anil79/"+subProcess+"_"+fileName+"_PtEtaBins_"+pred+"_"+treeName+".root"
##create a execution sequence
seq = getattr(process, stage+"_"+analyzers[0])
for i in range(1, len(analyzers)):
	subProcess = getattr(process,stage+"_"+analyzers[i])
	seq = seq + subProcess
	print "I am going to run following processes: \n\t", seq

process.myseq = cms.Sequence(seq)

process.fit = cms.Path(process.myseq)
