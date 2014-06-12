
import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

ptbins = (5,10)
#ptbins = (0,5,10,15,20,30,50)
#ptbins = (0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)
etabins = (1.4442,1.566)
#etabins = (0.0,0.8,1.4442,1.566,2.0,2.5)

fileName = "mrgdTnPTreeElec_DoubleElectronRun2012_June10"
stage = "GsfElectronToIdMedium"
pred = "probe_isWPMedium"

import os
#os.system("rm ZeeGenLevel.root")
#os.system("ln -s Templates/"+stage+".root ./ZeeGenLevel.root")

from PhysicsTools.TagAndProbe.FitConf.noPUrw_nPVMediumtnpEffiStudies_cfi import *
from PhysicsTools.TagAndProbe.FitConf.tnpFitModels_cfi import *

truth = 1
treeName = "treeForEffi"
if (truth):
	treeName = "fitter_tree"
else:
	treeName = "mcUnbias_tree"


##one more function is handy
def SetTreeAnalyzer( anlzr, fileName, pred, ptBin, etaBin):
	binVar = ptEtaBins.clone()
	binVar.BinnedVariables.event_nPV = ptBin#cms.vdouble(20,30)
	binVar.BinnedVariables.probe_sc_abseta =etaBin#cms.vdouble(0,1.44)
#        sigModel =bwConvCB.clone()# bWxCBandChebBackground.clone() 
#        sigModel =bwCBandChebBackground.clone()
#        sigModel =bwCBandCMSShapeBackground.clone()
        sigModel = surrogate.clone()
        SetSignalShape(sigModel,str(ptBin[0]),str(ptBin[1]),str(etaBin[0]),str(etaBin[1]))
#	sigModel = SYSsurrogate.clone()
	ConfigureTPAnalyzer(anlzr, fileName, pred, binVar, treeName)
	ActivateAndConfigureFit(anlzr, sigModel, True, 30)
	#print anlzr.PDFs.fittingFunction

# Below we create subprocesses, to perform ZMass fit: one job per pt-Eta 2D-Bin.
# The binning specifications will be stored in form of tuples. Process names will be stored in form of lists
# TUPLE: A tuple is an immutable, iterable sequence of comma separated objects called items, may or may not having same type. (x,y,z)
# LIST: A list is a mutable, sequence of comma separated objects of same or different type [x,y,z] 

analyzers = [];
#ptbins = (10,15,20)
#ptbins = (20,30,40)
#ptbins = (40,50,200)
#ptbins = (9,12,15,18)
#ptbins = (0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60)
#etabins = (0.0,1.4442,1.556,2.5)
#ptbins = (10,15,20,30,40,50,200)
#etabins = (0.0,0.8,1.4442,1.556,2.0,2.5)
for i in range(0,len(ptbins)-1):
	 pl = ptbins[i]
	 ph = ptbins[i+1]
	 for j in range(0,len(etabins)-1):
		 el = etabins[j]
		 eh = etabins[j+1]
		 nam = "TnP_"+str(pl)+"To"+str(ph)+"_"+str(el)+"To"+str(eh)
		 if ('.' in nam):
			 nam = nam.replace('.',"p")#at . does not work in subprocess name.
		 analyzers.append(nam)
		 setattr(process,nam,CommonDummyTree.clone())
		 SetTreeAnalyzer(getattr(process,nam),fileName,pred,ptBin=cms.vdouble(pl,ph),etaBin=cms.vdouble(el,eh))
  
#set input directory name, and outfile (anil: consider moving this to settreeanalzr funct)
for anlzr in analyzers:
	subProcess = stage+"_"+anlzr
	setattr(process,subProcess,getattr(process,anlzr))
	getattr(process,subProcess).InputDirectoryName=stage
	getattr(process,subProcess).OutputFileName="/scratch/kaur/nPVData/"+subProcess+"_"+fileName+"_PtEtaBins_"+pred+"_"+treeName+".root"


##create a execution sequence
seq = getattr(process, stage+"_"+analyzers[0])
for i in range(1, len(analyzers)):
	subProcess = getattr(process,stage+"_"+analyzers[i])
	seq = seq + subProcess
print "I am going to run following processes: \n\t", seq

process.myseq = cms.Sequence(seq)

process.fit = cms.Path(process.myseq)
