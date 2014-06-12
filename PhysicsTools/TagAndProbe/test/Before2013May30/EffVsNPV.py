import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

isMC = True
InputFileName = "/hdfs/store/user/anil79/mrgdTnPTreeElec_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraphSummer12_DR53X-PU_S10_START53_V7A-v1_Jan13.root"
OutputFilePrefix = "/scratch/kaur/MCMedIDEffiVsNVtx"

################################################
PDFName = "pdfSignalPlusBackground"
################################################

#specifies the binning of parameters
EfficiencyBins = cms.PSet(
event_nPV = cms.vdouble( 0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60),
    probe_sc_abseta = cms.vdouble( 0,1.4442,1.556,2.5 ),
#    run = cms.vdouble(160431, 167676, 173198, 176929, 178151, 180252),
 #   run = cms.vdouble(190000, 200000),
)

#### For data: except for HLT step
EfficiencyBinningSpecification = cms.PSet(
    #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
    UnbinnedVariables = cms.vstring("mass"),
    #specifies the binning of parameters
    BinnedVariables = cms.PSet(EfficiencyBins),
    #first string is the default followed by binRegExp - PDFname pairs
    BinToPDFmap = cms.vstring(PDFName)
)

#### For MC truth: do truth matching
EfficiencyBinningSpecificationMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(
event_nPV = cms.vdouble( 0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60),
    probe_sc_abseta = cms.vdouble( 0,0.8,1.4442,1.556,2.0,2.5 ),
#    run = cms.vdouble(160431, 167676, 173198, 176929, 178151, 180252),
#    run = cms.vdouble(190000, 200000),
    mcTrue = cms.vstring("true")
    ),
    BinToPDFmap = cms.vstring()  
)

##########################################################################################
############################################################################################
if isMC:
    mcTruthModules = cms.PSet(
        MCtruth_WPMedium = cms.PSet(
        EfficiencyBinningSpecificationMC,   
        EfficiencyCategoryAndState = cms.vstring("probe_isWPMedium","pass"),
        ),
        )
else:
    mcTruthModules = cms.PSet()
##########################################################################################

############################################################################################
####### GsfElectron->Id / selection efficiency 
############################################################################################

process.GsfElectronToIdMedium = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(InputFileName),
    InputDirectoryName = cms.string("GsfElectronToIdMedium"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(OutputFilePrefix+"GsfElectronToIdMedium.root"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    #fixVars = cms.vstring("mean"),
                                                 
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        event_nPV = cms.vstring("event_nPV", "0", "70", ""),
        probe_sc_abseta = cms.vstring("Probe |#eta|", "0.", "2.5", ""),
#                    run = cms.vstring("Run number", "160431", "180252", ""),
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        probe_isWPMedium = cms.vstring("probe_isWPMedium", "dummy[pass=1,fail=0]"),
        ),
    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        pdfSignalPlusBackground = cms.vstring(
"BreitWigner::physicsShape(mass,m0[91.1876,80,100],fwhm[2.495])",  #replace with floating + constriant later?
            "CBShape::resolution(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
            "FCONV::signal(mass,physicsShape,resolution)",
#            "Exponential::backgroundPass(mass, lp[0,-5,5])",
#            "Exponential::backgroundFail(mass, lf[0,-5,5])",
    "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
    "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]" 

##     "CBExGaussShape::signalRes(mass, mean[2.0946e-01], sigma[8.5695e-04],alpha[3.8296e-04], n[6.7489e+00], sigma_2[2.5849e+00], frac[6.5704e-01])",  
### the signal function goes here
#     "CBExGaussShape::signalResPass(mass, meanP[0.], sigmaP[8.5695e-04, 0., 3.],alphaP[3.8296e-04], nP[6.7489e+00], sigmaP_2[2.5849e+00], fracP[6.5704e-01])",  ### signal resolution for "pass" sample
 #    "CBExGaussShape::signalResFail(mass, meanF[2.0946e-01, -5., 5.], sigmaF[8.5695e-04, 0., 5.],alphaF[3.8296e-04], nF[6.7489e+00], sigmaF_2[2.5849e+00], fracF[6.5704e-01])",  ### signal resolution for "fail" sample     
  #  "ZGeneratorLineShape::signalPhy(mass)", ### NLO line shape
   # "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
    #"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
#    "FCONV::signalPass(mass, signalPhy, signalResPass)",
 #   "FCONV::signalFail(mass, signalPhy, signalResFail)",     
  #  "efficiency[0.9,0,1]",
#    "signalFractionInPassing[1.0]"     
   # "Gaussian::signal(mass, mean[91.2, 89.0, 93.0], sigma[2.3, 0.5, 10.0])",
#    "RooExponential::backgroundPass(mass, cPass[-0.02,-5,0])",
 #   "RooExponential::backgroundFail(mass, cFail[-0.02,-5,0])",
  #  "efficiency[0.9,0,1]",
#    "signalFractionInPassing[0.9]"
            #"BreitWigner::signal(mass,m0[91.1876,80,100],fwhm[2.495])",  #replace with floating + constriant later?
#"BreitWigner::physicsShape(mass,m0[91.1876,80,100],fwhm[2.495])",  #replace with floating + constriant later?
 #           "CBShape::resolution(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
  #          "FCONV::signal(mass,physicsShape,resolution)",
#            "Exponential::backgroundPass(mass, lp[0,-5,5])",
 #           "Exponential::backgroundFail(mass, lf[0,-5,5])",
  #          "efficiency[0.9,0,1]",
#            "signalFractionInPassing[0.9]" 

       ),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
    mcTruthModules,
##     #the name of the parameter set becomes the name of the directory
    #WPMedium = cms.PSet(
   # EfficiencyBinningSpecification,   
   # EfficiencyCategoryAndState = cms.vstring("probe_isWPMedium","pass"),
    ),
    )


process.GsfElectronToIdMedium.Variables.WeightVariable = cms.string("PUweight")
#process.GsfElectronToIdMedium.Efficiencies.WPMedium.BinToPDFmap= cms.vstring()


process.fit = cms.Path(
    process.GsfElectronToIdMedium   
 #  process.SCToGsfElectron  
#     process.WPMediumToHLTEle +
    # process.WPMediumToHLTEle17  
#     process.WPMediumToHLTEle8NotEle17 
    )
