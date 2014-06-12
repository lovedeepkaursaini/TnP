import FWCore.ParameterSet.Config as cms
surrogate = cms.PSet(
    fittingFunction   = cms.vstring(
#    "Gaussian::signalRes(mass, mean[91.2, 89.0, 93.0], sigma[2.3, 0.5, 10.0])",
    "RooCBShape::signalRes(mass, mean[0.,-5.0,1.0], sigma[2.0, 0.00001, 8.0], alpha[1.0,0.0002,5.0], n[6.7489,1.0,15.0])",  ### signal resolution for "pass" sample
#   "RooCBShape::signalRes(mass, mean[0.,-5.0,1.0], sigma[2.0, 0.00001, 8.0], alpha[1.0,-10.2,5.0], n[6.7489,-10.0,15.0])",  ### signal resolution for "pass" sample
#   "CBShape::signalRes(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
#   "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
    "ZGeneratorLineShape::signalPhy1(mass,0,30,2,2.5,1)",
    "ZGeneratorLineShape::signalPhy2(mass,20,30,2,2.5,0)",
#"RooPolynomial::backgroundPass(mass,{a0[0,-1000.,1000.],a1[0,-10000.,10000.],a2[0,-10000000.,10000000.]})",
#"RooPolynomial::backgroundFail(mass,{a0[0,-1000.,1000.],a1[0,-10000.,10000.],a2[0,-10000000.,10000000.]})",
#"RooPolynomial::backgroundPass(mass,{a0[0,-100.,100.],a1[0,-100000.,100000.]}) ",
#"RooPolynomial::backgroundFail(mass,{a0[0,-100.,100.],a1[0,-100000.,100000.]}) ",
#  "RooPolynomial::backgroundPass(mass,{a0[0,-10000.,10000.],a1[0,-10000.,10000.],a2[0,-1000000.,1000000.],a3[0,-1000000.,1000000.]}) ",
#  "RooPolynomial::backgroundFail(mass,{a0[0,-10000.,10000.],a1[0,-10000.,10000.],a2[0,-1000000.,1000000.],a3[0,-1000000.,1000000.]}) ",
#  "RooPolynomial::backgroundFail(mass,{a0f[0,-10000.,10000.],a1f[0,-10000.,10000.],a2f[0,-1000000.,1000000.],a3f[0,-1000000.,1000000.]}) ",
   "Exponential::backgroundPass(mass, lp[0,-1,1])",#10,15,20,40-500
   "Exponential::backgroundFail(mass, lf[0,-1,1])",#10,15,20,40-500
#   "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, -1.,0.5], betaPass, peakPass[90.0])",#for 20,30,40
#   "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, -1.,0.5], betaFail, peakFail[90.0])",#for 20,30,40
# "RooChebychev::backgroundPass(mass, {cPass[0,-2.5,2.5], cPass2[0,-12.5,12.5]})", #for 40-500
# "RooChebychev::backgroundFail(mass, {cFail[0,-2.5,2.5], cFail2[0,-12.5,12.5]})", #for 40-500
    "FCONV::signalPass(mass, signalPhy1, signalRes)",
    "FCONV::signalFail(mass, signalPhy2, signalRes)",
    "efficiency[0.7,0,1]",
    "signalFractionInPassing[0.7]"
    ),
    )
 

bwConvCB = cms.PSet(
    fittingFunction = cms.vstring(
        "BreitWigner::physicsShape(mass,m0[91.1876,80,100],fwhm[2.495])",
        "CBShape::resolution(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "FCONV::signal(mass,physicsShape,resolution)",
        "Exponential::backgroundPass(mass, lp[0,-5,5])",
        "Exponential::backgroundFail(mass, lf[0,-5,5])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
        ),
    )

bwCBandChebBackground =cms.PSet(
    fittingFunction= cms.vstring(
#      "RooCBShape::signalRes(mass, mean[0.,-5.0,1.0], sigma[2.0, 0.1, 8.0], alpha[1.0,0.2,5.0], n[6.7489,1.0,15.0])",  ### signal resolution for "pass" sample
      "RooCBShape::signalRes(mass, mean[0.,-5.0,1.0], sigma[2.0, 0.0001, 8.0], alpha[1.0,0.0002,5.0], n[6.7489,1.0,15.0])",  ### signal resolution for "pass" sample
      "BreitWigner::signalPhy(mass,ZMass[91.1876,80,100], ZGamma[2.4952])",
   "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",#for 20,30,40
   "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[91.0])",#for 20,30,40
#      "Exponential::backgroundPass(mass, lp[0,-5,5])",#10,15,20
#      "Exponential::backgroundFail(mass, lf[0,-5,5])",#10,15,20

#"RooPolynomial::backgroundPass(mass,{a0[0,-10.,10.],a1[0,-10.,10.]}) ",#10-15-20
#"RooPolynomial::backgroundFail(mass,{a0[0,-10.,10.],a1[0,-10.,10.]}) ",#10-15-2

#  "RooPolynomial::backgroundPass(mass,{a0[0,-10.,10.],a1[0,-10.,10.],a2[0,-10.,10.],a3[0,-10.,10.]}) ",
#  "RooPolynomial::backgroundFail(mass,{a0[0,-10.,10.],a1[0,-10.,10.],a2[0,-10.,10.],a3[0,-10.,10.]}) ",

#    "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[91.0])",#for 20,30,40
#  "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[91.0])",#for 20,30,40
#   "RooChebychev::backgroundPass(mass, {cPass[0,-5.5,5.5], cPass2[0,-5.5,5.5]})", #for 40-500
#   "RooPolynomial::backgroundPass(mass,{a0[0,-10.,10.],a1[0,-10.,10.],a2[0,-10.,10.],a3[0,-10.,10.],a4[0,-10.,10.]}) ",#20,30
 #  "RooPolynomial::backgroundFail(mass,{a0[0,-10.,10.],a1[0,-10.,10.],a2[0,-10.,10.],a3[0,-10.,10.],a4[0,-10.,10.]}) ",#20,30
#   "RooChebychev::backgroundFail(mass, {cFail[0,-5.5,5.5], cFail2[0,-5.5,5.5]})", #for 40-500
      "FCONV::signal(mass, signalPhy, signalRes)",
#      "Exponential::signal2(mass, lp[-0.1,-15.,15.])",
#      "RooLandau::signal2(mass,meanL[0.,-5.0,1.0],sigmaL[2.0, 0.1, 8.0])",
      "SUM::signal(vFrac[0.75,0,1]*signal2, signal1)",
      "efficiency[0.9,0.0,1.0]",
      "signalFractionInPassing[0.999,0.90,1.0]"
    ),
)

bwCBandCMSShapeBackground =cms.PSet(
    fittingFunction= cms.vstring(
      "RooCBShape::signalRes(mass, mean[0.,-5.0,1.0], sigma[2.0, 0.1, 8.0], alpha[1.0,0.2,5.0], n[6.7489,1.0,15.0])",  ### signal resolution for "pass" sample
      "BreitWigner::signalPhy(mass, ZMass[91.1876,80,100], ZGamma[2.4952])", ### Relativistic B-W
      "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",#for 20,30,40
      "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",#for 20,30,40
      "FCONV::signal(mass, signalPhy, signalRes)",
      "efficiency[0.9,0.0,1.0]",
      "signalFractionInPassing[0.999,0.90,1.0]"
    ),
)

## we need a separate PDF distribution for each bin.
## I create a mechanism below which may not be very
## elegant (because of an undelying C++ ugly class)
## but it works.

#    "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
 #   "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",

#    "CBShape::signalRes1(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
 #   "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
# CBExGaussShape::signalResPass(pair_mass, meanP[0.], sigmaP[8.5695e-04, 0., 3.],alphaP[3.8296e-04], nP[6.7489e+00], sigmaP_2[2.5849e+00], fracP[6.5704e-01])",  ### signal resolution for "pass" sample
  #    "CBExGaussShape::signalResFail(pair_mass, meanF[2.0946e-01, -5., 5.], sigmaF[8.5695e-04, 0., 5.],alphaF[3.8296e-04], nF[6.7489e+00], sigmaF_2[2.5849e+00], fracF[6.5704e-01])",  ### signal resolution for "fail" sample     
#for low pt bins, 10,15,20 use Exponential::backgroundPass(mass, lp[0,-5,5]),Exponential::backgroundFail(mass, lf[0,-5,5])   
 
##WARNING: if ur python is not good:  "DONT TOUCH THIS THING !!!"

SYSsurrogate = cms.PSet(
    fittingFunction   = cms.vstring(
   "RooCBShape::signalRes(mass, mean[0.,-5.0,1.0], sigma[2.0, 0.1, 8.0], alpha[1.0,0.2,5.0], n[6.7489,1.0,15.0])",  ### signal resolution for "pass" sample
#    "CBShape::signalRes(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
 #   "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
    "ZGeneratorLineShape::signalPhy1(mass,0,30,2,2.5,1)",
    "ZGeneratorLineShape::signalPhy2(mass,20,30,2,2.5,0)",
"RooPolynomial::backgroundPass(mass,{a0[0,-100.,100.],a1[0,-100000.,100000.]}) ",#10-15-20
"RooPolynomial::backgroundFail(mass,{a0f[0,-100.,100.],a1f[0,-100000.,100000.]}) ",#10-15-2
#"RooPolynomial::backgroundPass(mass,{a0[0,-1000.,1000.],a1[0,-10000.,10000.],a2[0,-10000000.,10000000.]})",
#"RooPolynomial::backgroundFail(mass,{a0f[0,-1000.,1000.],a1f[0,-10000.,10000.],a2f[0,-10000000.,10000000.]})",
#  "RooPolynomial::backgroundPass(mass,{a0[0,-10000.,10000.],a1[0,-10000.,10000.],a2[0,-1000000.,1000000.],a3[0,-1000000.,1000000.]}) ",
 # "RooPolynomial::backgroundFail(mass,{a0f[0,-10000.,10000.],a1f[0,-10000.,10000.],a2f[0,-1000000.,1000000.],a3f[0,-1000000.,1000000.]}) ",
#  "RooPolynomial::backgroundFail(mass,{a0f[0,-1000.,1000.],a1f[0,-1000.,1000.],a2f[0,-1000.,1000.],a3f[0,-1000.,1000.]}) ",
##    "Exponential::backgroundPass(mass, lp[0,-5,5])",#10,15,20
#    "Exponential::backgroundFail(mass, lf[0,-5,5])",#10,15,20
#    "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[91.0])",#for 20,30,40
#    "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[91.0])",#for 20,30,40
#   "RooChebychev::backgroundPass(mass, {cPass[0,-1.5,1.5], cPass2[0,-1.5,1.5]})", #for 40-500
#   "RooChebychev::backgroundFail(mass, {cFail[0,-1.5,1.5], cFail2[0,-1.5,1.5]})", #for 40-500
    "FCONV::signalPass(mass, signalPhy1, signalRes)",
    "FCONV::signalFail(mass, signalPhy2, signalRes)",
    "efficiency[0.9,0,1]",
    "signalFractionInPassing[0.9]"
    ),
    )

def SetSignalShape(pdfPset,ptLow,ptHi,etaLow,etaHi): 
    """
    WARNING: THIS FUNCTION IS EXTREMELY DEPENDENT
    ON HOW THE SURROGATE WAS DEFINED. SO THE
    DEFINITION OF \"surrorgate\" ABOVE SHOULD NOT BE
    ALTERED. THE CODE:
    
    pdfPset = Surrogate PSet
    token2 = {1,2,3,4}, 1= (20<M<30), 2 = (30<M<40), 3=(40<M<50) so on...
    token3 = {0,1}, 0 =  Endcap, 1 = Barrel
    
    """
    
    ##make string for passing signal shape
    passing = "ZGeneratorLineShape::signalPhy1(mass,"+ptLow+","+ptHi+","+etaLow+","+etaHi+",1)"
    ##make string for failing signal shape
    failing = "ZGeneratorLineShape::signalPhy2(mass,"+ptLow+","+ptHi+","+etaLow+","+etaHi+",0)"
        
    ##find index for the passing signal shape
    indexPass = pdfPset.fittingFunction.index("ZGeneratorLineShape::signalPhy1(mass,0,30,2,2.5,1)")
    ##find index for the failing signal shape
    indexFail = pdfPset.fittingFunction.index("ZGeneratorLineShape::signalPhy2(mass,20,30,2,2.5,0)")
    
    ##replace old pass string with new one
    pdfPset.fittingFunction[indexPass] = passing
    ###replace old fail string with new one
    pdfPset.fittingFunction[indexFail] = failing
    print passing, '\n', failing

genConvCBEx = surrogate.clone()
SetSignalShape(genConvCBEx, "40", "50","0","0.08")

#print genConvCBEx
    

genConvCBEx = cms.PSet(
    fittingFunction   = cms.vstring(
        "CBShape::signalRes1(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "ZGeneratorLineShape::signalPhy1(mass,1,1)",
        "ZGeneratorLineShape::signalPhy2(mass,0,1)",
        "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
        "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
        "FCONV::signalPass(mass, signalPhy1, signalRes1)",
        "FCONV::signalFail(mass, signalPhy2, signalRes2)",     
        "efficiency[0.9,0,1]", 
        "signalFractionInPassing[0.9]"
        ),
    )


genConvCBExHiPt= cms.PSet(
    fittingFunction   = cms.vstring(
        "CBShape::signalRes1(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "ZGeneratorLineShape::signalPhy1(mass,1,3)",
        "ZGeneratorLineShape::signalPhy2(mass,0,3)",
        "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
        "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
        "FCONV::signalPass(mass, signalPhy1, signalRes1)",
        "FCONV::signalFail(mass, signalPhy2, signalRes2)",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
        ),
    )


recConvCBEx = cms.PSet(
    fittingFunction   = cms.vstring(
        "CBExGaussShape::signalResPass(mass, meanP[0.], sigmaP[8.5695e-04, 0., 3.],alphaP[3.8296e-04], nP[6.7489e+00], sigmaP_2[2.5849e+00], fracP[6.5704e-01])",
        "CBExGaussShape::signalResFail(mass, meanF[2.0946e-01, -5., 5.], sigmaF[8.5695e-04, 0., 5.],alphaF[3.8296e-04], nF[6.7489e+00], sigmaF_2[2.5849e+00], fracF[6.5704e-01])",
        "ZGeneratorLineShape::passSignalPhy(mass, /afs/hep.wisc.edu/home/kaur/Effi/CMSSW_4_4_2/src/PhysicsTools/TagAndProbe/test/ScToPf.root, hMassPassHi)",
        "ZGeneratorLineShape::failSignalPhy(mass, ScToPf.root, hMassPassLow)",
        "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
        "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
        "FCONV::signalPass(mass, passSignalPhy, signalResPass)",
        "FCONV::signalFail(mass, failSignalPhy, signalResFail)",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
        ),
    )


gausPlusLinear = cms.PSet(
    fittingFunction = cms.vstring(
        "Voigtian::signal1(mass, mean1[90,80,100], width[2.495,-5,5],sigma1[2,-3,3])",
        "Voigtian::signal2(mass, mean2[90,80,100], width[2.495,-5,5],sigma2[4,-10,10])",
        "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
        "Exponential::backgroundPass(mass, lp[-0.1,-1,1])",
        "Exponential::backgroundFail(mass, lf[-0.1,-1,1])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
        ),
    )


ptBin=cms.vdouble(20,30)
print ptBin[0]
