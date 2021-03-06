// -*- C++ -*-
//
// Package:    ElectronTests
// Class:      ElectronTests
// 
/**\class ElectronTests ElectronTests.cc PhysicsTools/ElectronTests/src/ElectronTests.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lovedeep Kaur (Panjab U)
//         Created:  Thu Sep  5 09:39:49 CDT 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include <vector>

#include <TFile.h>
#include <TH1F.h>

//
// class declaration
//

class ElectronTests : public edm::EDAnalyzer {
   public:
      explicit ElectronTests(const edm::ParameterSet&);
      ~ElectronTests();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  // input tags
  edm::InputTag               electronsInputTag_;
  edm::InputTag               conversionsInputTag_;
  edm::InputTag               beamSpotInputTag_;
  edm::InputTag               rhoIsoInputTag;
  edm::InputTag               primaryVertexInputTag_;
  std::vector<edm::InputTag>  isoValInputTags_;
        // debug
        bool printDebug_;

        // histograms
        TH1F *h1_pt_;
        TH1F *h1_pt_medium_;
        TH1F *h1_pt_tight_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;

//
// constructors and destructor
//
ElectronTests::ElectronTests(const edm::ParameterSet& iConfig)

{
    // get input parameters
    electronsInputTag_      = iConfig.getParameter<edm::InputTag>("electronsInputTag");
    conversionsInputTag_    = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
    beamSpotInputTag_       = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
    rhoIsoInputTag          = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
    primaryVertexInputTag_  = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    isoValInputTags_        = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");


    // debug
    printDebug_             = iConfig.getParameter<bool>("printDebug");

    // output histograms
    edm::Service<TFileService> fs;

    h1_pt_               = fs->make<TH1F>("h1_pt",               "pt",              100, 0.0, 100.0);
    h1_pt_medium_        = fs->make<TH1F>("h1_pt_medium",        "pt (medium)",     100, 0.0, 100.0);
    h1_pt_tight_         = fs->make<TH1F>("h1_pt_tight",         "pt (tight)",      100, 0.0, 100.0);
   //now do what ever initialization is needed

}


ElectronTests::~ElectronTests()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronTests::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   // electrons
   edm::Handle<reco::GsfElectronCollection> els_h;
   iEvent.getByLabel(electronsInputTag_, els_h);

   // conversions
   edm::Handle<reco::ConversionCollection> conversions_h;
   iEvent.getByLabel(conversionsInputTag_, conversions_h);

   // iso deposits
   IsoDepositVals isoVals(isoValInputTags_.size());
   for (size_t j = 0; j < isoValInputTags_.size(); ++j) {
     iEvent.getByLabel(isoValInputTags_[j], isoVals[j]);
   }

   // beam spot
   edm::Handle<reco::BeamSpot> beamspot_h;
   iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
   const reco::BeamSpot &beamSpot = *(beamspot_h.product());

   // vertices
   edm::Handle<reco::VertexCollection> vtx_h;
   iEvent.getByLabel(primaryVertexInputTag_, vtx_h);

   // rho for isolation
   edm::Handle<double> rhoIso_h;
   iEvent.getByLabel(rhoIsoInputTag, rhoIso_h);
   double rhoIso = *(rhoIso_h.product());

 
   // loop on electrons
   unsigned int n = els_h->size();
   for(unsigned int i = 0; i < n; ++i) {

     // get reference to electron
     reco::GsfElectronRef ele(els_h, i);

     //
     // get particle flow isolation
     //

     double iso_ch =  (*(isoVals)[0])[ele];
     double iso_em = (*(isoVals)[1])[ele];
     double iso_nh = (*(isoVals)[2])[ele];

        bool medium     = 0;//EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
        bool tight      = 1;//EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
        h1_pt_->Fill(ele->pt());
        if (medium)     h1_pt_medium_       ->Fill(ele->pt());
        if (tight)      h1_pt_tight_        ->Fill(ele->pt());
        if (printDebug_) {
            printf("%u %u %u : ",       iEvent.id().run(), iEvent.luminosityBlock(), iEvent.id().event());
            printf("medium(%i), ",      medium);
            printf("tight(%i), ",       tight);
        }
     // get the ID variables from the electron object

     // kinematic variables
     bool isEB           = ele->isEB() ? true : false;
     float pt            = ele->pt();
     float eta           = ele->superCluster()->eta();

     // id variables
     float dEtaIn        = ele->deltaEtaSuperClusterTrackAtVtx();
     float dPhiIn        = ele->deltaPhiSuperClusterTrackAtVtx();
     float sigmaIEtaIEta = ele->sigmaIetaIeta();
     float hoe           = ele->hadronicOverEm();
     float ooemoop       = (1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy());

     // impact parameter variables
     float d0vtx         = 0.0;
     float dzvtx         = 0.0;
     if (vtx_h->size() > 0) {
       reco::VertexRef vtx(vtx_h, 0);    
       d0vtx = ele->gsfTrack()->dxy(vtx->position());
       dzvtx = ele->gsfTrack()->dz(vtx->position());
     } else {
       d0vtx = ele->gsfTrack()->dxy();
       dzvtx = ele->gsfTrack()->dz();
     }

     // conversion rejection variables
     bool vtxFitConversion = ConversionTools::hasMatchedConversion(*ele, conversions_h, beamSpot.position());
     float mHits = ele->gsfTrack()->trackerExpectedHitsInner().numberOfHits(); 

   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElectronTests::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronTests::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ElectronTests::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ElectronTests::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ElectronTests::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ElectronTests::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronTests::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronTests);
