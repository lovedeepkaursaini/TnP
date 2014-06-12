// -*- C++ -*-
//
// Package:    EGammaCutBasedEleIdAnalyzer
// Class:      EGammaCutBasedEleIdAnalyzer
// 
/**\class EGammaCutBasedEleIdAnalyzer EGammaCutBasedEleIdAnalyzer.cc EGamma/EGammaCutBasedEleIdAnalyzer/src/EGammaCutBasedEleIdAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Dave Evans,510 1-015,+41227679496,
//         Created:  Tue Apr 10 11:17:29 CEST 2012
// $Id: EGammaCutBasedEleIdAnalyzer.cc,v 1.2 2012/04/11 15:24:16 dlevans Exp $
//
//
using namespace std;

// system include files
#include <memory>

// user include files
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include <algorithm>
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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include <TFile.h>
#include <TH1F.h>
#include "TTree.h"
//
// class declaration
//

class EGammaCutBasedEleIdAnalyzer : public edm::EDAnalyzer {
    public:
        EGammaCutBasedEleIdAnalyzer(const edm::ParameterSet& cfg);
        ~EGammaCutBasedEleIdAnalyzer();

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
        TH1F *h1_pt_veto_;
        TH1F *h1_pt_loose_;
        TH1F *h1_pt_medium_;
        TH1F *h1_pt_tight_;
        TH1F *h1_pt_trig_;
        TH1F *h1_pt_fbremeopin_;
TTree *tree_;
  int nEle_;
vector<int>    eleCharge_;
};

//
// constants, enums and typedefs
//

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;

//
// static data member definitions
//

//
// constructors and destructor
//
EGammaCutBasedEleIdAnalyzer::EGammaCutBasedEleIdAnalyzer(const edm::ParameterSet& iConfig)
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
    h1_pt_veto_          = fs->make<TH1F>("h1_pt_veto",          "pt (veto)",       100, 0.0, 100.0);
    h1_pt_loose_         = fs->make<TH1F>("h1_pt_loose",         "pt (loose)",      100, 0.0, 100.0);
    h1_pt_medium_        = fs->make<TH1F>("h1_pt_medium",        "pt (medium)",     100, 0.0, 100.0);
    h1_pt_tight_         = fs->make<TH1F>("h1_pt_tight",         "pt (tight)",      100, 0.0, 100.0);
    h1_pt_trig_          = fs->make<TH1F>("h1_pt_trig",          "pt (trig)",       100, 0.0, 100.0); 
    h1_pt_fbremeopin_    = fs->make<TH1F>("h1_pt_fbremeopin",    "pt (fbremeopin)", 100, 0.0, 100.0);
//  tree_=tree;
   tree_ = fs->make<TTree>("tree","tree"); 
  tree_->Branch("nEle", &nEle_, "nEle/I");
  tree_->Branch("eleCharge", &eleCharge_);
}


EGammaCutBasedEleIdAnalyzer::~EGammaCutBasedEleIdAnalyzer()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
  eleCharge_.clear();
//  delete tree_;
}


//
// member functions
//

// ------------ method called for each event  ------------
    void
EGammaCutBasedEleIdAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nEle_ = 0;
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

    edm::Handle<edm::View<reco::PFCandidate> > pfCandidates;
  //event.getByLabel("particleFlow", pfCandidates);
    iEvent.getByLabel("pfNoPileUpIso",pfCandidates);


    // loop on electrons
    unsigned int n = els_h->size();
  nEle_=n;
  for(unsigned int i = 0; i < n; ++i) {

        // get reference to electron
        reco::GsfElectronRef ele(els_h, i);

   /*     float sumPt1(0.0);
        float sumPtChHad(0.0);float sumPtChEM(0.0);float sumPtNtEM(0.0); float sumPtMuon(0.0);float sumPtNtHad(0.0);

        for(unsigned ipf=0;ipf<pfCandidates->size();ipf++) {
          float dR1 = deltaR(ele->eta(),ele->phi(),(*pfCandidates)[ipf].eta(),(*pfCandidates)[ipf].phi());
         if (dR1 < 0.3 ) {
         std::cout<<"[ PFCandidate: Type: "<< (*pfCandidates)[ipf].particleId() <<"  et:  "<<(*pfCandidates)[ipf].et() <<"  Pt:  "<<(*pfCandidates)[ipf].pt()<<
             "  Eta: "<<(*pfCandidates)[ipf].eta()<<"  Phi: "<<(*pfCandidates)[ipf].phi()<<"],   [  delR: "<<dR1<<"]"<<std::endl;

//         sumPt1 += (*pfCandidates)[ipf].pt();
         if((*pfCandidates)[ipf].particleId() == 1)  sumPtChHad += (*pfCandidates)[ipf].pt();
         if((*pfCandidates)[ipf].particleId() == 2)  sumPtChEM += (*pfCandidates)[ipf].pt();
         if((*pfCandidates)[ipf].particleId() == 3)  sumPtMuon += (*pfCandidates)[ipf].pt();
         if((*pfCandidates)[ipf].particleId() == 4)  sumPtNtEM += (*pfCandidates)[ipf].pt();
         if((*pfCandidates)[ipf].particleId() == 5)  sumPtNtHad += (*pfCandidates)[ipf].pt();

//         std::cout<<"SumPtAll (electron not included):"<<sumPt1-ele->pt()<<std::endl;
         std::cout<<"SumPtCharHad :"<<sumPtChHad<<std::endl;
         std::cout<<"SumPtCharEM :"<<sumPtChEM<<std::endl;
         std::cout<<"SumPtMuon :"<<sumPtMuon<<std::endl;
         std::cout<<"SumPtNtEM :"<<sumPtNtEM<<std::endl;
         std::cout<<"SumPtNtHad :"<<sumPtNtHad<<std::endl;

         }


}*/
        //
        // get particle flow isolation
        //

        double iso_ch =  (*(isoVals)[0])[ele];
        double iso_em = (*(isoVals)[1])[ele];
        double iso_nh = (*(isoVals)[2])[ele];
//lvdp
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
 for(unsigned int i=0; i<vtx_h->size();i++)
 {
        reco::VertexRef vtx1(vtx_h, i);
//        std::cout<<i<<'\t'<<vtx1->position()<<'\t'<<vtx1->x()<<", isvalid: "<<vtx1->isValid()<<std::endl;
 }
}

    if (vtx_h->size() > 0) {
        reco::VertexRef vtx(vtx_h, 0);
//cout<<vtx->position()<<endl;
        d0vtx = ele->gsfTrack()->dxy(vtx->position());
        dzvtx = ele->gsfTrack()->dz(vtx->position());
    } else {
        d0vtx = ele->gsfTrack()->dxy();
        dzvtx = ele->gsfTrack()->dz();
    }
//cout<<d0vtx<<'\t'<<ele->gsfTrack()->pt()<<endl;
    // conversion rejection variables
    bool vtxFitConversion = ConversionTools::hasMatchedConversion(*ele, conversions_h, beamSpot.position());
    float mHits = ele->gsfTrack()->trackerExpectedHitsInner().numberOfHits();




    unsigned int idx = isEB ? 0 : 1;
double EffectiveArea=0.0;
    // effective area for isolation
          if (abs(eta) >= 0.0 && abs(eta) < 1.0 ) EffectiveArea = 0.130;
          if (abs(eta) >= 1.0 && abs(eta) < 1.479 ) EffectiveArea = 0.137;
          if (abs(eta) >= 1.479 && abs(eta) < 2.0 ) EffectiveArea = 0.067;
          if (abs(eta) >= 2.0 && abs(eta) < 2.2 ) EffectiveArea = 0.089;
          if (abs(eta) >= 2.2 && abs(eta) < 2.3 ) EffectiveArea = 0.107;
          if (abs(eta) >= 2.3 && abs(eta) < 2.4 ) EffectiveArea = 0.110;
          if (abs(eta) >= 2.4) EffectiveArea = 0.138;
//    float AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, eta, ElectronEffectiveArea::kEleEAData2012);

    // apply to neutrals
    double rhoPrime = std::max(rhoIso, 0.0);
    double iso_n = std::max(iso_nh + iso_em - rhoPrime * EffectiveArea, 0.0);

    // compute final isolation
    double iso = (iso_n + iso_ch) / pt;
//cout<< "ana:   hurrey :) isolation is: "<<iso <<" AEff, rhoPrime, iso_nh, iso_em (iso_n), iso_ch, iso:  "<< EffectiveArea<<'\t'<<rhoPrime<<'\t'<<iso_nh<<'\t'<<iso_em<<" ( "<<iso_n<<" ) "<<iso_ch<<'\t'<<iso<<endl;


//endlvdp
        //
        // test ID
        //
        // working points
        bool veto       = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
        bool loose      = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
        bool medium     = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
        bool tight      = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, ele, conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);

        // eop/fbrem cuts for extra tight ID
        bool fbremeopin = EgammaCutBasedEleId::PassEoverPCuts(ele);

        // cuts to match tight trigger requirements
        bool trigtight = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, ele);

        // for 2011 WP70 trigger
        bool trigwp70 = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERWP70, ele);

        //
        // fill histograms
        //

        h1_pt_->Fill(ele->pt());
        if (veto)       h1_pt_veto_         ->Fill(ele->pt());
        if (loose)      h1_pt_loose_        ->Fill(ele->pt());
        if (medium)     h1_pt_medium_       ->Fill(ele->pt());
        if (tight)      h1_pt_tight_        ->Fill(ele->pt());
        if (trigtight)  h1_pt_trig_         ->Fill(ele->pt());
        if (fbremeopin) h1_pt_fbremeopin_   ->Fill(ele->pt());
      eleCharge_ .push_back(ele->charge());
        //
        // print decisions
        //
        if (printDebug_) {
            printf("Run: %u , lumi: %u , event: %u ",       iEvent.id().run(), iEvent.luminosityBlock(), iEvent.id().event());
            printf("has electron with pt: %f , eta: %f, phi: %f ,",ele->pt(),ele->superCluster()->eta(),ele->phi());
        //    printf("veto(%i), ",        veto);
          //  printf("loose(%i), ",       loose);
            printf(" pass medium(%i), ",      medium);
            printf("pass tight(%i)\n ",       tight);
            //printf("trigtight(%i), ",   trigtight);
//            printf("trigwp70(%i), ",    trigwp70);
  //          printf("fbremeopin(%i)\n",  fbremeopin);
        }
  tree_->Fill();
    }

}


// ------------ method called once each job just before starting event loop  ------------
    void 
EGammaCutBasedEleIdAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
    void 
EGammaCutBasedEleIdAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
    void 
EGammaCutBasedEleIdAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
    void 
EGammaCutBasedEleIdAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
    void 
EGammaCutBasedEleIdAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
    void 
EGammaCutBasedEleIdAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EGammaCutBasedEleIdAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGammaCutBasedEleIdAnalyzer);
