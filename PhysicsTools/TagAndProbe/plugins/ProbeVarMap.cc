
// -*- C++ -*-
//
// Package:    ProbeVarMap
// Class:      ProbeVarMap
// 
/**\class ProbeVarMap ProbeVarMap.cc EGamma/ProbeVarMap/src/ProbeVarMap.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lovedeep Kaur (Panjab U)
//         Created:  Tue Oct 29 09:51:18 CDT 2013
// $Id$
//
//


// system include files
#include <memory>
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
// reco track and vertex 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include <map>
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
//
// class declaration
//

class ProbeVarMap : public edm::EDProducer {
   public:
      explicit ProbeVarMap(const edm::ParameterSet&);
      ~ProbeVarMap();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
        edm::InputTag probes_;            
  std::vector<edm::InputTag> inputTagIsoDepElectrons_;
  std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;
edm::InputTag rhoIsoInputTag_;

};

//typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
//typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;
//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
ProbeVarMap::ProbeVarMap(const edm::ParameterSet& iConfig):
  probes_(iConfig.getParameter<edm::InputTag>("probes")),
  inputTagIsoDepElectrons_(iConfig.getParameter< std::vector<edm::InputTag> >("IsoDepElectron")),
  inputTagIsoValElectronsPFId_(iConfig.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF")),
  rhoIsoInputTag_(iConfig.getParameter<edm::InputTag>("rhoIsoInputTag"))
    //isoValInputTags_(iConfig.getParameter<std::vector<edm::InputTag> >("IsoValElectronPF"))
{
//produces<edm::ValueMap<float> >();
produces<edm::ValueMap<float> >("getd0vtx");//.setBranchAlias("d0vtxs");
produces<edm::ValueMap<float> >("getdzvtx");//.setBranchAlias("dzvtxs");
produces<edm::ValueMap<float> >("getiso");
produces<edm::ValueMap<float> >("getisonh");
produces<edm::ValueMap<float> >("getisoem");
produces<edm::ValueMap<float> >("getisoch");
produces<edm::ValueMap<float> >("getvtxFitConversion");

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


ProbeVarMap::~ProbeVarMap()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ProbeVarMap::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
    // read input
    Handle<View<reco::GsfElectron> > probes;
    iEvent.getByLabel(probes_,  probes);
  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel("gsfElectrons", electrons);


    std::vector<float> d0vtxs;
    std::vector<float> dzvtxs;
    std::vector<float> viso;
    std::vector<float> viso_em;
    std::vector<float> viso_nh;
    std::vector<float> viso_ch;
    std::vector<float> vgetvtxFitConversion;     // fill
    // vertices
    edm::Handle<reco::VertexCollection> vtxs;
    iEvent.getByLabel("offlinePrimaryVertices", vtxs);

    // impact parameter variables
   float d0vtx         = 0.0;
   float dzvtx         = 0.0;
    View<reco::GsfElectron>::const_iterator probe, endprobes = probes->end();
    for (probe = probes->begin(); probe != endprobes; ++probe) {
  if (vtxs->size() > 0) {
        reco::VertexRef vtx(vtxs, 0);    
        d0vtx = probe->gsfTrack()->dxy(vtx->position());
        dzvtx = probe->gsfTrack()->dz(vtx->position());
    } else {
        d0vtx = probe->gsfTrack()->dxy();
        dzvtx = probe->gsfTrack()->dz();
    }
//std::cout<<d0vtx<<'\t'<<dzvtx<<std::endl;
    d0vtxs.push_back(d0vtx);
    dzvtxs.push_back(dzvtx);

}
    

    // conversions
    edm::Handle<reco::ConversionCollection> conversions_h;
    iEvent.getByLabel("allConversions", conversions_h);
//    const edm::Handle<reco::ConversionCollection> &conversions,

    // beam spot
    edm::Handle<reco::BeamSpot> beamspot_h;
    iEvent.getByLabel("offlineBeamSpot", beamspot_h);
    const reco::BeamSpot &beamSpot = *(beamspot_h.product());
//    const reco::BeamSpot &beamspot,

    //std::cout<<vtxFitConversion<<std::endl;
float iso=9999.0,iso_ch=9999.0,iso_em=9999.0,iso_nh=9999.0,pt=0.0,vtxFitConversion=0;
    // rho for isolation
    edm::Handle<double> rhoIso_h;
    iEvent.getByLabel(rhoIsoInputTag_, rhoIso_h);
    double rhoIso = *(rhoIso_h.product());
    
  // get the iso deposits. 3 (charged hadrons, photons, neutral hadrons)
  unsigned nTypes=3;
  IsoDepositMaps electronIsoDep(nTypes);

  for (size_t j = 0; j<inputTagIsoDepElectrons_.size(); ++j) {
    iEvent.getByLabel(inputTagIsoDepElectrons_[j], electronIsoDep[j]);
  }
  IsoDepositVals electronIsoValPFId(nTypes);
  const IsoDepositVals * electronIsoVals = &electronIsoValPFId;

  for (size_t j = 0; j<inputTagIsoValElectronsPFId_.size(); ++j) {
    iEvent.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoValPFId[j]);
  }


  for(edm::View<reco::GsfElectron>::const_iterator gsfIt = probes->begin();
      gsfIt != probes->end(); ++gsfIt){
//iso=0.0,iso_ch=0.0,iso_em=0.0,iso_nh=0.0;
  // Loop over electrons
    unsigned nele=electrons->size();
   for(unsigned iele=0; iele<nele;++iele) {
      reco::GsfElectronRef myElectronRef(electrons,iele);
   if (fabs(gsfIt->p4().Pt() - myElectronRef->pt()) < 0.0001 ){

    vtxFitConversion = ConversionTools::hasMatchedConversion(*myElectronRef, conversions_h, beamSpot.position());

      iso_ch =  (*(*electronIsoVals)[0])[myElectronRef];
      iso_em = (*(*electronIsoVals)[1])[myElectronRef];
      iso_nh = (*(*electronIsoVals)[2])[myElectronRef];
		float eta    = myElectronRef->superCluster()->eta();
 float AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, eta, ElectronEffectiveArea::kEleEAData2012);
    // apply to neutrals
    double rhoPrime = std::max(rhoIso, 0.0);
    double iso_n = std::max(iso_nh + iso_em - rhoPrime * AEff, 0.0);
    pt            = myElectronRef->pt();
    
    // compute final isolation
        iso = (iso_n + iso_ch) / pt;
}}
vgetvtxFitConversion.push_back(vtxFitConversion);
viso.push_back(iso);
viso_ch.push_back(iso_ch);
viso_nh.push_back(iso_nh);
viso_em.push_back(iso_em);
}
/*    View<reco::Candidate>::const_iterator probe, endprobes = probes->end();
    for (probe = probes->begin(); probe != endprobes; ++probe) {
    double pt = probe->pt();
    values.push_back(pt);
    }
*/
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    ValueMap<float>::Filler filler(*valMap);
    filler.insert(probes, d0vtxs.begin(), d0vtxs.end());
    filler.fill();
    iEvent.put(valMap,"getd0vtx");

    std::auto_ptr<ValueMap<float> > dzvalMap(new ValueMap<float>());
    ValueMap<float>::Filler dzfiller(*dzvalMap);
    dzfiller.insert(probes, dzvtxs.begin(), dzvtxs.end());
    dzfiller.fill();
    iEvent.put(dzvalMap,"getdzvtx");

    std::auto_ptr<ValueMap<float> > isoMap(new ValueMap<float>());
    ValueMap<float>::Filler isofiller(*isoMap);
    isofiller.insert(probes, viso.begin(), viso.end());
    isofiller.fill();
    iEvent.put(isoMap,"getiso");

    std::auto_ptr<ValueMap<float> > iso_chMap(new ValueMap<float>());
    ValueMap<float>::Filler iso_chfiller(*iso_chMap);
    iso_chfiller.insert(probes, viso_ch.begin(), viso_ch.end());
    iso_chfiller.fill();
    iEvent.put(iso_chMap,"getisoch");

    std::auto_ptr<ValueMap<float> > iso_nhMap(new ValueMap<float>());
    ValueMap<float>::Filler iso_nhfiller(*iso_nhMap);
    iso_nhfiller.insert(probes, viso_nh.begin(), viso_nh.end());
    iso_nhfiller.fill();
    iEvent.put(iso_nhMap,"getisonh");

    std::auto_ptr<ValueMap<float> > iso_emMap(new ValueMap<float>());
    ValueMap<float>::Filler iso_emfiller(*iso_emMap);
    iso_emfiller.insert(probes, viso_em.begin(), viso_em.end());
    iso_emfiller.fill();
    iEvent.put(iso_emMap,"getisoem");

  std::auto_ptr<ValueMap<float> > vtxFitConversionMap(new ValueMap<float>());
    ValueMap<float>::Filler vtxFitConversionfiller(*vtxFitConversionMap);
    vtxFitConversionfiller.insert(probes, vgetvtxFitConversion.begin(), vgetvtxFitConversion.end());
    vtxFitConversionfiller.fill();
    iEvent.put(vtxFitConversionMap,"getvtxFitConversion");

}

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 


// ------------ method called once each job just before starting event loop  ------------
void 
ProbeVarMap::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ProbeVarMap::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
ProbeVarMap::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ProbeVarMap::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ProbeVarMap::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ProbeVarMap::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ProbeVarMap::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProbeVarMap);
