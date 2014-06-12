// -*- C++ -*-
//
// Package:    PileupWeightProducer
// Class:      PileupWeightProducer
// 
/**\class PileupWeightProducer PileupWeightProducer.cc 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ricardo Vasquez Sierra,6 R-025,+41227672274,
//         Created:  Mon Nov 21 15:05:26 CET 2011
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include <vector>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//
// class declaration
//

class PileupWeightProducer : public edm::EDProducer {
   public:
      explicit PileupWeightProducer(const edm::ParameterSet&);
      ~PileupWeightProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  bool firsttime_;
  std::string pileupMC_;
  std::string pileupData_;
  edm::LumiReWeighting LumiWeights_;
  edm::Lumi3DReWeighting LumiWeightsNominal_;
  edm::Lumi3DReWeighting LumiWeightsUp_;
  edm::Lumi3DReWeighting LumiWeightsDown_;
std::vector< float > Data2012_;
std::vector<float> Summer2012_S10_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
PileupWeightProducer::PileupWeightProducer(const edm::ParameterSet& iConfig)
{
//   firsttime_= iConfig.existsAs<bool>("FirstTime") ? iConfig.getParameter<bool>("FirstTime") : true ;
//   pileupMC_ = iConfig.existsAs<std::string>("PileupMCFile") ? iConfig.getParameter<std::string>("PileupMCFile") : "PUMC_dist.root" ;
//   pileupData_ = iConfig.existsAs<std::string>("PileupDataFile") ? iConfig.getParameter<std::string>("PileupDataFile") : "PUData_dist.root" ;

  firsttime_ =  iConfig.getUntrackedParameter<bool>("FirstTime");

  //register your products
  
  produces<std::vector<float> >( "pileupWeights" ).setBranchAlias( "pileupWeights" );
  
  

  if ( firsttime_ )
    {
//      std::cout<< " Initializing with the following files MC: " << pileupMC_ << " data: " << pileupData_ << std::endl;
Double_t Summer2012_S10[60] = {
                         2.560E-06,
                         5.239E-06,
                         1.420E-05,
                         5.005E-05,
                         1.001E-04,
                         2.705E-04,
                         1.999E-03,
                         6.097E-03,
                         1.046E-02,
                         1.383E-02,
                         1.685E-02,
                         2.055E-02,
                         2.572E-02,
                         3.262E-02,
                         4.121E-02,
                         4.977E-02,
                         5.539E-02,
                         5.725E-02,
                         5.607E-02,
                         5.312E-02,
                         5.008E-02,
                         4.763E-02,
                         4.558E-02,
                         4.363E-02,
                         4.159E-02,
                         3.933E-02,
                         3.681E-02,
                         3.406E-02,
                         3.116E-02,
                         2.818E-02,
                         2.519E-02,
                         2.226E-02,
                         1.946E-02,
                         1.682E-02,
                         1.437E-02,
                         1.215E-02,
                         1.016E-02,
                         8.400E-03,
                         6.873E-03,
                         5.564E-03,
                         4.457E-03,
                         3.533E-03,
                         2.772E-03,
                         2.154E-03,
                         1.656E-03,
                         1.261E-03,
                         9.513E-04,
                         7.107E-04,
                         5.259E-04,
                         3.856E-04,
                         2.801E-04,
                         2.017E-04,
                         1.439E-04,
                         1.017E-04,
                         7.126E-05,
                         4.948E-05,
                         3.405E-05,
                         2.322E-05,
                         1.570E-05,
                         5.005E-06}; //from https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios#2012_Pileup_Scenario_s
double Data2012[60]={
12260.8,32850.4,92330.3,339464,618478,3.0497e+06,1.77215e+07,5.41421e+07,1.30521e+08,2.58981e+08,4.46344e+08,6.8564e+08,8.81642e+08,9.99085e+08,1.07862e+09,1.13797e+09,1.17211e+09,1.18207e+09,1.17701e+09,1.16108e+09,1.13609e+09,1.10481e+09,1.06807e+09,1.02107e+09,9.55582e+08,8.6706e+08,7.58729e+08,6.38851e+08,5.16436e+08,3.99862e+08,2.96257e+08,2.10055e+08,1.42404e+08,9.20546e+07,5.65387e+07,3.29089e+07,1.815e+07,9.51188e+06,4.76417e+06,2.29967e+06,1.08138e+06,501998,233744,111112,54826,28402.3,15490.1,8845.44,5236.34,3180.14,1964.06,1225.15,767.779,481.279,300.644,186.558,114.687,69.6938,41.7929,24.6979
};



 for( int i=0; i<60; ++i) {
      Data2012_.push_back(Data2012[i]);
      Summer2012_S10_.push_back(Summer2012_S10[i]);
   }
 LumiWeights_ = edm::LumiReWeighting(Summer2012_S10_,Data2012_);
}
//      LumiWeights_ = edm::LumiReWeighting(pileupMC_, pileupData_, "pileup", "pileup");}
/*      LumiWeightsNominal_.weight3D_set( pileupMC_, pileupData_, "pileup", "pileup");
      LumiWeightsUp_.weight3D_set( pileupMC_, pileupData_, "pileup", "pileup");
      LumiWeightsDown_.weight3D_set( pileupMC_, pileupData_, "pileup", "pileup");

      LumiWeightsNominal_.weight3D_init(1.0);
      LumiWeightsUp_.weight3D_init(1.08);
      LumiWeightsDown_.weight3D_init(0.92);
    }
  else 
    {
      std::cout<< " Initializing with Weight3D.root " << std::endl; 
      LumiWeightsNominal_.weight3D_init("Weight3D.root");
      LumiWeightsUp_.weight3D_init("Weight3DscaleUp.root");
      LumiWeightsDown_.weight3D_init("Weight3DscaleDown.root");
    }
*/


  
}


PileupWeightProducer::~PileupWeightProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called to produce the data  ------------
void
PileupWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   std::auto_ptr<std::vector<float> > pileupWeights( new std::vector<float> );

/*   edm::EventBase* iEventB = dynamic_cast<edm::EventBase*>(&iEvent);
  double MyWeight = LumiWeights_.weight( (*iEventB) );

   double nominalWeight3D = LumiWeightsNominal_.weight3D( (*iEventB) );
   double weight3DUp = LumiWeightsUp_.weight3D( (*iEventB) );
   double weight3DDown = LumiWeightsDown_.weight3D( (*iEventB) );

   pileupWeights->push_back( MyWeight );

   pileupWeights->push_back( nominalWeight3D );
   pileupWeights->push_back( weight3DUp );
   pileupWeights->push_back( weight3DDown );

   iEvent.put(pileupWeights, "pileupWeights");
*/   
///yar kuch kerna paina...mera apna style he theek e
Handle<std::vector< PileupSummaryInfo > >  PupInfo;
iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

std::vector<PileupSummaryInfo>::const_iterator PVI;

float Tnpv = -1;
for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

   int BX = PVI->getBunchCrossing();

   if(BX == 0) { 
     Tnpv = PVI->getTrueNumInteractions();
     continue;
   }

}
double MyWeight = LumiWeights_.weight( Tnpv );
   pileupWeights->push_back( MyWeight );
   iEvent.put(pileupWeights, "pileupWeights");

}

// ------------ method called once each job just before starting event loop  ------------
void 
PileupWeightProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PileupWeightProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PileupWeightProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PileupWeightProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PileupWeightProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PileupWeightProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PileupWeightProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PileupWeightProducer);
