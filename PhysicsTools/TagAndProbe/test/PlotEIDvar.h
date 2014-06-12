//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov  5 14:32:35 2013 by ROOT version 5.28/00c
// from TTree fitter_tree/fitter_tree
// found on file: TnPTree_MC.root
//////////////////////////////////////////////////////////

#ifndef PlotEIDvar_h
#define PlotEIDvar_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class PlotEIDvar {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         probe_gsfEle_abseta;
   Float_t         probe_gsfEle_dEtaIn;
   Float_t         probe_gsfEle_dPhiIn;
   Float_t         probe_gsfEle_e;
   Float_t         probe_gsfEle_et;
   Float_t         probe_gsfEle_eta;
   Float_t         probe_gsfEle_hoe;
   Float_t         probe_gsfEle_mHits;
   Float_t         probe_gsfEle_ooemoop;
   Float_t         probe_gsfEle_pt;
   Float_t         probe_gsfEle_q;
   Float_t         probe_gsfEle_reltrackiso;
   Float_t         probe_gsfEle_sigmaIEtaIEta;
   Float_t         probe_gsfEle_trackiso;
   Float_t         probe_sc_abseta;
   Float_t         probe_sc_energy;
   Float_t         probe_sc_et;
   Float_t         probe_sc_eta;
   Float_t         probe_d0vtx;
   Float_t         probe_dRTau;
   Float_t         probe_dzvtx;
   Float_t         probe_iso;
   Float_t         probe_iso_ch;
   Float_t         probe_iso_em;
   Float_t         probe_iso_nh;
   Float_t         probe_vtxFitConversion;
   Int_t           probe_isWPMedium;
   Float_t         PUweight;
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          event;
   Int_t           event_nPV;
   Float_t         event_met_calomet;
   Float_t         event_met_calosumet;
   Float_t         event_met_calometsignificance;
   Float_t         event_met_tcmet;
   Float_t         event_met_tcsumet;
   Float_t         event_met_tcmetsignificance;
   Float_t         event_met_pfmet;
   Float_t         event_met_pfsumet;
   Float_t         event_met_pfmetsignificance;
   Float_t         event_PrimaryVertex_x;
   Float_t         event_PrimaryVertex_y;
   Float_t         event_PrimaryVertex_z;
   Float_t         event_BeamSpot_x;
   Float_t         event_BeamSpot_y;
   Float_t         event_BeamSpot_z;
   Float_t         mass;
   Int_t           mcTrue;
   Float_t         mcMass;
   Float_t         tag_gsfEle_abseta;
   Float_t         tag_gsfEle_e;
   Float_t         tag_gsfEle_et;
   Float_t         tag_gsfEle_eta;
   Float_t         tag_gsfEle_pt;
   Float_t         tag_gsfEle_q;
   Float_t         tag_sc_abseta;
   Float_t         tag_sc_energy;
   Float_t         tag_sc_et;
   Float_t         tag_sc_eta;
   Float_t         tag_dRTau;
   Int_t           tag_isWPMedium;
   Int_t           tag_passingGsf;
   Float_t         mc_probe_abseta;
   Float_t         mc_probe_e;
   Float_t         mc_probe_et;
   Float_t         mc_probe_eta;
   Float_t         mc_probe_mass;
   Float_t         mc_probe_pt;
   Int_t           mc_probe_flag;
   Float_t         pair_abseta;
   Float_t         pair_eta;
   Float_t         pair_mass;
   Float_t         pair_pt;
   Int_t           pair_mass60to120;

   // List of branches
   TBranch        *b_probe_gsfEle_abseta;   //!
   TBranch        *b_probe_gsfEle_dEtaIn;   //!
   TBranch        *b_probe_gsfEle_dPhiIn;   //!
   TBranch        *b_probe_gsfEle_e;   //!
   TBranch        *b_probe_gsfEle_et;   //!
   TBranch        *b_probe_gsfEle_eta;   //!
   TBranch        *b_probe_gsfEle_hoe;   //!
   TBranch        *b_probe_gsfEle_mHits;   //!
   TBranch        *b_probe_gsfEle_ooemoop;   //!
   TBranch        *b_probe_gsfEle_pt;   //!
   TBranch        *b_probe_gsfEle_q;   //!
   TBranch        *b_probe_gsfEle_reltrackiso;   //!
   TBranch        *b_probe_gsfEle_sigmaIEtaIEta;   //!
   TBranch        *b_probe_gsfEle_trackiso;   //!
   TBranch        *b_probe_sc_abseta;   //!
   TBranch        *b_probe_sc_energy;   //!
   TBranch        *b_probe_sc_et;   //!
   TBranch        *b_probe_sc_eta;   //!
   TBranch        *b_probe_d0vtx;   //!
   TBranch        *b_probe_dRTau;   //!
   TBranch        *b_probe_dzvtx;   //!
   TBranch        *b_probe_iso;   //!
   TBranch        *b_probe_iso_ch;   //!
   TBranch        *b_probe_iso_em;   //!
   TBranch        *b_probe_iso_nh;   //!
   TBranch        *b_probe_vtxFitConversion;   //!
   TBranch        *b_probe_isWPMedium;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_mNPV;   //!
   TBranch        *b_mMET;   //!
   TBranch        *b_mSumET;   //!
   TBranch        *b_mMETSign;   //!
   TBranch        *b_mtcMET;   //!
   TBranch        *b_mtcSumET;   //!
   TBranch        *b_mtcMETSign;   //!
   TBranch        *b_mpfMET;   //!
   TBranch        *b_mpfSumET;   //!
   TBranch        *b_mpfMETSign;   //!
   TBranch        *b_mPVx;   //!
   TBranch        *b_mPVy;   //!
   TBranch        *b_mPVz;   //!
   TBranch        *b_mBSx;   //!
   TBranch        *b_mBSy;   //!
   TBranch        *b_mBSz;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_mcTrue;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_tag_gsfEle_abseta;   //!
   TBranch        *b_tag_gsfEle_e;   //!
   TBranch        *b_tag_gsfEle_et;   //!
   TBranch        *b_tag_gsfEle_eta;   //!
   TBranch        *b_tag_gsfEle_pt;   //!
   TBranch        *b_tag_gsfEle_q;   //!
   TBranch        *b_tag_sc_abseta;   //!
   TBranch        *b_tag_sc_energy;   //!
   TBranch        *b_tag_sc_et;   //!
   TBranch        *b_tag_sc_eta;   //!
   TBranch        *b_tag_dRTau;   //!
   TBranch        *b_tag_isWPMedium;   //!
   TBranch        *b_tag_passingGsf;   //!
   TBranch        *b_mc_probe_abseta;   //!
   TBranch        *b_mc_probe_e;   //!
   TBranch        *b_mc_probe_et;   //!
   TBranch        *b_mc_probe_eta;   //!
   TBranch        *b_mc_probe_mass;   //!
   TBranch        *b_mc_probe_pt;   //!
   TBranch        *b_mc_probe_flag;   //!
   TBranch        *b_pair_abseta;   //!
   TBranch        *b_pair_eta;   //!
   TBranch        *b_pair_mass;   //!
   TBranch        *b_pair_pt;   //!
   TBranch        *b_pair_mass60to120;   //!

   PlotEIDvar(TTree *tree=0);
   virtual ~PlotEIDvar();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PlotEIDvar_cxx
PlotEIDvar::PlotEIDvar(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TnPTree_MC.root");
      if (!f) {
         f = new TFile("TnPTree_MC.root");
         f->cd("TnPTree_MC.root:/GsfElectronToIdMedium");
      }
      tree = (TTree*)gDirectory->Get("fitter_tree");

   }
   Init(tree);
}

PlotEIDvar::~PlotEIDvar()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PlotEIDvar::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PlotEIDvar::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PlotEIDvar::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("probe_gsfEle_abseta", &probe_gsfEle_abseta, &b_probe_gsfEle_abseta);
   fChain->SetBranchAddress("probe_gsfEle_dEtaIn", &probe_gsfEle_dEtaIn, &b_probe_gsfEle_dEtaIn);
   fChain->SetBranchAddress("probe_gsfEle_dPhiIn", &probe_gsfEle_dPhiIn, &b_probe_gsfEle_dPhiIn);
   fChain->SetBranchAddress("probe_gsfEle_e", &probe_gsfEle_e, &b_probe_gsfEle_e);
   fChain->SetBranchAddress("probe_gsfEle_et", &probe_gsfEle_et, &b_probe_gsfEle_et);
   fChain->SetBranchAddress("probe_gsfEle_eta", &probe_gsfEle_eta, &b_probe_gsfEle_eta);
   fChain->SetBranchAddress("probe_gsfEle_hoe", &probe_gsfEle_hoe, &b_probe_gsfEle_hoe);
   fChain->SetBranchAddress("probe_gsfEle_mHits", &probe_gsfEle_mHits, &b_probe_gsfEle_mHits);
   fChain->SetBranchAddress("probe_gsfEle_ooemoop", &probe_gsfEle_ooemoop, &b_probe_gsfEle_ooemoop);
   fChain->SetBranchAddress("probe_gsfEle_pt", &probe_gsfEle_pt, &b_probe_gsfEle_pt);
   fChain->SetBranchAddress("probe_gsfEle_q", &probe_gsfEle_q, &b_probe_gsfEle_q);
   fChain->SetBranchAddress("probe_gsfEle_reltrackiso", &probe_gsfEle_reltrackiso, &b_probe_gsfEle_reltrackiso);
   fChain->SetBranchAddress("probe_gsfEle_sigmaIEtaIEta", &probe_gsfEle_sigmaIEtaIEta, &b_probe_gsfEle_sigmaIEtaIEta);
   fChain->SetBranchAddress("probe_gsfEle_trackiso", &probe_gsfEle_trackiso, &b_probe_gsfEle_trackiso);
   fChain->SetBranchAddress("probe_sc_abseta", &probe_sc_abseta, &b_probe_sc_abseta);
   fChain->SetBranchAddress("probe_sc_energy", &probe_sc_energy, &b_probe_sc_energy);
   fChain->SetBranchAddress("probe_sc_et", &probe_sc_et, &b_probe_sc_et);
   fChain->SetBranchAddress("probe_sc_eta", &probe_sc_eta, &b_probe_sc_eta);
   fChain->SetBranchAddress("probe_d0vtx", &probe_d0vtx, &b_probe_d0vtx);
   fChain->SetBranchAddress("probe_dRTau", &probe_dRTau, &b_probe_dRTau);
   fChain->SetBranchAddress("probe_dzvtx", &probe_dzvtx, &b_probe_dzvtx);
   fChain->SetBranchAddress("probe_iso", &probe_iso, &b_probe_iso);
   fChain->SetBranchAddress("probe_iso_ch", &probe_iso_ch, &b_probe_iso_ch);
   fChain->SetBranchAddress("probe_iso_em", &probe_iso_em, &b_probe_iso_em);
   fChain->SetBranchAddress("probe_iso_nh", &probe_iso_nh, &b_probe_iso_nh);
   fChain->SetBranchAddress("probe_vtxFitConversion", &probe_vtxFitConversion, &b_probe_vtxFitConversion);
   fChain->SetBranchAddress("probe_isWPMedium", &probe_isWPMedium, &b_probe_isWPMedium);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("event_nPV", &event_nPV, &b_mNPV);
   fChain->SetBranchAddress("event_met_calomet", &event_met_calomet, &b_mMET);
   fChain->SetBranchAddress("event_met_calosumet", &event_met_calosumet, &b_mSumET);
   fChain->SetBranchAddress("event_met_calometsignificance", &event_met_calometsignificance, &b_mMETSign);
   fChain->SetBranchAddress("event_met_tcmet", &event_met_tcmet, &b_mtcMET);
   fChain->SetBranchAddress("event_met_tcsumet", &event_met_tcsumet, &b_mtcSumET);
   fChain->SetBranchAddress("event_met_tcmetsignificance", &event_met_tcmetsignificance, &b_mtcMETSign);
   fChain->SetBranchAddress("event_met_pfmet", &event_met_pfmet, &b_mpfMET);
   fChain->SetBranchAddress("event_met_pfsumet", &event_met_pfsumet, &b_mpfSumET);
   fChain->SetBranchAddress("event_met_pfmetsignificance", &event_met_pfmetsignificance, &b_mpfMETSign);
   fChain->SetBranchAddress("event_PrimaryVertex_x", &event_PrimaryVertex_x, &b_mPVx);
   fChain->SetBranchAddress("event_PrimaryVertex_y", &event_PrimaryVertex_y, &b_mPVy);
   fChain->SetBranchAddress("event_PrimaryVertex_z", &event_PrimaryVertex_z, &b_mPVz);
   fChain->SetBranchAddress("event_BeamSpot_x", &event_BeamSpot_x, &b_mBSx);
   fChain->SetBranchAddress("event_BeamSpot_y", &event_BeamSpot_y, &b_mBSy);
   fChain->SetBranchAddress("event_BeamSpot_z", &event_BeamSpot_z, &b_mBSz);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("mcTrue", &mcTrue, &b_mcTrue);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("tag_gsfEle_abseta", &tag_gsfEle_abseta, &b_tag_gsfEle_abseta);
   fChain->SetBranchAddress("tag_gsfEle_e", &tag_gsfEle_e, &b_tag_gsfEle_e);
   fChain->SetBranchAddress("tag_gsfEle_et", &tag_gsfEle_et, &b_tag_gsfEle_et);
   fChain->SetBranchAddress("tag_gsfEle_eta", &tag_gsfEle_eta, &b_tag_gsfEle_eta);
   fChain->SetBranchAddress("tag_gsfEle_pt", &tag_gsfEle_pt, &b_tag_gsfEle_pt);
   fChain->SetBranchAddress("tag_gsfEle_q", &tag_gsfEle_q, &b_tag_gsfEle_q);
   fChain->SetBranchAddress("tag_sc_abseta", &tag_sc_abseta, &b_tag_sc_abseta);
   fChain->SetBranchAddress("tag_sc_energy", &tag_sc_energy, &b_tag_sc_energy);
   fChain->SetBranchAddress("tag_sc_et", &tag_sc_et, &b_tag_sc_et);
   fChain->SetBranchAddress("tag_sc_eta", &tag_sc_eta, &b_tag_sc_eta);
   fChain->SetBranchAddress("tag_dRTau", &tag_dRTau, &b_tag_dRTau);
   fChain->SetBranchAddress("tag_isWPMedium", &tag_isWPMedium, &b_tag_isWPMedium);
   fChain->SetBranchAddress("tag_passingGsf", &tag_passingGsf, &b_tag_passingGsf);
   fChain->SetBranchAddress("mc_probe_abseta", &mc_probe_abseta, &b_mc_probe_abseta);
   fChain->SetBranchAddress("mc_probe_e", &mc_probe_e, &b_mc_probe_e);
   fChain->SetBranchAddress("mc_probe_et", &mc_probe_et, &b_mc_probe_et);
   fChain->SetBranchAddress("mc_probe_eta", &mc_probe_eta, &b_mc_probe_eta);
   fChain->SetBranchAddress("mc_probe_mass", &mc_probe_mass, &b_mc_probe_mass);
   fChain->SetBranchAddress("mc_probe_pt", &mc_probe_pt, &b_mc_probe_pt);
   fChain->SetBranchAddress("mc_probe_flag", &mc_probe_flag, &b_mc_probe_flag);
   fChain->SetBranchAddress("pair_abseta", &pair_abseta, &b_pair_abseta);
   fChain->SetBranchAddress("pair_eta", &pair_eta, &b_pair_eta);
   fChain->SetBranchAddress("pair_mass", &pair_mass, &b_pair_mass);
   fChain->SetBranchAddress("pair_pt", &pair_pt, &b_pair_pt);
   fChain->SetBranchAddress("pair_mass60to120", &pair_mass60to120, &b_pair_mass60to120);
   Notify();
}

Bool_t PlotEIDvar::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PlotEIDvar::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PlotEIDvar::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PlotEIDvar_cxx
