//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 12 16:37:56 2012 by ROOT version 5.28/00c
// from TTree fitter_tree/fitter_tree
// found on file: setTightTag-FE34DF89-152A-E111-BE45-0015178C4D94.root
//////////////////////////////////////////////////////////

#ifndef SCtempMedium_h
#define SCtempMedium_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class SCtempMedium {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         probe_abseta;
   Float_t         probe_et;
   Float_t         probe_pt;
   Int_t           probe_passingGsf;
   Float_t         PUweight;
   Float_t         event_met_pfmet;

   Float_t         mass;
   // List of branches
   TBranch        *b_probe_abseta;   //!
   TBranch        *b_probe_et;   //!
   TBranch        *b_probe_pt;   //!

   TBranch        *b_probe_passingGsf;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_mpfMET;   //!
   TBranch        *b_mass;   //!

   SCtempMedium(TTree *tree=0);
   virtual ~SCtempMedium();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SCtempMedium_cxx
SCtempMedium::SCtempMedium(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/hdfs/store/user/kaur/mrgdTnPTreeElec_DYJetsToLLmadgraphSummer12_DR53X_SConly.root");
      if (!f) {
         f = new TFile("/hdfs/store/user/kaur/mrgdTnPTreeElec_DYJetsToLLmadgraphSummer12_DR53X_SConly.root");
         f->cd("/hdfs/store/user/kaur/mrgdTnPTreeElec_DYJetsToLLmadgraphSummer12_DR53X_SConly.root:/SuperClusterToGsfElectron");
      }
      tree = (TTree*)gDirectory->Get("fitter_tree");

   }
   Init(tree);
}

SCtempMedium::~SCtempMedium()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SCtempMedium::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SCtempMedium::LoadTree(Long64_t entry)
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

void SCtempMedium::Init(TTree *tree)
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

   fChain->SetBranchAddress("probe_abseta", &probe_abseta, &b_probe_abseta);
   fChain->SetBranchAddress("probe_pt", &probe_pt, &b_probe_pt);
   fChain->SetBranchAddress("probe_et", &probe_et, &b_probe_et);
   fChain->SetBranchAddress("probe_passingGsf", &probe_passingGsf, &b_probe_passingGsf);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("event_met_pfmet", &event_met_pfmet, &b_mpfMET);
   fChain->SetBranchAddress("mass", &mass, &b_mass);

   Notify();
}

Bool_t SCtempMedium::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SCtempMedium::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SCtempMedium::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SCtempMedium_cxx
