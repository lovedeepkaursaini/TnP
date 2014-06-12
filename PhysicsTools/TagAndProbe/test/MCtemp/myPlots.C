#define myPlots_cxx
#include "myPlots.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void myPlots::Loop()
{
    TFile* file  = new TFile("Data_nPVMedium.root","RECREATE");
    TH1F *MC_nPV=new TH1F("MC_nPV","MC_nPV",20,0,40);
    TH1F *MC_nPV_PU=new TH1F("MC_nPV_PU","MC_nPV_PU",20,0,40);
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      MC_nPV->Fill(event_nPV);
      MC_nPV_PU->Fill(event_nPV,PUweight);
   }
   //   MC_nPV->Draw();
   //   MC_nPV_PU->Draw();
   file->cd();
   file->Write();
}
