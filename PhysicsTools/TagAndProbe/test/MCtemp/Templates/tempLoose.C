#define tempLoose_cxx
#include "tempLoose.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <TFile.h>
#include <TH1D.h>
using std::string;

double pts[7] = {10,15,20,30,40,50,200};
double etas[6] = {0.0,0.8,1.4442,1.566,2.0,2.5};

std::string findBins(double size, double* array, double pT){
  std::string bin = "dump";
  for(int i=0; i!=size-1;i++){
    double low = array[i];
    double hi  = array[i+1];
    if((low <= pT)&&(hi>pT)){
	std::stringstream ss;
	ss<<low<<"To"<<hi;
	bin = ss.str();
      }
  }
  return bin;
}

void tempLoose::Loop()
{

   if (fChain == 0) return;
   
    TFile* file  = new TFile("Loose.root","RECREATE");
   std::map<std::string, TH1D*> histos;
   for(int i=0; i!=6; i++){
     for(int j=0; j!=5;j++){
       std::stringstream histNameSt;
       histNameSt<<"hMass_"<<pts[i]<<"To"<<pts[i+1]<<"_"
		 <<etas[j]<<"To"<<etas[j+1];
       std::string hp = histNameSt.str()+std::string("_Pass");
       std::string hf = histNameSt.str()+std::string("_Fail");
       TString p(hp.c_str());
       TString f(hf.c_str());
       histos[hp] = new TH1D(p,p,30,60,120);
       histos[hf] = new TH1D(f,f,30,60,120);
     }
   }
   histos["dump"]=new TH1D("dump","dump",1,0,1000);
  
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(!(jentry%100000))std::cout<<jentry<<std::endl;
      double probePt = probe_gsfEle_pt;
      double probeEta = probe_sc_abseta;
      std::string ptadd = findBins(7, pts, probePt);
      std::string etaadd = findBins(6, etas, probeEta);
      std::string key = "hMass_"+ptadd+"_"+etaadd;
      if(probe_isWPLoose) key = key+"_Pass";
      else key = key+"_Fail";
      if(key.find("dump")!=std::string::npos)key = "dump";
      //std::cout<<key<<std::endl;
      if( mcTrue && mass > 60 && mass <120 && tag_dRTau >0.2 && probe_dRTau > 0.2) (histos[key])->Fill(mass,PUweight);
//      (histos[key])->Fill(mass,PUweight);
      //(histos[key])->Fill(mass);
 
}
file->cd();
file->Write();

}
