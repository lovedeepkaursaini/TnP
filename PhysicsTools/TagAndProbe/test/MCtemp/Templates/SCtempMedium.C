#define SCtempMedium_cxx
#include "SCtempMedium.h"
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

double pts[6] = {10,20,30,40,50,200};
double etas[6] = {0.0,0.8,1.4442,1.556,2.0,2.5};
//double etas[4] = {0.0,1.4442,1.556,2.5};

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

void SCtempMedium::Loop()
{

   if (fChain == 0) return;
   
    TFile* file  = new TFile("SCMedium.root","RECREATE");
   std::map<std::string, TH1D*> histos;
   for(int i=0; i!=5; i++){
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
      double probePt = probe_et;
      double probeEta = probe_abseta;
      std::string ptadd = findBins(6, pts, probePt);
      std::string etaadd = findBins(6, etas, probeEta);
      std::string key = "hMass_"+ptadd+"_"+etaadd;
      if(probe_passingGsf) key = key+"_Pass";
      else key = key+"_Fail";
      if(key.find("dump")!=std::string::npos)key = "dump";
      //std::cout<<key<<std::endl;
//      if(event_met_pfmet<15.)
         (histos[key])->Fill(mass,PUweight);
      //(histos[key])->Fill(mass);
 
}
file->cd();
file->Write();

}
