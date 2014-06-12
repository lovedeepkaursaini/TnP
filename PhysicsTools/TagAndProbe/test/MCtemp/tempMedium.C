#define tempMedium_cxx
#include "tempMedium.h"
#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <TFile.h>
#include <fstream>
#include <iomanip>
using namespace std;

double pts[2] = {50,200};
double etas[2] = {0.8,1.4442};

void tempMedium::Loop()
{
   if (fChain == 0) return;
//   ofstream file;
std::ofstream fs("/scratch/kaur/rejectTau.txt",ios::app);
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
//      if(!(event==1296789))continue;
      // if (Cut(ientry) < 0) continue;
      if(!(jentry%100000))std::cout<<jentry<<std::endl;
if( mcTrue && mass > 60 && mass <120 && probe_sc_abseta < 2.0 && probe_sc_abseta >= 1.566 && probe_gsfEle_pt >=15. && probe_gsfEle_pt <20. && tag_dRTau >0.2 && probe_dRTau > 0.2)
//if(event==1296789)
{
//fs<<setw(14)<<mcTrue<<endl;
if(!probe_isWPMedium)
//fs<<lumi<<right<<setw(14)<<event<<endl;
      fs<<run<<right<<setw(14)<<event<<right<<setw(14)<<PUweight<<right<<setw(14)<<mass<<right<<setw(14)<<tag_gsfEle_pt<<right<<setw(14)<<tag_sc_abseta<<right<<setw(14)<<probe_gsfEle_pt<<right<<setw(14)<<probe_sc_abseta<<right<<setw(14)<<probe_isWPMedium<<std::endl;
//      cout<<run<<'\t'<<event<<'\t'<<PUweight<<'\t'<<mass<<'\t'<<tag_gsfEle_pt<<'\t'<<tag_sc_abseta<<'\t'<<probe_gsfEle_pt<<'\t'<<probe_sc_abseta<<'\t'<<probe_isWPMedium<<std::endl;
}
}
    fs.close();
cout<<"writing done"<<endl;
}
