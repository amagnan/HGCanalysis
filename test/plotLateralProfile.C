#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TProfile.h"

int plotLateralProfile(){
  TFile *fin = TFile::Open("Truth_Hgg_0pu.root");
  if (!fin) return 1;
  fin->cd("hgg");
  const unsigned nLayers = 30;
  
  std::string cut1 = "converted1==0 && eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetid1==0 && invalidNeighbour1==0";
  std::string cut2 = "converted2==0 && eTrue2/cosh(etaTrue2)>40 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetid2==0 && invalidNeighbour2==0";

  unsigned numSR=3;
  unsigned denSR=5;

  TTree *t = (TTree*)gDirectory->Get("tree");

  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }

  gStyle->SetOptStat("euo");

  const unsigned nC = 5;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "myc" << ic ;
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1500,1000);
    if (ic<5) myc[ic]->Divide(3,2);
  }

  TLatex lat;
  char buf[100];
  for (unsigned iL(0);iL<nLayers;++iL){
    myc[iL/6]->cd(iL%6+1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    std::ostringstream label;
    label.str("");
    label << "eSR" << numSR << "RecoLayer"<<iL<<"_1"<<"/eSR" << denSR << "RecoLayer"<<iL<<"_1";
    t->Draw(label.str().c_str(),cut1.c_str());
    std::ostringstream histname;
    histname << "hist"<<iL;
    TH1F *hist = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(histname.str().c_str());
    label.str("");
    label << "eSR" << numSR << "RecoLayer"<<iL<<"_2"<<"/eSR" << denSR << "RecoLayer"<<iL<<"_2>>"<<histname.str().c_str();
    t->Draw(label.str().c_str(),cut2.c_str());
    label.str("");
    label << ";eSR" << numSR << "/eSR" << denSR << ";photons";
    hist->SetTitle(label.str().c_str());
    hist->Draw();
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatex(0.8,0.9*hist->GetMaximum(),buf);
    myc[iL/6]->Update();

  }

  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "PLOTS/LatProf_Hgg_0pu_" << numSR << "over" << denSR << "_l"<<6*ic << "-" << 6*(ic+1)-1 << ".pdf";
    myc[ic]->Print(label.str().c_str());
  }

  return 0;
}
