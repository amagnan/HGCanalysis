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

int plotYvsX(){
  TFile *fin = TFile::Open("Truth_Hgg_0pu_maxFromTruth.root");
  if (!fin) return 1;
  fin->cd("hgg");
  const unsigned nLayers = 30;

  TH2F *hyvsx[nLayers];
  TH2F *hxvsz = (TH2F*)gDirectory->Get("hxvsz");
  TH2F *hyvsz = (TH2F*)gDirectory->Get("hyvsz");
  TH2F *hrvsz = (TH2F*)gDirectory->Get("hrvsz");


  for (unsigned iL(0);iL<nLayers;++iL){
    std::ostringstream label;
    label.str("");
    label << "hyvsx_" << iL;
    hyvsx[iL] = (TH2F*)gDirectory->Get(label.str().c_str());
  }

  gStyle->SetOptStat(0);

  const unsigned nC = 8;
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
    gPad->SetLogz(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    hyvsx[iL]->Draw("colz");
    //hyvsx[iL]->Draw("text");
    sprintf(buf,"Layer %d",iL);
    lat.DrawLatex(0,0,buf);
    myc[iL/6]->Update();
  }

  myc[5]->cd();
  hxvsz->Draw("colz");
  myc[5]->Update();

  myc[6]->cd();
  hyvsz->Draw("colz");
  myc[6]->Update();

  myc[7]->cd();
  hrvsz->Draw("colz");
  myc[7]->Update();


  return 0;
}
