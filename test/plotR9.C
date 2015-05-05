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
#include "TPaveStats.h"

double calibratedE(const double Etot, const double eta,
		   const unsigned cells=3){
  //calibration for signal region 2: 3*3 cm^2
  double pars[3] = {69.5,4.5,-0.8};
  //double pars[3] = {75.5,0,0};
  double paro[3] = {-34.4,0,0};
  //double pars[3] = {77,3.4,-0.50};
  //double paro[3] = {-11.6,-7.7,-8.8};
  //double paro[3] = {-5.3,-12.8,-6.9};  
  double offset = paro[0] + paro[1]*fabs(eta) + paro[2]*eta*eta;
  double slope = pars[0] + pars[1]*fabs(eta) + pars[2]*eta*eta;
  if (cells==5) slope *= 1.11;
  if (cells==7) slope *= 1.15;
  return (Etot-offset)/slope;
};

int plotR9() {

  std::string cut1pca = "eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetidPCA1==0 && invalidNeighbourPCA1==0 && noPhiCrackPCA1==1 && dRTrue1 < 0.05";
  std::string cut2pca = "eTrue2/cosh(etaTrue2)>40 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetidPCA2==0 && invalidNeighbourPCA2==0 && noPhiCrackPCA2==1 && dRTrue2 < 0.05";

  TFile *fin = TFile::Open("Truth_Hgg_0pu.root");
  if (!fin) return 1;
  fin->cd("hgg");

  TTree *t = (TTree*)gDirectory->Get("tree");

  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }

  const unsigned nC = 4;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "myc" << ic ;
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1);
  }


  myc[0]->cd();
  gStyle->SetOptStat("eMRuo");
  
  TH1F *unconv = new TH1F("unconv",";#frac{e_{3#times3}}{e_{SC}};Photons",100,0,2);
  unconv->Sumw2();
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eSC1>>unconv",(cut1pca+" && converted1==0").c_str());
  std::cout << " -- first photon " << unconv->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eSC2>>+unconv",(cut2pca+" && converted2==0").c_str());
  std::cout << " -- second photon " << unconv->GetEntries() << std::endl;

  TH1F *conv = new TH1F("conv",";#frac{e_{3#times3}}{e_{SC}};Photons",100,0,2);
  conv->Sumw2();
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eSC1>>conv",(cut1pca+" && converted1==1").c_str());
  std::cout << " -- first photon " << conv->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eSC2>>+conv",(cut2pca+" && converted2==1").c_str());
  std::cout << " -- second photon " << conv->GetEntries() << std::endl;
  unconv->Draw();
  gPad->Update();
  TPaveStats *stu = (TPaveStats*)unconv->FindObject("stats");
  if (stu){
    stu->SetLineColor(1);
    stu->SetTextColor(1);
    stu->SetX1NDC(0.63);
    stu->SetX2NDC(0.93);
    stu->SetY1NDC(0.72);
    stu->SetY2NDC(0.92);
    stu->SetLabel("Unconverted");
  }
  conv->SetLineColor(2);
  conv->Draw("sames");
  gPad->Update();
  TPaveStats *st = (TPaveStats*)conv->FindObject("stats");
  if (st){
    st->SetLineColor(2);
    st->SetTextColor(2);
    st->SetX1NDC(0.63);
    st->SetX2NDC(0.93);
    st->SetY1NDC(0.5);
    st->SetY2NDC(0.7);
    st->SetLabel("Converted");
  }

  myc[0]->Update();
  myc[0]->Print("e3x3overeSC_Hgg_140pu.pdf");
    
  TH1F *unconv_int = new TH1F("unconv_int",";#int{#frac{e_{3#times3}}{e_{SC}}};Photons",100,0,2);
  TH1F *conv_int = new TH1F("conv_int",";#int{#frac{e_{3#times3}}{e_{SC}}};Photons",100,0,2);
  const unsigned nbins = conv->GetNbinsX()+1;
  for (unsigned ib(1); ib<nbins;++ib){
    double err;
    unconv_int->SetBinContent(ib,unconv->IntegralAndError(0,ib,err));
    unconv_int->SetBinError(ib,err);
    conv_int->SetBinContent(ib,conv->IntegralAndError(0,ib,err));
    conv_int->SetBinError(ib,err);
 }


  gStyle->SetOptStat(0);
  myc[1]->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  unconv_int->GetXaxis()->SetRangeUser(0.7,1.1);
  unconv_int->Scale(1./unconv->Integral());
  unconv_int->SetMarkerStyle(23);
  unconv_int->Draw("PE");
  conv_int->SetLineColor(2);
  conv_int->SetMarkerColor(2);
  conv_int->SetMarkerStyle(22);
  conv_int->Scale(1./conv->Integral());
  conv_int->Draw("PEsame");


  TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
  leg->SetFillColor(10);
  leg->AddEntry(conv_int,"Converted photons","P");
  leg->AddEntry(unconv_int,"Unconverted photons","P");
  leg->Draw("same");

  myc[1]->Update();
  myc[1]->Print("e3x3overeSCintegrated_Hgg_140pu.pdf");


  myc[2]->cd();
  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);
  
  TH1F *htot = new TH1F("htot",";e_{reco}/e_{true};Photons",100,0,2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1>>htot",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1>0.9").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2>>+htot",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2>0.9").c_str());
  htot->Draw();
  htot->Fit("gaus","","same",0.95,1.05);

  std::cout << " -- first photon " << htot->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR5Reco1,etaPCA1,5)/eTrue1>>+htot",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1<=0.9").c_str());
  t->Draw("calibratedE(eSR5Reco2,etaPCA2,5)/eTrue2>>+htot",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2<=0.9").c_str());
  
  std::cout << " -- second photon " << htot->GetEntries() << std::endl;
  htot->Draw("same");

  htot->Fit("gaus","","same",0.95,1.05);

  myc[2]->Update();
  myc[2]->Print("eReco35overeTrue_Hgg_140pu.pdf");
    
  myc[3]->cd();
  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);
  
  TH1F *htot7 = new TH1F("htot7",";e_{reco}/e_{true};Photons",100,0,2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1>>htot7",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1>0.9").c_str());
  t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/eTrue1>>+htot7",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1<=0.9").c_str());

  std::cout << " -- first photon " << htot->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2>>+htot7",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2>0.9").c_str());
  t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/eTrue2>>+htot7",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2<=0.9").c_str());
  
  std::cout << " -- second photon " << htot->GetEntries() << std::endl;
  htot7->Draw();

  htot7->Fit("gaus","","same",0.95,1.05);

  myc[3]->Update();
  myc[3]->Print("eReco37overeTrue_Hgg_140pu.pdf");
    






  return 0;
}
