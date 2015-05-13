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
#include "effSigmaMacro.C"
#include "tdrstyle.C"

double DeltaPhi(const double & phi1, const double & phi2){
  double dphi = phi1 - phi2;
  if (dphi< (-1.*TMath::Pi())) dphi += 2*TMath::Pi();
  if (dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
  return dphi;
}

double DeltaR(const double & eta1, const double & eta2,
		const double & phi1, const double & phi2){
  double dphi = DeltaPhi(phi1,phi2);
  double deta = eta1-eta2;
  double dr = TMath::Sqrt(dphi*dphi+deta*deta);
  return dr;
}

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
  if (cells==1) {
    offset += slope*-4.3;//in original unit=mips...
    slope *= 0.55;
  }
  if (cells==5) slope *= 1.11;
  if (cells==7) slope *= 1.15;
  return (Etot-offset)/slope;
};

int plotR9() {

  bool doPu = true;
  setTDRStyle();

  std::string cut1pca = "eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetidPCA1==0 && invalidNeighbourPCA1==0 && dRTrue1 < 0.05 && noPhiCrackPCA1==1";
  std::string cut2pca = "eTrue2/cosh(etaTrue2)>40 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetidPCA2==0 && invalidNeighbourPCA2==0 && dRTrue2 < 0.05 && noPhiCrackPCA2==1";

  TFile *fin;
  if (doPu) fin = TFile::Open("Truth_Hgg_140pu.root");
  else fin = TFile::Open("Truth_Hgg_0pu.root");
  if (!fin) return 1;
  fin->cd("hgg");

  TTree *t = (TTree*)gDirectory->Get("tree");

  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }

  const unsigned nC = 16;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "myc" << ic ;
    //if (ic<9 || ic>10) myc[ic] = new TCanvas(label.str().c_str(),
    // label.str().c_str(),
    //1);
    //else 
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1000,600);
    //myc[ic]->SetRightMargin(0.05);
    //myc[ic]->SetTopMargin(0.05);
  }

  TLatex lat;
  TLatex latSmall;
  latSmall.SetTextSize(0.035);
  char buf[100];

  myc[0]->cd();
  //gStyle->SetOptStat("eMRuo");
  
  TH1F *unconv = new TH1F("unconv",";e_{3#times3}/e_{SC};Photons",50,0,1.5);
  unconv->Sumw2();
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eSC1>>unconv",(cut1pca+" && converted1==0").c_str());
  std::cout << " -- first photon " << unconv->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eSC2>>+unconv",(cut2pca+" && converted2==0").c_str());
  std::cout << " -- second photon " << unconv->GetEntries() << std::endl;

  TH1F *conv = new TH1F("conv",";e_{3#times3}/e_{SC};Photons",50,0,1.5);
  conv->Sumw2();
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eSC1>>conv",(cut1pca+" && converted1==1").c_str());
  std::cout << " -- first photon " << conv->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eSC2>>+conv",(cut2pca+" && converted2==1").c_str());
  std::cout << " -- second photon " << conv->GetEntries() << std::endl;

  TH1F *tot = (TH1F*)unconv->Clone("tot");
  tot->Add(conv);

  tot->SetLineColor(4);
  tot->SetFillColor(4);
  tot->SetFillStyle(3004);
  //tot->SetMarkerColor(4);
  //tot->SetMarkerStyle(20);
  tot->Draw("hist");
  unconv->SetMarkerStyle(23);
  unconv->Draw("same");
  gPad->Update();
  conv->SetLineColor(2);
  conv->SetMarkerColor(2);
  conv->SetMarkerStyle(22);
  conv->Draw("same");

  TLegend *leg = new TLegend(0.18,0.65,0.62,0.9);
  leg->SetFillStyle(0);
  leg->AddEntry(tot,"H#rightarrow#gamma#gamma p_{T}^{#gamma}>40 GeV","F");
  leg->AddEntry(conv,"Converted photons","P");
  leg->AddEntry(unconv,"Unconverted photons","P");
  leg->Draw("same");

  /*
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
    }*/

  if (doPu) latSmall.DrawLatex(0.2,1.07*tot->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU140");
  else latSmall.DrawLatex(0.2,1.07*tot->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU0");
  myc[0]->Update();
  if (doPu) myc[0]->Print("PLOTS_R9/e3x3overeSC_Hgg_140pu.pdf");
  else myc[0]->Print("PLOTS_R9/e3x3overeSC_Hgg_0pu.pdf");
    
  TH1F *unconv_int = new TH1F("unconv_int",";e_{3#times3}/e_{SC};Photon fraction e_{3#times3}/e_{SC}<val",50,0,1.5);
  TH1F *conv_int = new TH1F("conv_int",";e_{3#times3}/e_{SC};Photon fraction e_{3#times3}/e_{SC}<val",50,0,1.5);
  TH1F *tot_int = new TH1F("tot_int",";e_{3#times3}/e_{SC};Photon fraction e_{3#times3}/e_{SC}<val",50,0,1.5);
  const unsigned nbins = conv->GetNbinsX()+1;
  for (unsigned ib(1); ib<nbins;++ib){
    double err;
    unconv_int->SetBinContent(ib,unconv->IntegralAndError(0,ib,err));
    unconv_int->SetBinError(ib,err);
    conv_int->SetBinContent(ib,conv->IntegralAndError(0,ib,err));
    conv_int->SetBinError(ib,err);
    tot_int->SetBinContent(ib,tot->IntegralAndError(0,ib,err));
    tot_int->SetBinError(ib,err);
 }


  gStyle->SetOptStat(0);
  myc[1]->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  unconv_int->GetXaxis()->SetRangeUser(0.4,1.2);
  unconv_int->Scale(1./unconv->Integral());
  unconv_int->SetMarkerStyle(23);
  unconv_int->Draw("PE");
  conv_int->SetLineColor(2);
  conv_int->SetMarkerColor(2);
  conv_int->SetMarkerStyle(22);
  conv_int->Scale(1./conv->Integral());
  conv_int->Draw("PEsame");
  tot_int->SetLineColor(4);
  tot_int->SetMarkerColor(4);
  tot_int->SetMarkerStyle(20);
  tot_int->Scale(1./tot->Integral());
  tot_int->Draw("PEsame");

  leg->Draw("same");

  if (doPu) latSmall.DrawLatex(0.5,1.1*tot_int->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU140");
  else latSmall.DrawLatex(0.5,1.1*tot_int->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU0");
  myc[1]->Update();
  if (doPu) myc[1]->Print("PLOTS_R9/e3x3overeSCintegrated_Hgg_140pu.pdf");
  else myc[1]->Print("PLOTS_R9/e3x3overeSCintegrated_Hgg_0pu.pdf");


  myc[4]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  
  TH1F *htot3 = new TH1F("htot3",";e_{reco}/e_{true};Photons",100,0,1.2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1>>htot3",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1>0.9").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2>>+htot3",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2>0.9").c_str());
  TH1F *h333 = (TH1F*)htot3->Clone("h333");
  h333->SetFillColor(6);
  h333->SetFillStyle(3004);
  h333->SetLineColor(6);
  //h33->Fit("gaus","","same",0.95,1.05);

  std::cout << " -- unconv photon " << htot3->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1>>+htot3",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1<=0.9").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2>>+htot3",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2<=0.9").c_str());
  
  std::cout << " -- all photons " << htot3->GetEntries() << std::endl;
  htot3->Draw();

  //htot3->Fit("gaus","","same",0.95,1.05);
  h333->Draw("same");

  lat.DrawLatex(0.1,0.9*htot3->GetMaximum(),"#frac{e_{3#times3}}{e_{SC}}>0.9: e_{reco}=e_{3#times3}");
  lat.DrawLatex(0.1,0.75*htot3->GetMaximum(),"#frac{e_{3#times3}}{e_{SC}}<=0.9: e_{reco}=e_{3#times3}");

  sprintf(buf,"N_{#gamma}^{unconv}/N_{tot}= %3.2f",h333->GetEntries()/htot3->GetEntries());
  lat.DrawLatex(0.1,0.6*htot3->GetMaximum(),buf);

  double err=0;
  double sigeff = effSigmaMacro(h333,err);
  sprintf(buf,"#sigma_{eff}^{unconv}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(6);
  lat.DrawLatex(0.1,0.45*htot3->GetMaximum(),buf);
  sigeff = effSigmaMacro(htot3,err);
  sprintf(buf,"#sigma_{eff}^{all}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(1);
  lat.DrawLatex(0.1,0.3*htot3->GetMaximum(),buf);
 


  myc[4]->Update();
  if (doPu) myc[4]->Print("PLOTS_R9/eReco33overeTrue_Hgg_140pu.pdf");
  else myc[4]->Print("PLOTS_R9/eReco33overeTrue_Hgg_0pu.pdf");

  myc[2]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  
  TH1F *htot = new TH1F("htot",";e_{reco}/e_{true};Photons",100,0,1.2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1>>htot",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1>0.9").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2>>+htot",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2>0.9").c_str());
  TH1F *h33 = (TH1F*)htot->Clone("h33");
  h33->SetFillColor(6);
  h33->SetFillStyle(3004);
  h33->SetLineColor(6);
  //h33->Fit("gaus","","same",0.95,1.05);

  std::cout << " -- unconv photon " << htot->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR5Reco1,etaPCA1,5)/eTrue1>>+htot",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1<=0.9").c_str());
  t->Draw("calibratedE(eSR5Reco2,etaPCA2,5)/eTrue2>>+htot",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2<=0.9").c_str());
  
  std::cout << " -- all photons " << htot->GetEntries() << std::endl;
  htot->Draw();

  //htot->Fit("gaus","","same",0.95,1.05);
  h33->Draw("same");

  lat.DrawLatex(0.1,0.9*htot->GetMaximum(),"#frac{e_{3#times3}}{e_{SC}}>0.9: e_{reco}=e_{3#times3}");
  lat.DrawLatex(0.1,0.75*htot->GetMaximum(),"#frac{e_{3#times3}}{e_{SC}}<=0.9: e_{reco}=e_{5#times5}");

  sprintf(buf,"N_{#gamma}^{unconv}/N_{tot}= %3.2f",h33->GetEntries()/htot->GetEntries());
  lat.DrawLatex(0.1,0.6*htot->GetMaximum(),buf);

  sigeff = effSigmaMacro(h33,err);
  sprintf(buf,"#sigma_{eff}^{unconv}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(6);
  lat.DrawLatex(0.1,0.45*htot->GetMaximum(),buf);
  sigeff = effSigmaMacro(htot,err);
  sprintf(buf,"#sigma_{eff}^{all}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(1);
  lat.DrawLatex(0.1,0.3*htot->GetMaximum(),buf);
 
  myc[2]->Update();
  if (doPu) myc[2]->Print("PLOTS_R9/eReco35overeTrue_Hgg_140pu.pdf");
  else myc[2]->Print("PLOTS_R9/eReco35overeTrue_Hgg_0pu.pdf");
    
  myc[3]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  
  TH1F *htot7 = new TH1F("htot7",";e_{reco}/e_{true};Photons",100,0,1.2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1>>htot7",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1>0.9").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2>>+htot7",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2>0.9").c_str());
  TH1F *h337 = (TH1F*)htot7->Clone("h337");
  h337->SetFillColor(6);
  h337->SetFillStyle(3004);
  h337->SetLineColor(6);
  //h33->Fit("gaus","","same",0.95,1.05);

  std::cout << " -- unconv photon " << htot7->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/eTrue1>>+htot7",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/eSC1<=0.9").c_str());
  t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/eTrue2>>+htot7",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/eSC2<=0.9").c_str());
  
  std::cout << " -- all photons " << htot7->GetEntries() << std::endl;
  htot7->Draw();
  //htot7->Fit("gaus","","same",0.95,1.05);
  h337->Draw("same");

  lat.DrawLatex(0.1,0.9*htot7->GetMaximum(),"#frac{e_{3#times3}}{e_{SC}}>0.9: e_{reco}=e_{3#times3}");
  lat.DrawLatex(0.1,0.75*htot7->GetMaximum(),"#frac{e_{3#times3}}{e_{SC}}<=0.9: e_{reco}=e_{7#times7}");

  sprintf(buf,"N_{#gamma}^{unconv}/N_{tot}= %3.2f",h337->GetEntries()/htot7->GetEntries());
  lat.DrawLatex(0.1,0.6*htot7->GetMaximum(),buf);

  sigeff = effSigmaMacro(h337,err);
  sprintf(buf,"#sigma_{eff}^{unconv}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(6);
  lat.DrawLatex(0.1,0.45*htot7->GetMaximum(),buf);
  sigeff = effSigmaMacro(htot7,err);
  sprintf(buf,"#sigma_{eff}^{all}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(1);
  lat.DrawLatex(0.1,0.3*htot7->GetMaximum(),buf);
 
  myc[3]->Update();
  if (doPu) myc[3]->Print("PLOTS_R9/eReco37overeTrue_Hgg_140pu.pdf");
  else myc[3]->Print("PLOTS_R9/eReco37overeTrue_Hgg_0pu.pdf");
    

  myc[5]->cd();
  TH1F *heta1 = new TH1F("heta1",";|#eta|;Photons",15,1.5,3);
  heta1->Sumw2();
  t->Draw("etaTrue1>>heta1","converted1==1 && eTrue1/cosh(etaTrue1)>40");// && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetid1==0 && invalidNeighbour1==0");
  t->Draw("etaTrue2>>+heta1","converted2==1 && eTrue2/cosh(etaTrue2)>40");// && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetid2==0 && invalidNeighbour2==0");
  TH1F *heta2 = new TH1F("heta2",";|#eta|;Photons",15,1.5,3);
  heta2->Sumw2();
  t->Draw("etaTrue1>>heta2","eTrue1/cosh(etaTrue1)>40");// && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetid1==0 && invalidNeighbour1==0");
  t->Draw("etaTrue2>>+heta2","eTrue2/cosh(etaTrue2)>40");// && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetid2==0 && invalidNeighbour2==0");

  TH1F *heta = new TH1F("heta",";|#eta_{gen}^{#gamma}|;Photons",15,1.5,3);
  heta->Sumw2();
  heta->Divide(heta1,heta2);

  heta->SetMarkerStyle(1);
  gPad->SetGridy(1);
  gPad->SetGridx(1);
  heta->GetYaxis()->SetRangeUser(0.,1.0);
  heta->Draw("PE");
  lat.DrawLatex(1.6,0.9,"H#rightarrow#gamma#gamma p_{T}^{#gamma}>40 GeV");//, #Delta#phi_{crack}>2^{o}");
  //  lat.DrawLatex(-3,1.08*heta->GetMaximum(),"CMSSW_6_2_0_SLHC25_patch2 Relval Hgg PU140");

  myc[5]->Update();
  if (doPu) myc[5]->Print("PLOTS_R9/fractionTrueConvvseta_Hgg_140pu.pdf");
  else myc[5]->Print("PLOTS_R9/fractionTrueConvvseta_Hgg_0pu.pdf");


  myc[6]->cd();
  //gStyle->SetOptStat("eMRuo");
  
  TH1F *unconv2 = new TH1F("unconv2",";e_{3#times3}/e_{5#times5};Photons",50,0,1.1);
  unconv2->Sumw2();
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)>>unconv2",(cut1pca+" && converted1==0").c_str());
  std::cout << " -- first photon " << unconv2->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)>>+unconv2",(cut2pca+" && converted2==0").c_str());
  std::cout << " -- second photon " << unconv2->GetEntries() << std::endl;

  TH1F *conv2 = new TH1F("conv2",";e_{3#times3}/e_{5#times5};Photons",50,0,1.1);
  conv2->Sumw2();
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)>>conv2",(cut1pca+" && converted1==1").c_str());
  std::cout << " -- first photon " << conv2->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)>>+conv2",(cut2pca+" && converted2==1").c_str());
  std::cout << " -- second photon " << conv2->GetEntries() << std::endl;

  TH1F *tot2 = (TH1F*)unconv2->Clone("tot");
  tot2->Add(conv2);

  tot2->SetLineColor(4);
  tot2->SetFillColor(4);
  tot2->SetFillStyle(3004);
  //tot2->SetMarkerColor(4);
  //tot2->SetMarkerStyle(20);
  tot2->Draw("hist");
  unconv2->SetMarkerStyle(23);
  unconv2->Draw("same");
  gPad->Update();
  conv2->SetLineColor(2);
  conv2->SetMarkerColor(2);
  conv2->SetMarkerStyle(22);
  conv2->Draw("same");

  leg->Draw("same");

  if (doPu) latSmall.DrawLatex(0.2,1.07*tot2->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU140");
  else latSmall.DrawLatex(0.2,1.07*tot2->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU0");
  myc[6]->Update();
  if (doPu) myc[6]->Print("PLOTS_R9/e3x3overe5x5_Hgg_140pu.pdf");
  else myc[6]->Print("PLOTS_R9/e3x3overe5x5_Hgg_0pu.pdf");
    
  TH1F *unconv2_int = new TH1F("unconv2_int",";e_{3#times3}/e_{5#times5};Photon fraction e_{3#times3}/e_{5#times5}<val",50,0,1.1);
  TH1F *conv2_int = new TH1F("conv2_int",";e_{3#times3}/e_{5#times5};Photon fraction e_{3#times3}/e_{5#times5}<val",50,0,1.1);
  TH1F *tot2_int = new TH1F("tot2_int",";e_{3#times3}/e_{5#times5};Photon fraction e_{3#times3}/e_{5#times5}<val",50,0,1.1);
  const unsigned nbins2 = conv2->GetNbinsX()+1;
  for (unsigned ib(1); ib<nbins2;++ib){
    double err;
    unconv2_int->SetBinContent(ib,unconv2->IntegralAndError(0,ib,err));
    unconv2_int->SetBinError(ib,err);
    conv2_int->SetBinContent(ib,conv2->IntegralAndError(0,ib,err));
    conv2_int->SetBinError(ib,err);
    tot2_int->SetBinContent(ib,tot2->IntegralAndError(0,ib,err));
    tot2_int->SetBinError(ib,err);
 }


  gStyle->SetOptStat(0);
  myc[7]->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  unconv2_int->GetXaxis()->SetRangeUser(0.4,1.2);
  unconv2_int->Scale(1./unconv2->Integral());
  unconv2_int->SetMarkerStyle(23);
  unconv2_int->Draw("PE");
  conv2_int->SetLineColor(2);
  conv2_int->SetMarkerColor(2);
  conv2_int->SetMarkerStyle(22);
  conv2_int->Scale(1./conv2->Integral());
  conv2_int->Draw("PEsame");
  tot2_int->SetLineColor(4);
  tot2_int->SetMarkerColor(4);
  tot2_int->SetMarkerStyle(20);
  tot2_int->Scale(1./tot2->Integral());
  tot2_int->Draw("PEsame");


  leg->Draw("same");

  if (doPu) latSmall.DrawLatex(0.5,1.1*tot2_int->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU140");
  else latSmall.DrawLatex(0.5,1.1*tot2_int->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU0");
  myc[7]->Update();
  if (doPu) myc[7]->Print("PLOTS_R9/e3x3overee5x5integrated_Hgg_140pu.pdf");
  else myc[7]->Print("PLOTS_R9/e3x3overee5x5integrated_Hgg_0pu.pdf");


  myc[8]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  
  TH1F *htot27 = new TH1F("htot27",";e_{reco}/e_{true};Photons",100,0,1.2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/((1.017-1.29*calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)+8.5*TMath::Power(calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1),2))*eTrue1)>>htot27",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)>0.99").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/((1.017-1.29*calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2)/calibratedE(eSR3Reco2,etaPCA2)+8.5*TMath::Power(calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2)/calibratedE(eSR3Reco2,etaPCA2),2))*eTrue2)>>+htot27",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)>0.99").c_str());
  TH1F *h3327 = (TH1F*)htot27->Clone("h3327");
  h3327->SetFillColor(6);
  h3327->SetFillStyle(3004);
  h3327->SetLineColor(6);
  //h33->Fit("gaus","","same",0.95,1.05);

  std::cout << " -- unconv photon " << htot27->GetEntries() << std::endl;
  //t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/eTrue1>>+htot27",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)<=0.99").c_str());
  //t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/eTrue2>>+htot27",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)<=0.99").c_str());
  if (doPu) t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/((1.022-0.80*calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1))*eTrue1)>>+htot27",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)<=0.99").c_str());
  else t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/((0.993-0.60*calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1))*eTrue1)>>+htot27",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)<=0.99").c_str());
  if (doPu) t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/((1.022-0.80*calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2)/calibratedE(eSR3Reco2,etaPCA2))*eTrue2)>>+htot27",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)<=0.99").c_str());
  else t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/((0.993-0.60*calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2)/calibratedE(eSR3Reco2,etaPCA2))*eTrue2)>>+htot27",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)<=0.99").c_str());

  std::cout << " -- all photons " << htot27->GetEntries() << std::endl;
  htot27->Draw();
  //htot27->Fit("gaus","","same",0.95,1.05);
  h3327->Draw("same");

  lat.DrawLatex(0.1,0.9*htot27->GetMaximum(),"#frac{e_{3#times3}}{e_{5#times5}}>0.99: e_{reco}=e_{3#times3}");
  lat.DrawLatex(0.1,0.75*htot27->GetMaximum(),"#frac{e_{3#times3}}{e_{5#times5}}<=0.99: e_{reco}=e_{7#times7}");

  sprintf(buf,"N_{#gamma}^{unconv}/N_{tot}= %3.2f",h3327->GetEntries()/htot27->GetEntries());
  lat.DrawLatex(0.1,0.6*htot27->GetMaximum(),buf);

  sigeff = effSigmaMacro(h3327,err);
  sprintf(buf,"#sigma_{eff}^{unconv}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(6);
  lat.DrawLatex(0.1,0.45*htot27->GetMaximum(),buf);
  sigeff = effSigmaMacro(htot27,err);
  sprintf(buf,"#sigma_{eff}^{all}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(1);
  lat.DrawLatex(0.1,0.3*htot27->GetMaximum(),buf);
 
  myc[8]->Update();
  if (doPu) myc[8]->Print("PLOTS_R9/eReco37overeTrue_conv5x5_cor1-5_Hgg_140pu.pdf");
  else myc[8]->Print("PLOTS_R9/eReco37overeTrue_conv5x5_cor1-5_Hgg_0pu.pdf");
    

  myc[9]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  gPad->SetRightMargin(0.15);
  TH2F *h2Dtot27 = new TH2F("h2Dtot27",";|#Delta#Phi(Seed,SC)|;e_{reco}/e_{true};Photons",50,0,0.05,50,0,1.2);
  //t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1:TMath::Abs(DeltaPhi(phiSeed1,phiSC1))>>h2Dtot27",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)>0.99").c_str());
  //t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2:TMath::Abs(DeltaPhi(phiSeed2,phiSC2))>>+h2Dtot27",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)>0.99").c_str());
  //TH2F *h2D3327 = (TH2F*)h2Dtot27->Clone("h2D3327");

  //std::cout << " -- unconv photon " << h2Dtot27->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/eTrue1:TMath::Abs(DeltaPhi(phiSeed1,phiSC1))>>+h2Dtot27",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)<=0.99").c_str());
  t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/eTrue2:TMath::Abs(DeltaPhi(phiSeed2,phiSC2))>>+h2Dtot27",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)<=0.99").c_str());
  
  std::cout << " -- all photons " << h2Dtot27->GetEntries() << std::endl;
  gPad->SetLogz(1);
  h2Dtot27->GetZaxis()->SetTitleOffset(0.7);
  h2Dtot27->Draw("colz");
  (h2Dtot27->ProfileX())->Draw("same");
  TProfile *h27_pfx = (TProfile*)gDirectory->Get("h2Dtot27_pfx");
  h27_pfx->SetMarkerStyle(22);
  h27_pfx->SetMarkerColor(1);
  h27_pfx->Draw("PEsame");
  h27_pfx->Fit("pol1","RI","same",0,0.02);
  TF1 *fit = h27_pfx->GetFunction("pol1");
  //h2D3327->SetMarkerStyle(1);
  //h2D3327->SetMarkerColor(6);
  //h2D3327->SetLineColor(6);
  //h2D3327->Draw("Psame");
  //(h2D3327->ProfileX())->Draw("same");

  myc[9]->Update();
  if (doPu) myc[9]->Print("PLOTS_R9/eReco37overeTruevsdphi_conv5x5_Hgg_140pu.pdf");
  else myc[9]->Print("PLOTS_R9/eReco37overeTruevsdphi_conv5x5_Hgg_0pu.pdf");
    

  myc[10]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  gPad->SetRightMargin(0.15);
  TH2F *h2Dtot7_11o33 = new TH2F("h2Dtot7_11o33",";e_{1#times1}/e_{3#times3};e_{reco}/e_{true};Photons",30,0.,1.5,50,0,1.2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1:calibratedE(eSR1Reco1,etaPCA1,1)/calibratedE(eSR3Reco1,etaPCA1)>>h2Dtot7_11o33",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)>0.99").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2:calibratedE(eSR1Reco2,etaPCA2,1)/calibratedE(eSR3Reco2,etaPCA2,3)>>+h2Dtot7_11o33",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)>0.99").c_str());
  TH2F *h2D337_11o33 = (TH2F*)h2Dtot7_11o33->Clone("h2D337_11o33");

  std::cout << " -- unconv photon " << h2Dtot7_11o33->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/eTrue1:calibratedE(eSR1Reco1,etaPCA1,1)/calibratedE(eSR3Reco1,etaPCA1)>>+h2Dtot7_11o33",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)<=0.99").c_str());
  t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/eTrue2:calibratedE(eSR1Reco2,etaPCA2,1)/calibratedE(eSR3Reco2,etaPCA2,3)>>+h2Dtot7_11o33",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)<=0.99").c_str());
  
  std::cout << " -- all photons " << h2Dtot7_11o33->GetEntries() << std::endl;
  gPad->SetLogz(1);
  h2Dtot7_11o33->GetZaxis()->SetTitleOffset(0.7);
  h2Dtot7_11o33->Draw("colz");
  (h2Dtot7_11o33->ProfileX())->Draw("same");
  TProfile *h7_11o33_pfx = (TProfile*)gDirectory->Get("h2Dtot7_11o33_pfx");
  h7_11o33_pfx->SetMarkerStyle(22);
  h7_11o33_pfx->SetMarkerColor(1);
  h7_11o33_pfx->Draw("PEsame");
  h7_11o33_pfx->Fit("pol1","RI","same",0.4,0.9);
  fit = h7_11o33_pfx->GetFunction("pol1");
  fit->Draw("same");
  h2D337_11o33->SetMarkerStyle(1);
  h2D337_11o33->SetMarkerColor(6);
  h2D337_11o33->SetLineColor(6);
  h2D337_11o33->Draw("Psame");
  (h2D337_11o33->ProfileX())->Draw("same");

  myc[10]->Update();
  if (doPu) myc[10]->Print("PLOTS_R9/eReco37overeTruevse11oe33_conv5x5_Hgg_140pu.pdf");
  else myc[10]->Print("PLOTS_R9/eReco37overeTruevse11oe33_conv5x5_Hgg_0pu.pdf");
    



  myc[11]->cd();
  //gStyle->SetOptStat("eMRuo");
  
  TH1F *unconv3 = new TH1F("unconv3",";e_{3#times3}(L1-5)/e_{3#times3};Photons",50,0,0.1);
  unconv3->Sumw2();
  t->Draw("calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)>>unconv3",(cut1pca+" && converted1==0").c_str());
  std::cout << " -- first photon " << unconv3->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2)/calibratedE(eSR3Reco2,etaPCA2)>>+unconv3",(cut2pca+" && converted2==0").c_str());
  std::cout << " -- second photon " << unconv3->GetEntries() << std::endl;

  TH1F *conv3 = new TH1F("conv3",";e_{3#times3}(L1-5)/e_{3#times3};Photons",50,0,0.1);
  conv3->Sumw2();
  t->Draw("calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)>>conv3",(cut1pca+" && converted1==1").c_str());
  std::cout << " -- first photon " << conv3->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2)/calibratedE(eSR3Reco2,etaPCA2)>>+conv3",(cut2pca+" && converted2==1").c_str());
  std::cout << " -- second photon " << conv3->GetEntries() << std::endl;

  TH1F *tot3 = (TH1F*)unconv3->Clone("tot");
  tot3->Add(conv3);

  tot3->SetLineColor(4);
  tot3->SetFillColor(4);
  tot3->SetFillStyle(3004);
  //tot3->SetMarkerColor(4);
  //tot3->SetMarkerStyle(20);
  tot3->Draw("hist");
  unconv3->SetMarkerStyle(23);
  unconv3->Draw("same");
  gPad->Update();
  conv3->SetLineColor(2);
  conv3->SetMarkerColor(2);
  conv3->SetMarkerStyle(22);
  conv3->Draw("same");

  leg->Draw("same");

  if (doPu) latSmall.DrawLatex(0.2,1.07*tot3->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU140");
  else latSmall.DrawLatex(0.2,1.07*tot3->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU0");
  myc[11]->Update();
  if (doPu) myc[11]->Print("PLOTS_R9/e3x3L1-5overe3x3_Hgg_140pu.pdf");
  else myc[11]->Print("PLOTS_R9/e3x3L1-5overe3x3_Hgg_0pu.pdf");
    
  TH1F *unconv3_int = new TH1F("unconv3_int",";e_{3#times3}(L1-5)/e_{3#times3};Photon fraction e_{3#times3}(L1-5)/e_{3#times3}<val",50,0,0.1);
  TH1F *conv3_int = new TH1F("conv3_int",";e_{3#times3}(L1-5)/e_{3#times3};Photon fraction e_{3#times3}(L1-5)/e_{3#times3}<val",50,0,0.1);
  TH1F *tot3_int = new TH1F("tot3_int",";e_{3#times3}(L1-5)/e_{3#times3};Photon fraction e_{3#times3}(L1-5)/e_{3#times3}<val",50,0,0.1);
  const unsigned nbins3 = conv3->GetNbinsX()+1;
  for (unsigned ib(1); ib<nbins3;++ib){
    double err;
    unconv3_int->SetBinContent(ib,unconv3->IntegralAndError(0,ib,err));
    unconv3_int->SetBinError(ib,err);
    conv3_int->SetBinContent(ib,conv3->IntegralAndError(0,ib,err));
    conv3_int->SetBinError(ib,err);
    tot3_int->SetBinContent(ib,tot3->IntegralAndError(0,ib,err));
    tot3_int->SetBinError(ib,err);
 }


  gStyle->SetOptStat(0);
  myc[12]->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  unconv3_int->GetXaxis()->SetRangeUser(0.4,1.2);
  unconv3_int->Scale(1./unconv3->Integral());
  unconv3_int->SetMarkerStyle(23);
  unconv3_int->Draw("PE");
  conv3_int->SetLineColor(2);
  conv3_int->SetMarkerColor(2);
  conv3_int->SetMarkerStyle(22);
  conv3_int->Scale(1./conv3->Integral());
  conv3_int->Draw("PEsame");
  tot3_int->SetLineColor(4);
  tot3_int->SetMarkerColor(4);
  tot3_int->SetMarkerStyle(20);
  tot3_int->Scale(1./tot3->Integral());
  tot3_int->Draw("PEsame");


  leg->Draw("same");

  if (doPu) latSmall.DrawLatex(0.5,1.1*tot3_int->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU140");
  else latSmall.DrawLatex(0.5,1.1*tot3_int->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU0");
  myc[12]->Update();
  if (doPu) myc[12]->Print("PLOTS_R9/e3x3L1-5overee3x3integrated_Hgg_140pu.pdf");
  else myc[12]->Print("PLOTS_R9/e3x3L1-5overee3x3integrated_Hgg_0pu.pdf");


  myc[13]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  
  TH1F *htot37 = new TH1F("htot37",";e_{reco}/e_{true};Photons",100,0,1.2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1>>htot37",(cut1pca+" && calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)<0.045").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2>>+htot37",(cut2pca+" && calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2)/calibratedE(eSR3Reco2,etaPCA2)<0.045").c_str());
  TH1F *h3337 = (TH1F*)htot37->Clone("h3337");
  h3337->SetFillColor(6);
  h3337->SetFillStyle(3004);
  h3337->SetLineColor(6);
  //h33->Fit("gaus","","same",0.95,1.05);

  std::cout << " -- unconv photon " << htot37->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/eTrue1>>+htot37",(cut1pca+" && calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)>=0.045").c_str());
  t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/eTrue2>>+htot37",(cut2pca+" && calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2)/calibratedE(eSR3Reco2,etaPCA2)>=0.045").c_str());
  
  std::cout << " -- all photons " << htot37->GetEntries() << std::endl;
  htot37->Draw();
  //htot37->Fit("gaus","","same",0.95,1.05);
  h3337->Draw("same");

  lat.DrawLatex(0.1,0.9*htot37->GetMaximum(),"#frac{e_{3#times3}(L1-5)}{e_{3#times3}}<0.045: e_{reco}=e_{3#times3}");
  lat.DrawLatex(0.1,0.75*htot37->GetMaximum(),"#frac{e_{3#times3}(L1-5)}{e_{3#times3}}>=0.045: e_{reco}=e_{7#times7}");

  sprintf(buf,"N_{#gamma}^{unconv}/N_{tot}= %3.2f",h3337->GetEntries()/htot37->GetEntries());
  lat.DrawLatex(0.1,0.6*htot37->GetMaximum(),buf);

  sigeff = effSigmaMacro(h3337,err);
  sprintf(buf,"#sigma_{eff}^{unconv}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(6);
  lat.DrawLatex(0.1,0.45*htot37->GetMaximum(),buf);
  sigeff = effSigmaMacro(htot37,err);
  sprintf(buf,"#sigma_{eff}^{all}=%3.2f #pm %3.2f %%",sigeff*100,err*100);
  lat.SetTextColor(1);
  lat.DrawLatex(0.1,0.3*htot37->GetMaximum(),buf);
 
  myc[13]->Update();
  if (doPu) myc[13]->Print("PLOTS_R9/eReco37overeTrue_convFirst5Layers_Hgg_140pu.pdf");
  else myc[13]->Print("PLOTS_R9/eReco37overeTrue_convFirst5Layers_Hgg_0pu.pdf");
    


  myc[14]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  gPad->SetRightMargin(0.15);
  TH2F *h2Dtot7_33l13o33 = new TH2F("h2Dtot7_33l13o33",";e_{3#times3}(L1-5)/e_{3#times3};e_{reco}/e_{true};Photons",50,0.,0.1,50,0,1.2);
  t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1:calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)>>h2Dtot7_33l13o33",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)>0.99").c_str());
  t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2:calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2,3)/calibratedE(eSR3Reco2,etaPCA2)>>+h2Dtot7_33l13o33",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)>0.99").c_str());
  //TH2F *h2D337_33l13o33 = (TH2F*)h2Dtot7_33l13o33->Clone("h2D337_33l13o33");

  //std::cout << " -- unconv photon " << h2Dtot7_33l13o33->GetEntries() << std::endl;
  //t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/eTrue1:calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)>>+h2Dtot7_33l13o33",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)<=0.99").c_str());
  //t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/eTrue2:calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2,3)/calibratedE(eSR3Reco2,etaPCA2)>>+h2Dtot7_33l13o33",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)<=0.99").c_str());
  
  std::cout << " -- all photons " << h2Dtot7_33l13o33->GetEntries() << std::endl;
  gPad->SetLogz(1);
  h2Dtot7_33l13o33->GetZaxis()->SetTitleOffset(0.7);
  h2Dtot7_33l13o33->Draw("colz");
  (h2Dtot7_33l13o33->ProfileX("h2Dtot7_33l13o33_pfx",h2Dtot7_33l13o33->GetYaxis()->FindBin(0.9),h2Dtot7_33l13o33->GetYaxis()->FindBin(1.2)))->Draw("same");
  TProfile *h7_33l13o33_pfx = (TProfile*)gDirectory->Get("h2Dtot7_33l13o33_pfx");
  h7_33l13o33_pfx->SetMarkerStyle(22);
  h7_33l13o33_pfx->SetMarkerColor(1);
  h7_33l13o33_pfx->Draw("PEsame");
  h7_33l13o33_pfx->Fit("pol2","RI","same",0.,0.1);
  //fit = h7_33l13o33_pfx->GetFunction("pol1");
  //fit->Draw("same");
  //h2D337_33l13o33->SetMarkerStyle(1);
  //h2D337_33l13o33->SetMarkerColor(6);
  //h2D337_33l13o33->SetLineColor(6);
  // h2D337_33l13o33->Draw("P");//same");
  //(h2D337_33l13o33->ProfileX())->Draw("same");
  //TProfile *h7uc_33l13o33_pfx = (TProfile*)gDirectory->Get("h2D337_33l13o33_pfx");
  //h7uc_33l13o33_pfx->SetMarkerStyle(23);
  //h7uc_33l13o33_pfx->SetMarkerColor(6);
  //h7uc_33l13o33_pfx->Draw("PEsame");
  //h7uc_33l13o33_pfx->Fit("pol2","RI","same",0.,0.1);

  TF1 *fitnopu = new TF1("fitnopu","1.017-1.29*x+8.5*x*x",0,0.1);
  fitnopu->SetLineColor(7);
  fitnopu->SetLineWidth(2);
  fitnopu->Draw("same");

  myc[14]->Update();
  if (doPu) myc[14]->Print("PLOTS_R9/eReco37overeTruevse33L1-5oe33_unconv5x5_Hgg_140pu.pdf");
  else myc[14]->Print("PLOTS_R9/eReco37overeTruevse33L1-5oe33_unconv5x5_Hgg_0pu.pdf");
    

  myc[15]->cd();
  //gStyle->SetOptStat("eMRuo");
  //gStyle->SetOptFit(1111);
  gPad->SetRightMargin(0.15);
  TH2F *h2Dtot7_33l13o33conv = new TH2F("h2Dtot7_33l13o33conv",";e_{3#times3}(L1-5)/e_{3#times3};e_{reco}/e_{true};Photons",50,0.,0.1,50,0,1.2);
  //t->Draw("calibratedE(eSR3Reco1,etaPCA1)/eTrue1:calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)>>h2Dtot7_33l13o33conv",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)>0.99").c_str());
  //t->Draw("calibratedE(eSR3Reco2,etaPCA2)/eTrue2:calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2,3)/calibratedE(eSR3Reco2,etaPCA2)>>+h2Dtot7_33l13o33conv",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)>0.99").c_str());
  //TH2F *h2D337_33l13o33conv = (TH2F*)h2Dtot7_33l13o33conv->Clone("h2D337_33l13o33conv");

  //std::cout << " -- unconv photon " << h2Dtot7_33l13o33conv->GetEntries() << std::endl;
  t->Draw("calibratedE(eSR7Reco1,etaPCA1,7)/eTrue1:calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaPCA1)/calibratedE(eSR3Reco1,etaPCA1)>>+h2Dtot7_33l13o33conv",(cut1pca+" && calibratedE(eSR3Reco1,etaPCA1)/calibratedE(eSR5Reco1,etaPCA1,5)<=0.99").c_str());
  t->Draw("calibratedE(eSR7Reco2,etaPCA2,7)/eTrue2:calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaPCA2,3)/calibratedE(eSR3Reco2,etaPCA2)>>+h2Dtot7_33l13o33conv",(cut2pca+" && calibratedE(eSR3Reco2,etaPCA2)/calibratedE(eSR5Reco2,etaPCA2,5)<=0.99").c_str());
  
  std::cout << " -- all photons " << h2Dtot7_33l13o33conv->GetEntries() << std::endl;
  gPad->SetLogz(1);
  h2Dtot7_33l13o33conv->GetZaxis()->SetTitleOffset(0.7);
  h2Dtot7_33l13o33conv->Draw("colz");
  (h2Dtot7_33l13o33conv->ProfileX("h2Dtot7_33l13o33conv_pfx",h2Dtot7_33l13o33conv->GetYaxis()->FindBin(0.9),h2Dtot7_33l13o33conv->GetYaxis()->FindBin(1.2)))->Draw("same");
  TProfile *h7_33l13o33conv_pfx = (TProfile*)gDirectory->Get("h2Dtot7_33l13o33conv_pfx");
  h7_33l13o33conv_pfx->SetMarkerStyle(22);
  h7_33l13o33conv_pfx->SetMarkerColor(1);
  h7_33l13o33conv_pfx->Draw("PEsame");
  h7_33l13o33conv_pfx->Fit("pol1","RI","same",0.,0.1);
  //fit = h7_33l13o33conv_pfx->GetFunction("pol1");
  //fit->Draw("same");
  //h2D337_33l13o33conv->SetMarkerStyle(1);
  //h2D337_33l13o33conv->SetMarkerColor(6);
  //h2D337_33l13o33conv->SetLineColor(6);
  //h2D337_33l13o33conv->Draw("P");//same");
  //(h2D337_33l13o33conv->ProfileX())->Draw("same");
  //TProfile *h7uc_33l13o33conv_pfx = (TProfile*)gDirectory->Get("h2D337_33l13o33conv_pfx");
  //h7uc_33l13o33conv_pfx->SetMarkerStyle(23);
  //h7uc_33l13o33conv_pfx->SetMarkerColor(6);
  //h7uc_33l13o33conv_pfx->Draw("PEsame");
  //h7uc_33l13o33conv_pfx->Fit("pol2","RI","same",0.,0.1);
  TF1 *fitnopuconv = new TF1("fitnopuconv","0.993-0.60*x",0,0.1);
  fitnopuconv->SetLineColor(7);
  fitnopuconv->SetLineWidth(2);
  fitnopuconv->Draw("same");

  myc[15]->Update();
  if (doPu) myc[15]->Print("PLOTS_R9/eReco37overeTruevse33L1-5oe33_conv5x5_Hgg_140pu.pdf");
  else myc[15]->Print("PLOTS_R9/eReco37overeTruevse33L1-5oe33_conv5x5_Hgg_0pu.pdf");
    




  return 0;
}
