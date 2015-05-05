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

int plotReso() {

  std::string outfile = "InfoPlots_0_7x7";
  //std::string outfile = "InfoPlots_140";
  //std::string outfile = "InfoPlots_singleGamma";

  //std::string cut1 = "converted1==0 && eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && showerMax1>5 && showerMax1<25 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11)";// && eSR3Reco1/eTrue1>0.85";// && TMath::Abs(etaReco1-etaTrue1)<0.1 && eSeed1/eSC1>0.95 && eSeed1/eSC1<1.01 && dRTrue1<0.05";

  std::string cut1 = "converted1==0 && eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetid1==0 && invalidNeighbour1==0";
  std::string cut2 = "converted2==0 && eTrue2/cosh(etaTrue2)>40 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetid2==0 && invalidNeighbour2==0";
  std::string cut1pca = "calibratedE(eSR3Reco1,etaPCA1)/eSC1>0.9 && eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetidPCA1==0 && invalidNeighbourPCA1==0 && dRTrue1 < 0.05";
  std::string cut2pca = "calibratedE(eSR3Reco2,etaPCA2)/eSC2>0.9 && eTrue2/cosh(etaTrue2)>40 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetidPCA2==0 && invalidNeighbourPCA2==0 && dRTrue2 < 0.05";

  TFile *fin = TFile::Open("Truth_Hgg_0pu.root");
  //TFile *fin = TFile::Open("test.root");
  if (!fin) return 1;
  fin->cd("hgg");

  TTree *t = (TTree*)gDirectory->Get("tree");

  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }

  const unsigned nV = 30;
  std::string var[nV] = {
    "xvtxTrue","yvtxTrue","zvtxTrue",
    "dRTrue",
    "eSR7RecoOvereTrue","eSR7RecoOvereTrue",
    "etaSCMinusetadetTrue","etaRecoMinusetaTrue","etaPCAMinusetaTrue",
    "phiSCMinusphiTrue","phiRecoMinusphiTrue","phiPCAMinusphiTrue",
    "etaWidth","phiWidth","clustersSize","nClusters09",
    "eSeedOvereSC","showerMax",
    "eSR3RecoOvereTrue",
    "xPCA","yPCA","zPCA","eSR7RecoOvereTrue","eSR7RecoPCAOvereTrue","eSR7RecoPCAOvereTrue",
    "isValid","invalidDetid","invalidNeighbour","noPhiCrack","noPhiCrackPCA"
  };
  std::string var1[nV];
  std::string var2[nV];
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    size_t over=var[iV].find("Over");
    size_t minus=var[iV].find("Minus");
    size_t end=var[iV].npos;
    size_t size=var[iV].size();
    //std::cout << var[iV] << " pos(over) " << over << " pos(minus) " << minus << " end " << size << std::endl;
    if (over==end && 
	minus==end){
      var1[iV] = var[iV]+"1";
      var2[iV] = var[iV]+"2";
    }
    else {
      if (over!=end){
	var1[iV]=var[iV].substr(0,over)+"1/"+var[iV].substr(over+4,end)+"1";
	var2[iV]=var[iV].substr(0,over)+"2/"+var[iV].substr(over+4,end)+"2";
      }
      else if (minus!=end){
	var1[iV]=var[iV].substr(0,minus)+"1-"+var[iV].substr(minus+5,end)+"1";
	var2[iV]=var[iV].substr(0,minus)+"2-"+var[iV].substr(minus+5,end)+"2";
      }
    }
    std::cout << " Vars " << var[iV] << " " << var1[iV] << " " << var2[iV] << std::endl;
  }

  //z=r*costheta
  //x=r*sintheta*cosphi
  //y=r*sintheta*sinphi
  //theta = 2*TMath::ATan(exp(-1.*etaTrue1))
  var1[19] = var1[19]+"-zPCA1/TMath::Cos(2*TMath::ATan(exp(-1.*etaTrue1)))*TMath::Sin(2*TMath::ATan(exp(-1.*etaTrue1)))*TMath::Cos(phiTrue1)";
  var2[19] = var2[19]+"-zPCA2/TMath::Cos(2*TMath::ATan(exp(-1.*etaTrue2)))*TMath::Sin(2*TMath::ATan(exp(-1.*etaTrue2)))*TMath::Cos(phiTrue2)";
  var1[20] = var1[20]+"-zPCA1/TMath::Cos(2*TMath::ATan(exp(-1.*etaTrue1)))*TMath::Sin(2*TMath::ATan(exp(-1.*etaTrue1)))*TMath::Sin(phiTrue1)";
  var2[20] = var2[20]+"-zPCA2/TMath::Cos(2*TMath::ATan(exp(-1.*etaTrue2)))*TMath::Sin(2*TMath::ATan(exp(-1.*etaTrue2)))*TMath::Sin(phiTrue2)";

  const unsigned nMore = 6;
  const unsigned nC = nV+nMore;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "myc" << ic ;
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1);
  }
  myc[0]->cd();
  myc[0]->Print((outfile+".pdf[").c_str());

  gStyle->SetOptStat("e");
  gStyle->SetOptFit(1111);

  //TH2F *h2 = new TH2F("h2",";egen (GeV);ereco 3x3 (mips); photons",60,0,600,1000,0,50000);
  //t->Draw("eSR3Reco1:eTrue1>>h2",cut1.c_str());
  //t->Draw("eSR3Reco2:eTrue2>>+h2",cut2.c_str());
  TH2F *h2 = new TH2F("h2",";egen (GeV);ereco 3x3 (mips); photons",60,0,600,60,0,600);
  t->Draw("calibratedE(eSR7Reco1,etaTrue1):eTrue1>>h2",cut1.c_str());
  t->Draw("calibratedE(eSR7Reco2,etaTrue2):eTrue2>>+h2",cut2.c_str());
  myc[0]->Clear();
  myc[0]->Divide(1,2);

  h2->Draw("colz");

  h2->ProfileX();
  myc[0]->Update();
  TProfile *h2_pfx = (TProfile*)gDirectory->Get("h2_pfx");
  h2_pfx->SetMarkerStyle(22);
  h2_pfx->SetMarkerColor(1);
  h2_pfx->Draw("PEsame");
  h2_pfx->Fit("pol1","RI","same",100,400);

  TF1 *fit3 = h2_pfx->GetFunction("pol1");
  //if (!fit3) return 1;
  double slope3 = fit3?fit3->GetParameter(1):0;
  double offset3 = 0;//fit3->GetParameter(0);

  TLatex lat;
  char buf[100];
  double max = h2->GetYaxis()->GetBinLowEdge(h2->GetYaxis()->GetNbins());
  sprintf(buf,"a = %3.2f #pm %3.2f",fit3?fit3->GetParameter(0):0,fit3?fit3->GetParError(0):0);
  lat.DrawLatex(50,max*0.8,"Ereco = a + b #times Egen");
  lat.DrawLatex(50,max*0.7,buf);
  sprintf(buf,"b = %3.2f #pm %3.2f",slope3,fit3?fit3->GetParError(1):0);
  lat.DrawLatex(50,max*0.6,buf);

  myc[0]->Update();
  myc[0]->Print((outfile+".pdf").c_str());

  myc[1]->cd();
  TH2F *h2eta = new TH2F("h2eta",";#eta;ereco/egen; photons",60,-3,3,100,0,2);

  std::ostringstream cor;
  cor.str("");
  //cor << "(eSR7Reco1-" << offset3 << ")/(" << slope3 << "*eTrue1):TMath::Abs(etaTrue1)>>h2eta";
  //cor << "(calibratedE(eSR7Reco1,etaTrue1)-" << offset3 << ")/(" << slope3 << "*eTrue1):TMath::Abs(etaTrue1)>>h2eta";
   cor << "calibratedE(eSR7Reco1,etaTrue1)/eTrue1:etaTrue1>>h2eta";
  t->Draw(cor.str().c_str(),cut1.c_str());
  cor.str("");
  //cor << "(eSR7Reco2-" << offset3 << ")/(" << slope3 << "*eTrue2):TMath::Abs(etaTrue2)>>+h2eta";
  //cor << "(calibratedE(eSR7Reco2,etaTrue2)-" << offset3 << ")/(" << slope3 << "*eTrue2):TMath::Abs(etaTrue2)>>+h2eta";
  cor << "calibratedE(eSR7Reco2,etaTrue2)/eTrue2:etaTrue2>>+h2eta";
  t->Draw(cor.str().c_str(),cut2.c_str());
  //h2eta->GetYaxis()->SetRangeUser(0.8,1.2);
  h2eta->Draw("colz");

  h2eta->ProfileX();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  myc[1]->Update();
  TProfile *h2eta_pfx = (TProfile*)gDirectory->Get("h2eta_pfx");
  h2eta_pfx->SetMarkerStyle(22);
  h2eta_pfx->SetMarkerColor(1);
  h2eta_pfx->Draw("PEsame");
  h2eta_pfx->Fit("pol1","","same");

  TF1 *fiteta3 = h2eta_pfx->GetFunction("pol1");
  //if (!fiteta3) return 1;

  double offset3eta = fiteta3?fiteta3->GetParameter(0):0;
  double slope3eta = fiteta3?fiteta3->GetParameter(1):1;
  //double p23 = fit3->GetParameter(2);
  max = 1.2;//h2eta->GetYaxis()->GetBinLowEdge(h2eta->GetYaxis()->GetNbins());
  sprintf(buf,"a = %3.2f #pm %3.2f",fiteta3?fiteta3->GetParameter(0):0,fiteta3?fiteta3->GetParError(0):0);
  lat.DrawLatex(50,max*0.8,"Ereco = a + b #times Egen");
  lat.DrawLatex(50,max*0.7,buf);
  sprintf(buf,"b = %3.2f #pm %3.2f",slope3eta,fiteta3?fiteta3->GetParError(1):0);
  lat.DrawLatex(50,max*0.6,buf);

  myc[1]->Update();
  myc[1]->Print((outfile+".pdf").c_str());

  myc[2]->cd();
  TH2F *h2pca = new TH2F("h2pca",";egen (GeV);ereco 3x3 (mips); photons",60,0,600,60,0,600);
  t->Draw("calibratedE(eSR7RecoPCA1,etaPCA1):eTrue1>>h2pca",cut1pca.c_str());
  t->Draw("calibratedE(eSR7RecoPCA2,etaPCA2):eTrue2>>+h2pca",cut2pca.c_str());
  //TH2F *h2pca = new TH2F("h2pca",";egen (GeV);ereco 3x3 (mips); photons",60,0,600,1000,0,50000);
  //t->Draw("eSR7RecoPCA1:eTrue1>>h2pca",cut1.c_str());
  //t->Draw("eSR7RecoPCA2:eTrue2>>+h2pca",cut2.c_str());
  h2pca->Draw("colz");

  h2pca->ProfileX();
  myc[2]->Update();
  TProfile *h2pca_pfx = (TProfile*)gDirectory->Get("h2pca_pfx");
  h2pca_pfx->SetMarkerStyle(22);
  h2pca_pfx->SetMarkerColor(1);
  h2pca_pfx->Draw("PEsame");
  h2pca_pfx->Fit("pol1","RI","same",100,400);

  TF1 *fit3pca = h2pca_pfx->GetFunction("pol1");
  //if (!fit3pca) return 1;
  //double slope3pca = fit3pca->GetParameter(1);
  //double offset3pca = 0;//fit3pca->GetParameter(0);

  //max = h2pca->GetYaxis()->GetBinLowEdge(h2pca->GetYaxis()->GetNbins());
  //sprintf(buf,"a = %3.2f #pm %3.2f",fit3pca->GetParameter(0),fit3pca->GetParError(0));
  //lat.DrawLatex(50,max*0.8,"Ereco = a + b #times Egen");
  //lat.DrawLatex(50,max*0.7,buf);
  //sprintf(buf,"b = %3.2f #pm %3.2f",slope3pca,fit3pca->GetParError(1));
  // lat.DrawLatex(50,max*0.6,buf);

  myc[2]->Update();
  myc[2]->Print((outfile+".pdf").c_str());


  myc[3]->cd();
  TH2F *h2pcaeta = new TH2F("h2pcaeta",";#eta;ereco/egen; photons",60,-3,3,100,0,2);

  cor.str("");
  //cor << "(eSR7RecoPCA1-" << offset3pca << ")/(" << slope3pca << "*eTrue1):TMath::Abs(etaTrue1)>>h2pcaeta";
  cor << "calibratedE(eSR7RecoPCA1,etaPCA1)/eTrue1:etaTrue1>>h2pcaeta";
  t->Draw(cor.str().c_str(),cut1pca.c_str());
  cor.str("");
  //cor << "(eSR7RecoPCA2-" << offset3pca << ")/(" << slope3pca << "*eTrue2):TMath::Abs(etaTrue2)>>+h2pcaeta";
  cor << "calibratedE(eSR7RecoPCA2,etaPCA2)/eTrue2:etaTrue2>>+h2pcaeta";
  t->Draw(cor.str().c_str(),cut2pca.c_str());
  h2pcaeta->Draw("colz");

  h2pcaeta->ProfileX();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  myc[3]->Update();
  TProfile *h2pcaeta_pfx = (TProfile*)gDirectory->Get("h2pcaeta_pfx");
  h2pcaeta_pfx->SetMarkerStyle(22);
  h2pcaeta_pfx->SetMarkerColor(1);
  h2pcaeta_pfx->Draw("PEsame");
  h2pcaeta_pfx->Fit("pol1","","same");

  //TF1 *fiteta3pca = h2pcaeta_pfx->GetFunction("pol1");
  //if (!fiteta3pca) return 1;

  //double p03 = fit3->GetParameter(0);
  //double p13 = fit3->GetParameter(1);
  //double p23 = fit3->GetParameter(2);

  myc[3]->Update();
  myc[3]->Print((outfile+".pdf").c_str());


  /*
  TH2F *h25 = new TH2F("h25",";egen (GeV);ereco 5x5 (mips); photons",100,0,1000,100,0,1000);
  t->Draw("eSR7Reco1:eTrue1>>h25",cut1.c_str());
  t->Draw("eSR7Reco2:eTrue2>>+h25",cut2.c_str());
  h25->Draw("colz");

  h25->ProfileX();
  myc[2]->Update();
  TProfile *h25_pfx = (TProfile*)gDirectory->Get("h25_pfx");
  h25_pfx->SetMarkerStyle(22);
  h25_pfx->SetMarkerColor(1);
  h25_pfx->Draw("PEsame");
  h25_pfx->Fit("pol1","RL","same",80,200);

  TF1 *fit5 = h25_pfx->GetFunction("pol1");
  if (!fit5) return 1;
  double slope5 = fit5->GetParameter(1);
  double offset5 = fit5->GetParameter(0);

  myc[2]->Update();
  myc[2]->Print((outfile+".pdf").c_str());

  myc[3]->cd();

  TH2F *h25eta = new TH2F("h25eta",";#eta;ereco 5x5/egen; photons",30,1.5,3,100,0,2);
  cor.str("");
  cor << "(eSR7Reco1-" << offset5 << ")/(" << slope5 << "*eTrue1):TMath::Abs(etaTrue1)>>h25eta";
  t->Draw(cor.str().c_str(),cut1.c_str());
  cor.str("");
  cor << "(eSR7Reco2-" << offset5 << ")/(" << slope5 << "*eTrue2):TMath::Abs(etaTrue2)>>+h25eta";
  t->Draw(cor.str().c_str(),cut2.c_str());
  //  t->Draw("eSR7Reco1/eTrue1:TMath::Abs(etaTrue1)>>h25eta",cut1.c_str());
  //t->Draw("eSR7Reco2/eTrue2:TMath::Abs(etaTrue2)>>htmp5eta",cut2.c_str());
  h25eta->Draw("colz");

  h25eta->ProfileX();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  myc[3]->Update();
  TProfile *h25eta_pfx = (TProfile*)gDirectory->Get("h25eta_pfx");
  h25eta_pfx->SetMarkerStyle(22);
  h25eta_pfx->SetMarkerColor(1);
  h25eta_pfx->Draw("PEsame");
  h25eta_pfx->Fit("pol2","","same");

  TF1 *fiteta5 = h25eta_pfx->GetFunction("pol2");
  if (!fiteta5) return 1;
  double p05 = fit5->GetParameter(0);
  double p15 = fit5->GetParameter(1);
  double p25 = fit5->GetParameter(2);
  myc[3]->Update();
  myc[3]->Print((outfile+".pdf").c_str());
  */

  myc[4]->cd();
  //gPad->SetGridx(1);
  TH2F *h2phi = new TH2F("h2phi",";#phi;ereco 3x3/egen; photons",360,-3.1416,3.1416,100,0,2);
  //t->Draw("eSR7Reco1/eTrue1:phiTrue1>>h2phi",cut1.c_str());
  //t->Draw("eSR7Reco2/eTrue2:phiTrue2>>+h2phi",cut2.c_str());
  cor.str("");
  //cor << "((eSR7Reco1-" << offset3 << ")/(" << slope3 << "*eTrue1)-" << offset3eta << ")/" << slope3eta << ":phiTrue1>>h2phi";
  //cor << "(calibratedE(eSR7Reco1,etaTrue1)-" << offset3 << ")/(" << slope3 << "*eTrue1):phiTrue1>>h2phi";
  cor << "calibratedE(eSR7Reco1,etaTrue1)/eTrue1:phiTrue1>>h2phi";
  t->Draw(cor.str().c_str(),cut1.c_str());
  cor.str("");
  //cor << "((eSR7Reco2-" << offset3 << ")/(" << slope3 << "*eTrue2)-" << offset3eta << ")/" << slope3eta << ":phiTrue2>>+h2phi";
  //cor << "(calibratedE(eSR7Reco2,etaTrue2)-" << offset3 << ")/(" << slope3 << "*eTrue2):phiTrue2>>+h2phi";
  cor << "calibratedE(eSR7Reco2,etaTrue2)/eTrue2:phiTrue2>>+h2phi";
  t->Draw(cor.str().c_str(),cut2.c_str());

  h2phi->Draw("colz");
  h2phi->ProfileX();
  TProfile *h2phi_pfx = (TProfile*)gDirectory->Get("h2phi_pfx");
  h2phi_pfx->SetMarkerStyle(22);
  h2phi_pfx->SetMarkerColor(1);
  h2phi_pfx->Draw("PEsame");

  TLine *l[18];
  for (unsigned is(0); is<18;++is){
    int crack = -170+20*is;
    l[is] = new TLine(crack*2*TMath::Pi()/360.,0,crack*2*TMath::Pi()/360.,2);
    l[is]->Draw();
  }

  myc[4]->Update();
  myc[4]->Print((outfile+".pdf").c_str());
  myc[4]->Print((outfile+"_erecooveregenvsphi.pdf").c_str());

  myc[5]->cd();
  TH2F *h2phipca = new TH2F("h2phipca",";#phi;ereco 3x3/egen; photons",360,-3.1416,3.1416,100,0,2);
  cor.str("");
  cor << "calibratedE(eSR7RecoPCA1,etaPCA1)/eTrue1:phiTrue1>>h2phipca";
  t->Draw(cor.str().c_str(),cut1pca.c_str());
  cor.str("");
  cor << "calibratedE(eSR7RecoPCA2,etaPCA2)/eTrue2:phiTrue2>>+h2phipca";
  t->Draw(cor.str().c_str(),cut2pca.c_str());

  h2phipca->Draw("colz");
  h2phipca->ProfileX();
  TProfile *h2phipca_pfx = (TProfile*)gDirectory->Get("h2phipca_pfx");
  h2phipca_pfx->SetMarkerStyle(22);
  h2phipca_pfx->SetMarkerColor(1);
  h2phipca_pfx->Draw("PEsame");

  for (unsigned is(0); is<18;++is){
    l[is]->Draw();
  }

  myc[5]->Update();
  myc[5]->Print((outfile+".pdf").c_str());

  cor.str("");
  //cor << "(eSR7RecoPCA1-" << offset3pca << ")/(" << slope3pca << "*eTrue1)";
  cor << "calibratedE(eSR7Reco1,etaTrue1)/eTrue1";
  var1[5] = cor.str();
  
  cor.str("");
  //cor << "((eSR7Reco1-" << offset3 << ")/(" << slope3 << "*eTrue1)-" << offset3eta << ")/" << slope3eta;
  //cor << "(calibratedE(eSR7Reco1,etaTrue1)-" << offset3 << ")/(" << slope3 << "*eTrue1)";
  cor << "calibratedE(eSR7Reco1,etaTrue1)/eTrue1";
  var1[4] = cor.str();
  
  cor.str("");
  //cor << "((eSR7Reco1-" << offset3 << ")/(" << slope3 << "*eTrue1)-" << offset3eta << ")/" << slope3eta;
  //cor << "(calibratedE(eSR7Reco1,etaTrue1)-" << offset3 << ")/(" << slope3 << "*eTrue1)";
  cor << "calibratedE(eSR7RecoPCA1,etaPCA1)/eTrue1";
  var1[23] = cor.str();
  
  cor.str("");
  //cor << "(eSR7RecoPCA2-" << offset3pca << ")/(" << slope3pca << "*eTrue2)";
  cor << "calibratedE(eSR7Reco2,etaTrue2)/eTrue2";
  var2[5] = cor.str();

  cor.str("");
  //cor << "((eSR7Reco2-" << offset3 << ")/(" << slope3 << "*eTrue2)-" << offset3eta << ")/" << slope3eta;
  //cor << "(calibratedE(eSR7Reco2,etaTrue2)-" << offset3 << ")/(" << slope3 << "*eTrue2)";
  cor << "calibratedE(eSR7Reco2,etaTrue2)/eTrue2";
  var2[4] = cor.str();

  cor.str("");
  //cor << "((eSR7Reco2-" << offset3 << ")/(" << slope3 << "*eTrue2)-" << offset3eta << ")/" << slope3eta;
  //cor << "(calibratedE(eSR7Reco2,etaTrue2)-" << offset3 << ")/(" << slope3 << "*eTrue2)";
  cor << "calibratedE(eSR7RecoPCA2,etaPCA2)/eTrue2";
  var2[23] = cor.str();

  var[4]+="Calib";
  var[5]+="Calib";
  var[23]+="Calib";

  //cor.str("");
  //cor << "(eSR7Reco2-" << offset5 << ")/(" << slope5 << "*eTrue2)";
  //var2[5] = cor.str();

  TH1F *hist[nV];

  gStyle->SetOptStat("eMRuo");
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    myc[nMore+iV]->cd();
    if (var1[iV].find("PCA")!=var1[iV].npos) t->Draw(var1[iV].c_str(),cut1pca.c_str());
    else t->Draw(var1[iV].c_str(),cut1.c_str());
    hist[iV] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(var[iV].c_str());
    std::cout << " -- first photon " << hist[iV]->GetEntries() << std::endl;
    cor.str("");
    cor << var2[iV] << ">>+" << var[iV];
    std::cout << "debug " << cor.str() << std::endl;
    if (var1[iV].find("PCA")!=var1[iV].npos) t->Draw(cor.str().c_str(),cut2pca.c_str());
    else t->Draw(cor.str().c_str(),cut2.c_str());
    std::cout << " -- second photon " << hist[iV]->GetEntries() << std::endl;
    std::cout << " -- Hist " << var[iV] << " " << hist[iV]->GetEntries() << " " << hist[iV]->GetMean() << " " << hist[iV]->GetRMS() << std::endl;
    hist[iV]->SetTitle("");
    hist[iV]->GetXaxis()->SetTitle(var[iV].c_str());
    if (iV==19) hist[iV]->GetXaxis()->SetTitle("xPCA-xTrue");
    if (iV==20) hist[iV]->GetXaxis()->SetTitle("yPCA-yTrue");
    hist[iV]->GetYaxis()->SetTitle("Photons");
    myc[nMore+iV]->Clear();
    myc[nMore+iV]->cd();
    //if (iV==4) hist[iV]->Rebin(2); 
    hist[iV]->Draw();
    gStyle->SetOptFit(1111);
    if (iV==4 || iV==5 || iV==23) hist[iV]->Fit("gaus","LR","same",0.95,1.05);
    myc[nMore+iV]->Update();
    myc[nMore+iV]->Print((outfile+".pdf").c_str());
    
    //if (iV==3) return 1;
  }//loop on variables
  
  myc[0]->Print((outfile+".pdf]").c_str());


  return 0;
}
