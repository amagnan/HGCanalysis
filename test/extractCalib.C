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

int fitEnergy(TTree* t, TCanvas * & myc,const std::vector<int> & evec,
	       const std::string & cut1,const std::string & cut2,
	      const std::string & var,
	      const double & offset, const double & slope,
	      TGraphErrors* & grE3,TGraphErrors* & grE3reso,
	      bool calib){


  TH1F *hist[evec.size()];
  myc->Divide(evec.size()/3,3);
  unsigned ip = 0;
  for (unsigned ie(0); ie<evec.size();++ie){
    myc->cd(ie+1);
    std::ostringstream label;
    label << "hist" << ie;
    //hist[ie] = new TH2F(label.str().c_str(),";#eta;ereco 3x3 (mips); photons",100,-3,3,1000,0,50000);
    //hist[ie] = new TH1F(label.str().c_str(),";ereco 3x3 (mips); photons",1000,0,30000);
    std::ostringstream cut;
    cut << cut1 << "&& TMath::Nint(eTrue1)==" << evec[ie];
    //t->Draw(("eSR3Reco1:etaTrue1>>"+label.str()).c_str(),cut.str().c_str());
    std::ostringstream varname;
    varname << "(" << var << "1-" << offset << ")/" << slope;
    t->Draw(varname.str().c_str(),cut.str().c_str());
    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
    if (!htemp) {
      std::cout << " -- energy " << evec[ie] << " no entries found." << std::endl;
      continue;
    }
    hist[ie] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(label.str().c_str());
    cut.str(""); 
    cut << cut2 << "&& TMath::Nint(eTrue2)==" << evec[ie];
    //t->Draw(("eSR3Reco2:etaTrue2>>+"+label.str()).c_str(),cut.str().c_str());
    varname.str("");
    varname << "(" << var << "2-" << offset << ")/" << slope << ">>+" << label.str();

    t->Draw(varname.str().c_str(),cut.str().c_str());
    if (hist[ie]->GetEntries()==0) {
      std::cout << " -- energy " << evec[ie] << " no entries found." << std::endl;
      continue;
    }
    hist[ie]->Fit("gaus","LIQ");
    TF1 *fit3 = hist[ie]->GetFunction("gaus");
    if (!fit3) {
      std::cout << " -- fit failed for " << var << " E=" << evec[ie] << std::endl;
      continue;
    }
    hist[ie]->Draw();
    hist[ie]->Fit("gaus","LIRQ","same",
		  fit3->GetParameter(1)-2*fit3->GetParameter(2),
		  fit3->GetParameter(1)+2*fit3->GetParameter(2));
    fit3 = hist[ie]->GetFunction("gaus");
    if (!fit3) {
      std::cout << " -- fit failed for " << var << " E=" << evec[ie] << std::endl;
      continue;
    }
    //std::cout << evec[ie] << " " << fit3->GetParameter(1) << " " << fit3->GetParameter(2) << std::endl;
    double reso = fit3->GetParameter(2)/fit3->GetParameter(1);
    grE3reso->SetPoint(ip,evec[ie],reso);
    double err = reso*sqrt(pow(fit3->GetParError(1)/fit3->GetParameter(1),2)+pow(fit3->GetParError(2)/fit3->GetParameter(2),2));
    grE3reso->SetPointError(ip,0,err);
    if (!calib) {
      grE3->SetPoint(ip,evec[ie],fit3->GetParameter(1));
      grE3->SetPointError(ip,0,fit3->GetParError(1));
    } else {
      double res = (fit3->GetParameter(1)-evec[ie])/evec[ie];
      double reserr = fit3->GetParError(1)/evec[ie];
      grE3->SetPoint(ip,evec[ie],res);
      grE3->SetPointError(ip,0,reserr);
    }
    ip++;
    //hist[ie]->Draw("colz");
  }

  return 0;
};
int calibAndReso(TTree* t, TCanvas ** myc,
		 unsigned idx,
		 const std::vector<int> & evec,
		 const std::string & cut1,const std::string & cut2,
		 const std::string & var,
		 double & offset, double & slope,
		 bool calib)
{
  TGraphErrors *grE3 = new TGraphErrors();
  grE3->SetName("grE3");

  TGraphErrors *grE3reso = new TGraphErrors();
  grE3reso->SetName("grE3reso");

  if (fitEnergy(t,myc[idx],evec,cut1,cut2,var,offset,slope,grE3,grE3reso,calib)) return 1;
  myc[idx]->Update();

  myc[idx+1]->cd();
  grE3->SetMarkerStyle(22);
  grE3->SetMarkerColor(2);
  grE3->SetLineColor(2);
  grE3->Draw("APE");
  grE3->GetXaxis()->SetTitle("Egen (GeV)");
  if (!calib) {
    grE3->GetYaxis()->SetTitle("Ereco (Mips)");
  }
  else {
    grE3->GetYaxis()->SetTitle("(Ereco-Egen)/Egen");
    grE3->GetYaxis()->SetRangeUser(-0.1,0.1);
  }
  if (!calib) grE3->Fit("pol1","IR","same",20,300);
  else grE3->Fit("pol1","I","same");
  TF1 *fit = grE3->GetFunction("pol1");
  if (!fit) {
    std::cout << " -- fit failed." << std::endl;
  } else {
    slope = fit->GetParameter(1);
    offset = fit->GetParameter(0);
  }
  myc[idx+1]->Update();

  myc[idx+2]->cd();
  grE3reso->SetMarkerStyle(22);
  grE3reso->SetMarkerColor(2);
  grE3reso->SetLineColor(2);
  grE3reso->Draw("APE");
  grE3reso->GetXaxis()->SetTitle("Egen (GeV)");
  grE3reso->GetYaxis()->SetTitle("sigma(E)/E");
  //grE3->Fit("pol1","I","same");
  //TF1 *fit = grE3->GetFunction("pol1");
  myc[idx+2]->Update();

  return 0;
};

int extractCalib() {

  std::string outfile = "CalibPlots";

  std::string cut1base = "isValid1==1 && converted1==0 && noPhiCrack1==1 && etaTrue1>0 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11)";// && eSR3Reco1/eTrue1>0.85";// && TMath::Abs(etaReco1-etaTrue1)<0.1 && eSeed1/eSC1>0.95 && eSeed1/eSC1<1.01";
  std::string cut2base = "isValid2==1 && converted2==0 && noPhiCrack2==1 && etaTrue2>0 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11)";// && eSR3Reco2/eTrue2>0.85";// && TMath::Abs(etaReco2-etaTrue2)<0.1 && eSeed2/eSC2>0.95 && eSeed2/eSC2<1.01";

  TFile *fin = TFile::Open("Calib_singleGamma_0pu_all.root");
  if (!fin) {
    std::cout << " -- Input file not found!" << std::endl;
    return 1;
  }
  fin->cd("hgg");
  
  TTree *t = (TTree*)gDirectory->Get("tree");
  
  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }
  
  const unsigned nC = 15;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "myc" << ic ;
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1500,1000);
    if (ic>0) {
      myc[ic]->SetGridx(1);
      myc[ic]->SetGridy(1);
    }
  }
  myc[0]->cd();
  myc[0]->Print((outfile+".pdf[").c_str());

  //get list of true energies
  TH1I *hE = new TH1I("hE",";eGen (GeV)",600,0,600);
  t->Draw("TMath::Nint(eTrue1)>>hE",cut1base.c_str());
  t->Draw("TMath::Nint(eTrue2)>>+hE",cut2base.c_str());

  myc[0]->Print((outfile+".pdf").c_str());


  std::vector<int> evec;

  for (int bin(1); bin<hE->GetNbinsX()+1;++bin){
    if (hE->GetBinContent(bin)>0) {
      std::cout << " -- energy " << hE->GetXaxis()->GetBinLowEdge(bin) << " GeV with " << hE->GetBinContent(bin) << " photons." << std::endl;
      evec.push_back(hE->GetXaxis()->GetBinLowEdge(bin));
    }
  }

  std::cout << " -- Processing " << evec.size() << " energy points." << std::endl;

  gStyle->SetOptStat("e");
  gStyle->SetOptFit(1111);
  gStyle->SetStatW(0.25);
  gStyle->SetStatX(0.6);
  gStyle->SetStatY(1);


  const unsigned neta = 6;
  const double etamin = 1.6;
  const double step = 0.2;
  double eta[neta];
  double slope[neta];
  double offset[neta];
  double slopecor[neta];
  double offsetcor[neta];

  for (unsigned ieta(0); ieta<neta;++ieta){
    eta[ieta] = etamin+step*ieta+step/2.;

    std::cout << " -- Processing eta " << eta[ieta] << std::endl;

    slope[ieta] = 1;
    offset[ieta] = 0;
 
    std::ostringstream cut1;
    cut1 << cut1base << " && TMath::Abs(etaTrue1)>" << eta[ieta]-step/2. << " && TMath::Abs(etaTrue1)<=" << eta[ieta]+step/2.;
   std::ostringstream cut2;
    cut2 << cut2base << " && TMath::Abs(etaTrue2)>" << eta[ieta]-step/2. << " && TMath::Abs(etaTrue2)<=" << eta[ieta]+step/2.;

    if (calibAndReso(t,myc,1,evec,cut1.str(),cut2.str(),"eSR3Reco",offset[ieta],slope[ieta],false)) return 1;
    
    myc[1]->Update();
    myc[1]->Print((outfile+".pdf").c_str());
    
    myc[2]->Update();
    myc[2]->Print((outfile+".pdf").c_str());
    
    myc[3]->Update();
    myc[3]->Print((outfile+".pdf").c_str());
    
    std::cout << " -- No eta cor: offset and slope = " << offset[ieta] << " " << slope[ieta] << std::endl;
    
    double offsettmp = offset[ieta];
    double slopetmp = slope[ieta];
    if (calibAndReso(t,myc,4,evec,cut1.str(),cut2.str(),"eSR3Reco",offsettmp,slopetmp,true)) return 1;
    
    std::cout << " -- No eta cor after calibration: offset and slope = " << offsettmp << " " << slopetmp << std::endl;
    
    myc[4]->Update();
    myc[4]->Print((outfile+".pdf").c_str());
    
    myc[5]->Update();
    myc[5]->Print((outfile+".pdf").c_str());
    
    myc[6]->Update();
    myc[6]->Print((outfile+".pdf").c_str());
    
    slopecor[ieta] = 1;
    offsetcor[ieta] = 0;
    
    if (calibAndReso(t,myc,7,evec,cut1.str(),cut2.str(),"eSR4Reco",offsetcor[ieta],slopecor[ieta],false)) return 1;
    
    myc[7]->Update();
    myc[7]->Print((outfile+".pdf").c_str());
    
    myc[8]->Update();
    myc[8]->Print((outfile+".pdf").c_str());
    
    myc[9]->Update();
    myc[9]->Print((outfile+".pdf").c_str());
    
    std::cout << " -- Eta cor: offset and slope = " << offsetcor[ieta] << " " << slopecor[ieta] << std::endl;

    offsettmp = offsetcor[ieta];
    slopetmp = slopecor[ieta];

    if (calibAndReso(t,myc,10,evec,cut1.str(),cut2.str(),"eSR4Reco",offsettmp,slopetmp,true)) return 1;
    
    std::cout << " -- Eta cor after calibration: offset and slope = " << offsettmp << " " << slopetmp << std::endl;
    
    myc[10]->Update();
    myc[10]->Print((outfile+".pdf").c_str());
    
    myc[11]->Update();
    myc[11]->Print((outfile+".pdf").c_str());
    
    myc[12]->Update();
    myc[12]->Print((outfile+".pdf").c_str());
    
  }//loop on eta

  myc[13]->cd();
  TGraphErrors *grSlope = new TGraphErrors(neta,eta,slope);
  grSlope->SetMarkerStyle(22);
  grSlope->SetMarkerColor(2);
  grSlope->SetLineColor(2);
  grSlope->Draw("APE");
  grSlope->GetXaxis()->SetTitle("|#eta|");
  grSlope->GetYaxis()->SetTitle("slope (Mips/GeV)");
  TGraphErrors *grSlopeCor = new TGraphErrors(neta,eta,slopecor);
  grSlopeCor->SetMarkerStyle(23);
  grSlopeCor->SetMarkerColor(4);
  grSlopeCor->SetLineColor(4);
  grSlopeCor->Draw("PEsame");

  myc[13]->Update();
  myc[13]->Print((outfile+".pdf").c_str());


  myc[14]->cd();
  TGraphErrors *grOffset = new TGraphErrors(neta,eta,offset);
  grOffset->SetMarkerStyle(22);
  grOffset->SetMarkerColor(2);
  grOffset->SetLineColor(2);
  grOffset->Draw("APE");
  grOffset->GetXaxis()->SetTitle("|#eta|");
  grOffset->GetYaxis()->SetTitle("offset (Mips/GeV)");
  TGraphErrors *grOffsetCor = new TGraphErrors(neta,eta,offsetcor);
  grOffsetCor->SetMarkerStyle(23);
  grOffsetCor->SetMarkerColor(4);
  grOffsetCor->SetLineColor(4);
  grOffsetCor->Draw("PEsame");

  myc[14]->Update();
  myc[14]->Print((outfile+".pdf").c_str());



  myc[0]->Print((outfile+".pdf]").c_str());
  return 0;
  
}
