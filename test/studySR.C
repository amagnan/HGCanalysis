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

#include "utilities.C"

double calibratedE(const double Etot, const double eta,
		   const unsigned cells=0){
  //calibration for signal region 2: 3*3 cm^2
  double pars[3] = {69.5,4.5,-0.8};
  //double pars[3] = {75.5,0,0};
  double paro[3] = {-34.4,0,0};
  //double pars[3] = {77,3.4,-0.50};
  //double paro[3] = {-11.6,-7.7,-8.8};
  //double paro[3] = {-5.3,-12.8,-6.9};  
  double offset = paro[0] + paro[1]*fabs(eta) + paro[2]*eta*eta;
  double slope = pars[0] + pars[1]*fabs(eta) + pars[2]*eta*eta;
  if (cells==0) return (Etot-offset)/slope;
  if (cells==1) {
    offset += slope*-4.4;//in original unit=mips...
    slope *= 0.549;
  }
  if (cells==3) {
    offset += slope*-2.09;//in original unit=mips...
    slope *= 1.009;
  }
  if (cells==35) {
    //offset += slope*-1.8;//in original unit=mips...
    slope *= 1.032;
  }
  if (cells==37) {
    //offset += slope*-1.5;//in original unit=mips...
    slope *= 1.035;
  }
  if (cells==39) {
    //offset += slope*-1.4;//in original unit=mips...
    slope *= 1.037;
  }
  if (cells==5) slope *= 1.10;
  if (cells==57) slope *= 1.116;
  if (cells==59) slope *= 1.120;
  if (cells==7) slope *= 1.143;
  if (cells==79) slope *= 1.149;
  return (Etot-offset)/slope;
};


int studySR() {

  bool doPu = false;
  setTDRStyle();

  std::string dir;
  if (doPu) dir="PLOTS_SR_pu140";
  else dir="PLOTS_SR_pu0";

  //true for calib
  std::string cut1 = "isValid1==1 && converted1==0 && eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetid1==0 && invalidNeighbour1==0";
  std::string cut2 = "isValid2==1 && converted2==0 && eTrue2/cosh(etaTrue2)>40 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetid2==0 && invalidNeighbour2==0";


  std::string cut1pca = "eTrue1/cosh(etaTrue1)>40 && TMath::Abs(etaTrue1)>1.6 && TMath::Abs(etaTrue1)<2.8 && (TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue1)/(2*TMath::Pi())*360.)%20>11) && invalidDetidPCA1==0 && dRTrue1 < 0.05 && invalidNeighbourPCA1==0";// && noPhiCrackPCA1==1";
  std::string cut2pca = "eTrue2/cosh(etaTrue2)>40 && TMath::Abs(etaTrue2)>1.6 && TMath::Abs(etaTrue2)<2.8 && (TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(phiTrue2)/(2*TMath::Pi())*360.)%20>11) && invalidDetidPCA2==0 && dRTrue2 < 0.05 && invalidNeighbourPCA2==0";// && noPhiCrackPCA2==1";

  TFile *fin;
  if (doPu) fin = TFile::Open("PCA_topoFix_Hgg_140pu.root");
  else fin = TFile::Open("PCA_topoFix_Hgg_0pu.root");
  //if (doPu) fin = TFile::Open("PCA_Hgg_140pu.root");
  //else fin = TFile::Open("PCA_Hgg_0pu.root");
  if (!fin) return 1;
  fin->cd("hgg");

  TTree *t = (TTree*)gDirectory->Get("tree");

  if (!t) {
    std::cout << " Error, tree not found!" << std::endl;
    return 1;
  }

  const unsigned nC = 10;
  TCanvas *myc[nC];
  for (unsigned ic(0); ic<nC; ++ic){
    std::ostringstream label;
    label << "myc" << ic ;
    myc[ic] = new TCanvas(label.str().c_str(),
			  label.str().c_str(),
			  1000,600);
  }

  TLatex lat;
  TLatex latSmall;
  latSmall.SetTextSize(0.035);
  char buf[100];

  unsigned SR[nC] = {1,3,5,7,35,37,39,57,59,79};

  std::ostringstream printStr;

  for (unsigned ic(0);ic<nC;++ic){
    myc[ic]->cd();

    plotEnergyCalib(t,myc[ic],SR[ic],true,cut1,cut2);
    sprintf(buf,"Region %d",SR[ic]);
    lat.DrawLatex(50,600,buf);
    //if (doPu) latSmall.DrawLatex(0.2,1.07*tot2->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU140");
    //else latSmall.DrawLatex(0.2,1.07*tot2->GetMaximum(),"6_2_0_SLHC25_p2 Relval Hgg PU0");
    myc[ic]->Update();
    printStr.str("");
    printStr << dir << "/energyCalibration_SR" << SR[ic] << "_trueDir_trueUnconv.pdf";
    myc[ic]->Print(printStr.str().c_str());
  }

  unsigned proceed = 0;
  std::cout << " -- Continue ? " << std::endl;
  std::cin >> proceed;
  if (!proceed) return 0;

  for (unsigned ic(0);ic<nC;++ic){
    myc[ic]->cd();
    plotEnergyCalibvsEta(t,myc[ic],SR[ic],false,cut1pca+" && converted1==0",cut2pca+" && converted2==0");
    sprintf(buf,"Region %d",SR[ic]);
    lat.DrawLatex(-2.9,2,buf);
    myc[ic]->Update();
    printStr.str("");
    printStr << dir << "/energyPCACalibrationvseta_SR" << SR[ic] << "_pcaDir_trueUnconv.pdf";
    myc[ic]->Print(printStr.str().c_str());
  }

  std::cout << " -- Continue ? " << std::endl;
  std::cin >> proceed;
  if (!proceed) return 0;

  for (unsigned ic(0);ic<nC;++ic){
    myc[ic]->cd();
    plotEnergyCalibvsPhi(t,myc[ic],SR[ic],false,cut1pca+" && converted1==0",cut2pca+" && converted2==0");
    sprintf(buf,"Region %d",SR[ic]);
    lat.DrawLatex(-2.9,2,buf);
    myc[ic]->Update();
    printStr.str("");
    printStr << dir << "/energyPCACalibrationvsphi_SR" << SR[ic] << "_pcaDir_trueUnconv.pdf";
    myc[ic]->Print(printStr.str().c_str());
  }

  std::cout << " -- Continue ? " << std::endl;
  std::cin >> proceed;
  if (!proceed) return 0;


  plotConversionVeto(t,myc[0],myc[1],cut1pca,cut2pca,
		     "e_{3#times3}/e_{5#times5}",
		     "calibratedE(eSR3RecoPCA1,etaPCA1,3)/calibratedE(eSR5RecoPCA1,etaPCA1,5)",
		     "calibratedE(eSR3RecoPCA2,etaPCA2,3)/calibratedE(eSR5RecoPCA2,etaPCA2,5)",
		     50,0,1.1);

  myc[0]->Update();
  printStr.str("");
  printStr << dir << "/e3x3overe5x5_pcaDir.pdf";
  myc[0]->Print(printStr.str().c_str());

  myc[1]->Update();
  printStr.str("");
  printStr << dir << "/e3x3overe5x5integrated_pcaDir.pdf";
  myc[1]->Print(printStr.str().c_str());


  double paruc[3];

  plotCalibVsFirstLayerFraction(t,myc[2],0,
				cut1,
				cut2,
				"e_{3#times3}(L1-5)/e_{3#times3}",
				"calibratedE(eSR3Reco1,etaTrue1,3)/eTrue1:calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaTrue1,3)/calibratedE(eSR3Reco1,etaTrue1,3)",
				"calibratedE(eSR3Reco2,etaTrue2,3)/eTrue2:calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaTrue2,3)/calibratedE(eSR3Reco2,etaTrue2,3)",
				50,0,0.1,50,0,1.2,
				0.9,1.2,
				paruc[0],paruc[1],paruc[2]
				);

  myc[2]->Update();
  TF1 *fitnopu = new TF1("fitnopu","1.017-1.29*x+8.5*x*x",0,0.1);
  fitnopu->SetLineColor(7);
  fitnopu->SetLineWidth(2);
  fitnopu->Draw("same");
  printStr.str("");
  printStr << dir << "/eReco33overeTruevse33L1-5oe33_trueUnconv_trueDir.pdf";
  myc[2]->Print(printStr.str().c_str());

  plotCalibVsFirstLayerFraction(t,myc[3],10,
				cut1pca+" && calibratedE(eSR3RecoPCA1,etaPCA1,3)/calibratedE(eSR5RecoPCA1,etaPCA1,5)>0.99",
				cut2pca+" && calibratedE(eSR3RecoPCA2,etaPCA2,3)/calibratedE(eSR5RecoPCA2,etaPCA2,5)>0.99",
				"e_{3#times3}(L1-5)/e_{3#times3}",
				"calibratedE(eSR3RecoPCA1,etaPCA1,3)/eTrue1:calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaTrue1,3)/calibratedE(eSR3Reco1,etaTrue1,3)",
				"calibratedE(eSR3RecoPCA2,etaPCA2,3)/eTrue2:calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaTrue2,3)/calibratedE(eSR3Reco2,etaTrue2,3)",
				50,0,0.1,50,0,1.2,
				0.9,1.2,
				paruc[0],paruc[1],paruc[2]
				);

  myc[3]->Update();
  fitnopu->Draw("same");
  printStr.str("");
  printStr << dir << "/eReco33overeTruevse33L1-5oe33_unconv5x5_pcaDir.pdf";
  myc[3]->Print(printStr.str().c_str());

  std::cout << " -- Continue ? " << std::endl;
  std::cin >> proceed;
  if (!proceed) return 0;

  double parc[nC][3];

  for (unsigned ic(0);ic<nC;++ic){
    myc[ic]->cd();
    std::ostringstream var1,var2;
    var1 << "calibratedE(eSR" << SR[ic] << "RecoPCA1,etaPCA1," << SR[ic] << ")/eTrue1";
    var2 << "calibratedE(eSR" << SR[ic] << "RecoPCA2,etaPCA2," << SR[ic] << ")/eTrue2";
    plotCalibVsFirstLayerFraction(t,myc[ic],SR[ic],
				  cut1pca+" && calibratedE(eSR3RecoPCA1,etaPCA1,3)/calibratedE(eSR5RecoPCA1,etaPCA1,5)<=0.99",
				  cut2pca+" && calibratedE(eSR3RecoPCA2,etaPCA2,3)/calibratedE(eSR5RecoPCA2,etaPCA2,5)<=0.99",
				  "e_{3#times3}(L1-5)/e_{3#times3}",
				  (var1.str()+":calibratedE(eSR3RecoLayer0_1+eSR3RecoLayer1_1+eSR3RecoLayer2_1+eSR3RecoLayer3_1+eSR3RecoLayer4_1,etaTrue1,3)/calibratedE(eSR3Reco1,etaTrue1,3)").c_str(),
				  (var2.str()+":calibratedE(eSR3RecoLayer0_2+eSR3RecoLayer1_2+eSR3RecoLayer2_2+eSR3RecoLayer3_2+eSR3RecoLayer4_2,etaTrue2,3)/calibratedE(eSR3Reco2,etaTrue2,3)").c_str(),
				  50,0,0.1,50,0,1.2,
				  0.9,1.2,
				  parc[ic][0],parc[ic][1],parc[ic][2]
				  );

    myc[ic]->Update();
    fitnopu->Draw("same");
    sprintf(buf,"Region %d",SR[ic]);
    lat.DrawLatex(0,1.21,buf);
    printStr.str("");
    printStr << dir << "/eReco" << SR[ic] << "overeTruevse33L1-5oe33_conv5x5_pcaDir.pdf";
    myc[ic]->Print(printStr.str().c_str());
  }
  std::cout << " -- Continue ? " << std::endl;
  std::cin >> proceed;
  if (!proceed) return 0;

  for (unsigned ic(0);ic<nC;++ic){
    myc[ic]->cd();
    plotEnergyResolution(t,myc[ic],SR[ic],
			 cut1pca,cut2pca,
			 //1.017,-1.29,8.5,
			 paruc[0],paruc[1],paruc[2],
			 parc[ic][0],parc[ic][1],parc[ic][2]
			 );
    //1.02141,-0.861073,3.39355,
    //0.993,-0.60,0.);//pu140:1.022,-0.8,0.
    myc[ic]->Update();
    sprintf(buf,"Region %d",SR[ic]);
    lat.DrawLatex(0,20,buf);
    printStr.str("");
    //printStr << dir << "/eReco" << SR[ic] << "overeTrue_conv5x5_cor1-5_pcaDir.pdf";
    printStr << dir << "/eReco" << SR[ic] << "overeTrue_conv5x5_pcaDir.pdf";
    myc[ic]->Print(printStr.str().c_str());
  }
 

  return 0;
}
