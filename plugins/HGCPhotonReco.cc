// -*- C++ -*-
//
// Package:    HGCPhotonReco
// Class:      HGCPhotonReco
// 
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_GenEvent.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "TTree.h"
#include "TMath.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include "TH1.h"
#include "TH2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"

#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "TFile.h"
#include "TGraphErrors.h"

  
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"

#include "PCAShowerAnalysis.h"
using namespace std;



//
// class declaration
//

// information to be loaded into TTree
struct info_t {
  bool isValid;
  bool invalidDetid;
  bool invalidNeighbour;
  bool invalidDetidPCA;
  bool invalidNeighbourPCA;
  float etaTrue;
  float phiTrue;
  float etadetTrue;
  float dRTrue;
  float eTrue;
  float xvtxTrue;
  float yvtxTrue;
  float zvtxTrue;
  float etaSC;
  float phiSC;
  float etaSeed;
  float phiSeed;
  float etaReco;
  float phiReco;
  float etadetReco;
  float phiWidth;
  float etaWidth;
  int clustersSize;
  int nClusters09;
  float eSeed;
  float eSC;
  float eReco[3];
  float eRecoPCA[3];
  int converted;
  unsigned showerMax;
  double zShowerMax;
  bool noPhiCrack;
  bool noPhiCrackPCA;
  float xPCA;
  float yPCA;
  float zPCA;
  float etaPCA;
  float phiPCA;

  void correctForEtaDep(const ROOT::Math::XYZVector & posSC, const bool isPCA=false){
    if (!isPCA) {
      ROOT::Math::XYZVector vtx = ROOT::Math::XYZVector(xvtxTrue,yvtxTrue,zvtxTrue);
      etadetReco = posSC.eta();
      ROOT::Math::XYZVector photon = posSC-vtx;
      etaReco = photon.eta();
      phiReco = photon.phi();
      if (fabs(etaReco)>5) return;
    }
  };

  void fillSC(const edm::Ptr<reco::SuperCluster> & aSC, PCAShowerAnalysis  & pcaShowerAnalysis){
    if (aSC.isNull()) return;
    etaSC = aSC->eta();
    phiSC = aSC->phi();
    eSC = aSC->energy();
    etaSeed  =aSC->seed()->eta();
    phiSeed  =aSC->seed()->phi();
    eSeed  =aSC->seed()->energy();
    clustersSize = static_cast<int>(aSC->clustersSize());
    etaWidth=aSC->etaWidth();
    phiWidth=aSC->phiWidth();
    double etot = 0;
    //std::cout << " --- Found " << aSC->clusters().size() << " clusters in SC" << std::endl;
    for (unsigned int ic =0 ; ic < aSC->clusters().size() ; ic++){
      etot += aSC->clusters()[ic]->energy();
      //std::cout << " -- clus " << ic << " E=" << aSC->clusters()[ic]->energy() << " etot = " << etot << std::endl;
      if(etot/eSC <=0.9) nClusters09++;
      //if(etot/eSC >0.9) break;
    }

    //get shower position and direction
    GlobalPoint pcaShowerPos;
    GlobalVector pcaShowerDir;
    pcaShowerAnalysis.showerParameters(&(*aSC),pcaShowerPos,pcaShowerDir);
    xPCA = pcaShowerPos.x();
    yPCA = pcaShowerPos.y();
    zPCA = pcaShowerPos.z();
    etaPCA =  pcaShowerDir.eta();
    phiPCA =  pcaShowerDir.phi();

  };

  void fillTruth(const edm::Ptr<reco::GenParticle> & aPhoton, const edm::Handle<edm::SimTrackContainer> SimTk, const edm::Handle<edm::SimVertexContainer> SimVtx){
    if (aPhoton.isNull()) return;
    etaTrue = aPhoton->eta();
    phiTrue = aPhoton->phi();
    eTrue = aPhoton->energy();
    xvtxTrue = (aPhoton->vertex()).X();
    yvtxTrue = (aPhoton->vertex()).Y();
    zvtxTrue = (aPhoton->vertex()).Z();

    math::XYZVectorD hitPos=getInteractionPositionLC(SimTk.product(),SimVtx.product(), aPhoton->pt()).pos;
    const double z = std::fabs(hitPos.z());
    converted = (unsigned)(z < 317 && z > 1e-3);
  };

  void initialise(){
    isValid=false;
    invalidDetid=false;
    invalidNeighbour=false;
    invalidDetidPCA=false;
    invalidNeighbourPCA=false;
    etaTrue=-5.;
    etadetTrue=-5.;
    phiTrue=-5.;
    dRTrue=-5;
    eTrue=-1.;
    xvtxTrue = 0;
    yvtxTrue = 0;
    zvtxTrue = 0;
    etaSC=-5.;
    phiSC=-5.;
    etaSeed=-5.;
    phiSeed=-5.;
    etaReco = -5;
    etadetReco = -5;
    phiReco = -5;
    for (unsigned i(0);i<3;++i){
      eReco[i]  = 0;         
      eRecoPCA[i]  = 0;         
    }
    eSeed     =-1.;
    eSC = -1;
    etaWidth=-1.;
    phiWidth=-1.;
    clustersSize=-1;
    nClusters09=-1;
    converted=-1;
    showerMax = 0;
    zShowerMax = 0;
    noPhiCrack=true;
    noPhiCrackPCA=true;
    xPCA = 0;
    yPCA = 0;
    zPCA = 0;
    etaPCA = 0;
    phiPCA = 0;
  };
};

// .h class info
class HGCPhotonReco : public edm::EDAnalyzer {
public:
  explicit HGCPhotonReco(const edm::ParameterSet&);
  ~HGCPhotonReco();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  double DeltaPhi(const double & phi1, const double & phi2);

  void fillzPositions(const edm::PtrVector<HGCRecHit>& rechitvec);

  void getPhotonEnergy(const edm::PtrVector<HGCRecHit>& rechitvec,info_t & info);
  void getPhotonEnergy(const edm::PtrVector<HGCRecHit>& rechitvec,info_t & info,const double & phimax, const double & etamax, const bool isPCA=false);
  void getMaximumCell(const edm::PtrVector<HGCRecHit>& rechitvec,const double & phimax,const double & etamax,std::vector<HGCEEDetId> & detidmax);
  //void getTotalEnergy(const edm::PtrVector<HGCRecHit>& rechitvec,
  //const std::vector<HGCEEDetId> & eventPos,
  //info_t & info);
  /*void getEnergyWeightedPosition(const edm::PtrVector<HGCRecHit>& rechitvec,
				 const std::vector<HGCEEDetId> & xmax,
				 std::vector<HGCEEDetId> & recoPos,
				 std::vector<double> & recoE);*/

  float resumEmEnergyTest( const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters);
   bool isInFid(const edm::Ptr<reco::GenParticle> & aPhoton);
  double getW0(const unsigned layer);
  double absWeight(const unsigned layer, const bool dedx=false);
  //double calibratedE(const double Etot, const double eta);
  void fillNeighbours(const HGCEEDetId & detidmax,
		      std::vector<double> & Exy,
		      info_t & info,
		      bool isPCA);
private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


  edm::EDGetTokenT<edm::View<HGCRecHit> >endcapRecHitCollection_      ; 
  edm::EDGetTokenT<edm::View<reco::SuperCluster> >endcapSuperClusterCollection_;
  edm::EDGetTokenT<edm::View<reco::PFCluster> >endcapClusterCollection_     ;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >genParticlesCollection_     ;

  std::string g4TracksSource_, g4VerticesSource_; 
  std::vector<std::string> geometrySource_;
  const HGCalGeometry * hgcEEGeom_;

  edm::Handle<HGCRecHitCollection> recHits_;	
  TTree *tree_;
  info_t info1_;
  info_t info2_;
  
  const TGraphErrors * _hgcOverburdenParam;
  const TGraphErrors * _hgcLambdaOverburdenParam;
  const std::vector<double> _weights_ee;

  unsigned nNoClusters_;
  unsigned nLayers_;
  double cellSize_;
  bool doLogWeight_;
  double mipE_;
  unsigned nSR_;
  unsigned debug_;
  bool singleGamma_;
  std::map<unsigned,double> zPos_;
  ROOT::Math::XYZPoint truthVtx_;

  //histograms
  TH2F *hEvsLayer_;
  //TH2F *hetavsphi_[30];
  TH2F *hyvsx_[30];
  TH2F *hxvsz_;
  TH2F *hyvsz_;
  TH2F *hrvsz_;
  //TH2F *hphisecvscellid_[30];
  // TH2F *hyvsxzoom_[30];

  TH1F *dxMaxTruth_;
  TH1F *dyMaxTruth_;
  TH1F *dzMaxTruth_;
  TH2F *dphiMaxTruthvsLayer_;
  TH1F *dxMaxNeigh_[9];
  TH1F *dyMaxNeigh_[9];
  TH1F *dzMaxNeigh_[9];
  TH1F *dRmin_[30];


};

// constructor
HGCPhotonReco::HGCPhotonReco(const edm::ParameterSet& iConfig):
  endcapRecHitCollection_(consumes <edm::View<HGCRecHit> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapRecHitCollection",     edm::InputTag("HGCalRecHit:HGCEERecHits")))),
  endcapSuperClusterCollection_(consumes <edm::View<reco::SuperCluster> >(iConfig.getUntrackedParameter<edm::InputTag>("endcapSuperClusterCollection",edm::InputTag("particleFlowSuperClusterHGCEE")))),
  endcapClusterCollection_(consumes <edm::View<reco::PFCluster> > (iConfig.getUntrackedParameter<edm::InputTag>("endcapClusterCollection",edm::InputTag("particleFlowClusterHGCEE")))),
  genParticlesCollection_(consumes <edm::View<reco::GenParticle> > (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag",edm::InputTag("genParticles")))),
  _hgcOverburdenParam(nullptr),
  _hgcLambdaOverburdenParam(nullptr),
  _weights_ee(iConfig.getParameter<std::vector<double> >("weights_ee")),
  debug_(iConfig.getParameter<unsigned>("debug")),
  singleGamma_(iConfig.getParameter<bool>("singleGamma"))
{ 
  g4TracksSource_           = iConfig.getUntrackedParameter<std::string>("g4TracksSource");
  g4VerticesSource_         = iConfig.getUntrackedParameter<std::string>("g4VerticesSource");
	geometrySource_ = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
	edm::Service<TFileService> fs_;

	tree_ = fs_->make<TTree>("tree","");
	tree_->Branch("isValid1"              ,&info1_.isValid             ,"isValid1/B");
	tree_->Branch("invalidDetid1"              ,&info1_.invalidDetid             ,"invalidDetid1/B");
	tree_->Branch("invalidNeighbour1"              ,&info1_.invalidNeighbour             ,"invalidNeighbour1/B");
	tree_->Branch("invalidDetidPCA1"              ,&info1_.invalidDetidPCA             ,"invalidDetidPCA1/B");
	tree_->Branch("invalidNeighbourPCA1"              ,&info1_.invalidNeighbourPCA             ,"invalidNeighbourPCA1/B");
	tree_->Branch("etaTrue1"              ,&info1_.etaTrue             ,"etaTrue1/F");
	tree_->Branch("etadetTrue1"              ,&info1_.etadetTrue             ,"etadetTrue1/F");
	tree_->Branch("phiTrue1"              ,&info1_.phiTrue             ,"phiTrue1/F");
	tree_->Branch("dRTrue1"              ,&info1_.dRTrue            ,"dRTrue1/F");
	tree_->Branch("eTrue1"              ,&info1_.eTrue            ,"eTrue1/F");
	tree_->Branch("xvtxTrue1"              ,&info1_.xvtxTrue            ,"xvtxTrue1/F");
	tree_->Branch("yvtxTrue1"              ,&info1_.yvtxTrue            ,"yvtxTrue1/F");
	tree_->Branch("zvtxTrue1"              ,&info1_.zvtxTrue            ,"zvtxTrue1/F");

	tree_->Branch("etaSC1"              ,&info1_.etaSC             ,"etaSC1/F");
	tree_->Branch("phiSC1"              ,&info1_.phiSC             ,"phiSC1/F");
	tree_->Branch("etaSeed1"              ,&info1_.etaSeed             ,"etaSeed1/F");
	tree_->Branch("phiSeed1"              ,&info1_.phiSeed             ,"phiSeed1/F");
	tree_->Branch("etaReco1"              ,&info1_.etaReco             ,"etaReco1/F");
	tree_->Branch("etadetReco1"              ,&info1_.etadetReco             ,"etadetReco1/F");
	tree_->Branch("phiReco1"              ,&info1_.phiReco             ,"phiReco1/F");
	tree_->Branch("etaWidth1"              ,&info1_.etaWidth             ,"etaWidth1/F");
	tree_->Branch("phiWidth1"              ,&info1_.phiWidth             ,"phiWidth1/F");
	tree_->Branch("clustersSize1"              ,&info1_.clustersSize           ,"clustersSize1/I");
	tree_->Branch("nClusters091"              ,&info1_.nClusters09           ,"nClusters091/I");


	tree_->Branch("eSeed1"              ,&info1_.eSeed            ,"eSeed1/F");
	tree_->Branch("eSC1"              ,&info1_.eSC            ,"eSC1/F");

	tree_->Branch("eSR3Reco1"              ,&info1_.eReco[0]            ,"eSR3Reco1/F");
	tree_->Branch("eSR5Reco1"              ,&info1_.eReco[1]            ,"eSR5Reco1/F");
	tree_->Branch("eSR7Reco1"              ,&info1_.eReco[2]            ,"eSR7Reco1/F");
	tree_->Branch("eSR3RecoPCA1"              ,&info1_.eRecoPCA[0]            ,"eSR3RecoPCA1/F");
	tree_->Branch("eSR5RecoPCA1"              ,&info1_.eRecoPCA[1]            ,"eSR5RecoPCA1/F");
	tree_->Branch("eSR7RecoPCA1"              ,&info1_.eRecoPCA[2]            ,"eSR7RecoPCA1/F");
	tree_->Branch("converted1"              ,&info1_.converted            ,"converted1/I");
	tree_->Branch("showerMax1"              ,&info1_.showerMax            ,"showerMax1/I");
	tree_->Branch("noPhiCrack1"              ,&info1_.noPhiCrack             ,"noPhiCrack1/B");
	tree_->Branch("noPhiCrackPCA1"              ,&info1_.noPhiCrackPCA             ,"noPhiCrackPCA1/B");

	tree_->Branch("xPCA1"              ,&info1_.xPCA            ,"xPCA1/F");
	tree_->Branch("yPCA1"              ,&info1_.yPCA            ,"yPCA1/F");
	tree_->Branch("zPCA1"              ,&info1_.zPCA            ,"zPCA1/F");
	tree_->Branch("etaPCA1"              ,&info1_.etaPCA            ,"etaPCA1/F");
	tree_->Branch("phiPCA1"              ,&info1_.phiPCA            ,"phiPCA1/F");


	if (!singleGamma_){
	  //photon 2
	  tree_->Branch("isValid2"              ,&info2_.isValid             ,"isValid2/B");
	tree_->Branch("invalidDetid2"              ,&info2_.invalidDetid             ,"invalidDetid2/B");
	tree_->Branch("invalidNeighbour2"              ,&info2_.invalidNeighbour             ,"invalidNeighbour2/B");
	tree_->Branch("invalidDetidPCA2"              ,&info2_.invalidDetidPCA             ,"invalidDetidPCA2/B");
	tree_->Branch("invalidNeighbourPCA2"              ,&info2_.invalidNeighbourPCA             ,"invalidNeighbourPCA2/B");
	  tree_->Branch("etaTrue2"              ,&info2_.etaTrue             ,"etaTrue2/F");
	  tree_->Branch("etadetTrue2"              ,&info2_.etadetTrue             ,"etadetTrue2/F");
	  tree_->Branch("phiTrue2"              ,&info2_.phiTrue             ,"phiTrue2/F");
	  tree_->Branch("dRTrue2"              ,&info2_.dRTrue            ,"dRTrue2/F");
	  tree_->Branch("eTrue2"              ,&info2_.eTrue            ,"eTrue2/F");
	  tree_->Branch("xvtxTrue2"              ,&info2_.xvtxTrue            ,"xvtxTrue2/F");
	  tree_->Branch("yvtxTrue2"              ,&info2_.yvtxTrue            ,"yvtxTrue2/F");
	  tree_->Branch("zvtxTrue2"              ,&info2_.zvtxTrue            ,"zvtxTrue2/F");
	  
	  tree_->Branch("etaSC2"              ,&info2_.etaSC             ,"etaSC2/F");
	  tree_->Branch("phiSC2"              ,&info2_.phiSC             ,"phiSC2/F");
	  tree_->Branch("etaReco2"              ,&info2_.etaReco             ,"etaReco2/F");
	  tree_->Branch("etadetReco2"              ,&info2_.etadetReco             ,"etadetReco2/F");
	  tree_->Branch("phiReco2"              ,&info2_.phiReco             ,"phiReco2/F");
	  tree_->Branch("etaSeed2"              ,&info2_.etaSeed             ,"etaSeed2/F");
	  tree_->Branch("phiSeed2"              ,&info2_.phiSeed             ,"phiSeed2/F");
	  tree_->Branch("etaWidth2"              ,&info2_.etaWidth             ,"etaWidth2/F");
	  tree_->Branch("phiWidth2"              ,&info2_.phiWidth             ,"phiWidth2/F");
	  tree_->Branch("clustersSize2"              ,&info2_.clustersSize           ,"clustersSize2/I");
	  tree_->Branch("nClusters092"              ,&info2_.nClusters09           ,"nClusters092/I");
	  
	  
	  tree_->Branch("eSeed2"              ,&info2_.eSeed            ,"eSeed2/F");
	  tree_->Branch("eSC2"              ,&info2_.eSC            ,"eSC2/F");
	  tree_->Branch("eSR3Reco2"              ,&info2_.eReco[0]            ,"eSR3Reco2/F");
	  tree_->Branch("eSR5Reco2"              ,&info2_.eReco[1]            ,"eSR5Reco2/F");
	  tree_->Branch("eSR7Reco2"              ,&info2_.eReco[2]            ,"eSR7Reco2/F");
	  tree_->Branch("eSR3RecoPCA2"              ,&info2_.eRecoPCA[0]            ,"eSR3RecoPCA2/F");
	  tree_->Branch("eSR5RecoPCA2"              ,&info2_.eRecoPCA[1]            ,"eSR5RecoPCA2/F");
	  tree_->Branch("eSR7RecoPCA2"              ,&info2_.eRecoPCA[2]            ,"eSR7RecoPCA2/F");

	  tree_->Branch("converted2"              ,&info2_.converted            ,"converted2/I");
	  tree_->Branch("showerMax2"              ,&info2_.showerMax            ,"showerMax2/I");
	  tree_->Branch("noPhiCrack2"              ,&info2_.noPhiCrack             ,"noPhiCrack2/B");
	  tree_->Branch("noPhiCrackPCA2"              ,&info2_.noPhiCrackPCA             ,"noPhiCrackPCA2/B");

	tree_->Branch("xPCA2"              ,&info2_.xPCA            ,"xPCA2/F");
	tree_->Branch("yPCA2"              ,&info2_.yPCA            ,"yPCA2/F");
	tree_->Branch("zPCA2"              ,&info2_.zPCA            ,"zPCA2/F");
	tree_->Branch("etaPCA2"              ,&info2_.etaPCA            ,"etaPCA2/F");
	tree_->Branch("phiPCA2"              ,&info2_.phiPCA            ,"phiPCA2/F");
	}
	if(iConfig.exists("hgcOverburdenParamFile"))
	  {
	    edm::FileInPath fp = iConfig.getParameter<edm::FileInPath>("hgcOverburdenParamFile");
	    TFile *fIn=TFile::Open(fp.fullPath().c_str());
	    if(fIn)
	      {
		_hgcOverburdenParam=(const TGraphErrors *) fIn->Get("x0Overburden");
		_hgcLambdaOverburdenParam = (const TGraphErrors *) fIn->Get("lambdaOverburden");
		fIn->Close();
	      }
	  }

	//@TODO
	//get this from geometry !!
	nLayers_ = 30;
	cellSize_ = 1;
	doLogWeight_ = true;
	mipE_ = 0.0000551;
	nSR_ = 5;
	truthVtx_ = ROOT::Math::XYZPoint(0,0,0);
	hEvsLayer_ = fs_->make<TH2F>("hEvsLayer_",";layer;E;photons",
				     nLayers_,0,nLayers_,
				     1000,0,10000);

	dxMaxTruth_ = fs_->make<TH1F>("dxMaxTruth_",";#Deltax(truth,cell) (mm);cells",100,-50,50);
	dyMaxTruth_ = fs_->make<TH1F>("dyMaxTruth_",";#Deltay(truth,cell) (mm);cells",100,-50,50);
	dzMaxTruth_ = fs_->make<TH1F>("dzMaxTruth_",";#Deltaxz(truth,cell) (mm);cells",100,-1,1);
	dphiMaxTruthvsLayer_ = fs_->make<TH2F>("dphiMaxTruthvsLayer_",";layer;#Delta#phi(truth,cell);cells",nLayers_,0,nLayers_,100,-1,1);

	for (unsigned idx(0);idx<9;++idx){
	  std::ostringstream label;
	  label << "dxMaxNeigh_" << idx;
	  dxMaxNeigh_[idx] = fs_->make<TH1F>(label.str().c_str(),";#Deltax(neigh,cell) (mm);cells",100,-50,50);
	  label.str("");
	  label << "dyMaxNeigh_" << idx;
	  dyMaxNeigh_[idx] = fs_->make<TH1F>(label.str().c_str(),";#Deltay(neigh,cell) (mm);cells",100,-50,50);
	  label.str("");
	  label << "dzMaxNeigh_" << idx;
	  dzMaxNeigh_[idx] = fs_->make<TH1F>(label.str().c_str(),";#Deltaxz(neigh,cell) (mm);cells",100,-1,1);
	}

	for (unsigned iL(0);iL<nLayers_;++iL){
	  std::ostringstream label;
	  label << "dRmin_" << iL;
	  dRmin_[iL] = fs_->make<TH1F>(label.str().c_str(),";dRmin;showers",
				       100,0,0.5);
	}
	
	for (unsigned iL(0);iL<nLayers_;++iL){
	  std::ostringstream label;
	  //label << "hetavsphi_" << iL;
	  //hetavsphi_[iL] = fs_->make<TH2F>(label.str().c_str(),
	  //";#phi;#eta;hits",
	  //360,-3.1416,3.1416,
	  //30,1.5,3.0);
	  label.str("");
	  label << "hyvsx_" << iL;
	  hyvsx_[iL] = fs_->make<TH2F>(label.str().c_str(),
				       ";x (cm);y (cm); hits",
				       340,-170,170,
				       340,-170,170);
	  //label.str("");
	  //label << "hphisecvscellid_" << iL;
	  //hphisecvscellid_[iL] = fs_->make<TH2F>(label.str().c_str(),
	  //";cellid;phi sector; hits",
	  //2120,0,2120,
	  //20,0,20);
	  //label.str("");
	  //label << "hyvsxzoom_" << iL;
	  //hyvsxzoom_[iL] = fs_->make<TH2F>(label.str().c_str(),
	  //";x (cm);y (cm); hits",
	  //1200,-60,60,
	  //1200,-60,60);
	}
	hxvsz_ = fs_->make<TH2F>("hxvsz",
				 ";z (cm);x (cm); hits",
				 600,-30,30,
				 340,-170,170);
	hyvsz_ = fs_->make<TH2F>("hyvsz",
				 ";z (cm);y (cm); hits",
				 600,-30,30,
				 340,-170,170);
	hrvsz_ = fs_->make<TH2F>("hrvsz",
				 ";z (cm);r (cm); hits",
				 600,-30,30,
				 170,0,170);

}

// destructor
HGCPhotonReco::~HGCPhotonReco()
{

}


//
// member functions
//

double HGCPhotonReco::DeltaPhi(const double & phi1, const double & phi2){
  double dphi = phi1 - phi2;
  if (dphi< (-1.*TMath::Pi())) dphi += 2*TMath::Pi();
  if (dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
  return dphi;
}

void HGCPhotonReco::fillzPositions(const edm::PtrVector<HGCRecHit>& rechitvec){
  for (unsigned iH(0); iH<rechitvec.size(); ++iH){
    const HGCRecHit & lHit = *(rechitvec[iH]);
    const HGCEEDetId & hgcid = lHit.detid();
    unsigned layer = hgcid.layer()-1;
    if (layer >= nLayers_) {
      std::cout << " -- Warning! Wrong layer number: " << layer << " max is set to " << nLayers_ << std::endl;
      continue;
    }
    GlobalPoint cellPos = hgcEEGeom_->getPosition(hgcid);
    zPos_.insert(std::pair<unsigned,double>(layer,fabs(cellPos.z())));
    if (zPos_.size()==nLayers_) break;
  }
}

void
HGCPhotonReco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  static bool first = true;

  if (debug_) std::cout << " -- Processing event " << iEvent.run() << " " << iEvent.luminosityBlock() << " " << iEvent.id().event() << std::endl;

  using namespace edm;
  
  Handle<edm::View<reco::SuperCluster> > HGCEESCs;
  iEvent.getByToken(endcapSuperClusterCollection_,HGCEESCs);
  const edm::PtrVector<reco::SuperCluster>& sclusters = HGCEESCs->ptrVector();
  
  Handle<edm::View<reco::PFCluster> > HGCEEClusters;
  iEvent.getByToken(endcapClusterCollection_,HGCEEClusters);
  const edm::PtrVector<reco::PFCluster>& clusters = HGCEEClusters->ptrVector();
  
  Handle<edm::View<reco::GenParticle> > genParts;
  iEvent.getByToken(genParticlesCollection_,genParts);
  const edm::PtrVector<reco::GenParticle>& gens = genParts->ptrVector();
  
  Handle<edm::View<HGCRecHit> > eeRecHits;
  iEvent.getByToken(endcapRecHitCollection_, eeRecHits);
  const edm::PtrVector<HGCRecHit>& rechitvec = eeRecHits->ptrVector();
  iEvent.getByLabel(edm::InputTag("HGCalRecHit:HGCEERecHits"),recHits_);

  if (zPos_.size()!=nLayers_) {
    fillzPositions(rechitvec);
    std::cout << " -- z Positions of layers: " << std::endl;
    for (unsigned iL(0); iL<nLayers_;++iL){
      std::cout << " ----- " << iL << " " << zPos_[iL] << std::endl;
    }
  }


  PCAShowerAnalysis pcaShowerAnalysis1(iEvent,iSetup);
  PCAShowerAnalysis pcaShowerAnalysis2(iEvent,iSetup);

  	
  //fill hit histos
  /*  for (unsigned iH(0); iH<rechitvec.size(); ++iH){
    const HGCRecHit & lHit = *(rechitvec[iH]);
    const HGCEEDetId & hgcid = lHit.detid();
    unsigned layer = hgcid.layer()-1;
    if (layer >= nLayers_) {
      std::cout << " -- Warning! Wrong layer number: " << layer << " max is set to " << nLayers_ << std::endl;
      continue;
    }
    GlobalPoint cellPos = hgcEEGeom_->getPosition(hgcid);
    double posx = cellPos.x();
    double posy = cellPos.y();
    //double posz = cellPos.z();
    //correction for deposited energy: more with more length in Si
    //double costheta = fabs(posz)/sqrt(posz*posz+posx*posx+posy*posy);
    //double energy = lHit.energy()/mipE_*costheta;//in MIP
    hyvsx_[layer]->Fill(posx,posy);//,energy);
    //hetavsphi_[layer]->Fill(cellPos.phi(),cellPos.eta());//,energy);
    // hphisecvscellid_[layer]->Fill(hgcid.cell(),hgcid.sector());
    //hyvsxzoom_[layer]->Fill(posx,posy);//,energy);
    }*/
  

  //Geant4 collections
  edm::Handle<edm::SimTrackContainer> SimTk;
  iEvent.getByLabel(g4TracksSource_,SimTk);
  edm::Handle<edm::SimVertexContainer> SimVtx;
  iEvent.getByLabel(g4VerticesSource_,SimVtx);
  
  
  //if (rechitvec.size());
  
  // initialise tree entries
  //photon 1
  info1_.initialise();
  //photon 2
  if (!singleGamma_) info2_.initialise();
  
  assert(gens.size() >0); // only the case for the electron gun sample
  
  edm::Ptr<reco::GenParticle> photon1;
  edm::Ptr<reco::GenParticle> photon2;
  unsigned nPhotons = 0;

  if (debug_) std::cout << "[debug] Number of gen particles: " << gens.size() << std::endl;

  for (unsigned int igp =0; igp < gens.size() ; igp++) { // loop over gen particles to fill truth-level tree
    
    if (debug_>1) {
      std::cout << "[debug] gen pdgid " << gens[igp]->pdgId() << ", status " << gens[igp]->status() << ", e " << gens[igp]->energy() << ", pt " << gens[igp]->pt() << ", eta = " << gens[igp]->eta();
      if (gens[igp]->mother()) std::cout << " mother " << (gens[igp]->mother()->pdgId());
      std::cout <<  std::endl;
    }
    if ((gens.size()>2 && (gens[igp]->pdgId() != 22 || gens[igp]->status() != 3)) || (gens.size()==2 && (gens[igp]->pdgId() != 22 || gens[igp]->status() != 1)) ) continue;
    //std::cout << " theta check: " << gens[igp]->theta() << " " << 2*atan(exp(-1.*gens[igp]->eta())) << std::endl;

    if (photon1.isNull()) {
      photon1 = gens[igp];
      truthVtx_ = gens[igp]->vertex();
    }
    else if (!singleGamma_) {
      photon2 = gens[igp];
      if (gens[igp]->vertex()!=truthVtx_){
	std::cout << " -- Error, photons don't have the same vertex !" << std::endl;
	exit(1);
      }
    }
    nPhotons++;
  }

  if ((!singleGamma_ && nPhotons!=2) || (singleGamma_ && nPhotons!=1)) {
    std::cout << " -- Error ! Found " << nPhotons << " photons status 3." <<  std::endl;
    exit(1);
  }
  
  //get closest supercluster
  float dRBest1 =999;
  float dRBest2 =999;
  int idx1 = -1;
  int idx2 = -1;
  
  if (sclusters.size()==0) {
    if (debug_) std::cout << " -- no clusters ! Skipping" << std::endl;
    nNoClusters_++;
    return;
  }
  
  for (unsigned int isc =0; isc < sclusters.size() ; isc++){ //subloop over sc's to find matches
    
    // calculate dR...
    float dE = sclusters[isc]->eta() - photon1->eta();
    dE =dE*dE;
    float dP = DeltaPhi(sclusters[isc]->phi(),photon1->phi());
    dP =dP*dP;
    float dR = sqrt(dE +dP);
    
    if (dR < dRBest1) { // only true if dR is both below limit value and smaller than previous best dR.
      dRBest1 = dR;
      idx1 = isc;
    }
    if (!singleGamma_){
      dE = sclusters[isc]->eta() - photon2->eta();
      dE =dE*dE;
      dP = DeltaPhi(sclusters[isc]->phi(),photon2->phi());
      dP =dP*dP;
      dR = sqrt(dE +dP);
      
      if (dR < dRBest2) { // only true if dR is both below limit value and smaller than previous best dR.
	dRBest2 = dR;
	idx2 = isc;
      }
    }
  }
  
  info1_.dRTrue = dRBest1;
  if (!singleGamma_) info2_.dRTrue = dRBest2;
  
  if (isInFid(photon1)){
    //process photons
    if (debug_) std::cout << " -------------------------- Photon1: " << std::endl;
    //fill SC info
    info1_.isValid = true;
    info1_.fillTruth(photon1, SimTk, SimVtx);  
    info1_.fillSC(sclusters[idx1],pcaShowerAnalysis1);

    bool pass = info1_.converted==0 && info1_.eTrue/cosh(info1_.etaTrue)>40 && TMath::Abs(info1_.etaTrue)>1.6 && TMath::Abs(info1_.etaTrue)<2.8 && (TMath::Nint(TMath::Abs(info1_.phiTrue)/(2*TMath::Pi())*360.)%20<9 || TMath::Nint(TMath::Abs(info1_.phiTrue)/(2*TMath::Pi())*360.)%20>11)  && info1_.noPhiCrackPCA==1;
    if (pass && first){
      const reco::CaloCluster* clus = &(*(sclusters[idx1]->seed()));
      for (unsigned int ih=0;ih<clus->hitsAndFractions().size();++ih) {
	const DetId & id = (clus->hitsAndFractions())[ih].first ;
	HGCRecHitCollection::const_iterator theHit = recHits_->find(id);    
	if (id.det()!=DetId::Forward || id.subdetId()!=HGCEE) continue;
	GlobalPoint cellPos = hgcEEGeom_->getPosition(HGCEEDetId(id));
	hyvsx_[HGCEEDetId(id).layer()-1]->Fill(cellPos.x(),cellPos.y(),theHit->energy()/0.0000551);//energy in mips
	//std::cout << cellPos.z() << std::endl;
	hxvsz_->Fill(cellPos.z()-(cellPos.z()>0?320:-320.),cellPos.x(),theHit->energy()/0.0000551);//energy in mips
	hyvsz_->Fill(cellPos.z()-(cellPos.z()>0?320:-320.),cellPos.y(),theHit->energy()/0.0000551);//energy in mips
	hrvsz_->Fill(cellPos.z()-(cellPos.z()>0?320:-320.),sqrt(pow(cellPos.x(),2)+pow(cellPos.y(),2)),theHit->energy()/0.0000551);//energy in mips
	
      }
      //fill expected position from PCA
      double theta = 2*atan(exp(-1.*info1_.etaPCA));
      hxvsz_->Fill(info1_.zPCA-(info1_.zPCA>0?320:-320.),info1_.xPCA,20000);
      hyvsz_->Fill(info1_.zPCA-(info1_.zPCA>0?320:-320.),info1_.yPCA,20000);
      for (unsigned iL(0);iL<nLayers_;++iL){
	double r = (info1_.etaPCA/fabs(info1_.etaPCA)*zPos_[iL]-info1_.zPCA)/cos(theta);
	hyvsx_[iL]->Fill(r*sin(theta)*cos(info1_.phiPCA)+info1_.xPCA,r*sin(theta)*sin(info1_.phiPCA)+info1_.yPCA,10000);
	hxvsz_->Fill(info1_.etaPCA/fabs(info1_.etaPCA)*(zPos_[iL]-320),r*sin(theta)*cos(info1_.phiPCA)+info1_.xPCA,10000);//energy in mips
	hyvsz_->Fill(info1_.etaPCA/fabs(info1_.etaPCA)*(zPos_[iL]-320),r*sin(theta)*sin(info1_.phiPCA)+info1_.yPCA,10000);//energy in mips
	hrvsz_->Fill(info1_.etaPCA/fabs(info1_.etaPCA)*(zPos_[iL]-320),sqrt(pow(r*sin(theta)*cos(info1_.phiPCA)+info1_.xPCA,2)+pow(r*sin(theta)*sin(info1_.phiPCA)+info1_.yPCA,2)),10000);
      }
      std::cout << " - Found photon to save!" << std::endl;
      first = false;
    }

    info1_.eSC = resumEmEnergyTest(sclusters[idx1],clusters);
    getPhotonEnergy(rechitvec,info1_);
  }
  if (!singleGamma_ && isInFid(photon2)){
    if (debug_) std::cout << " -------------------------- Photon2: " << std::endl;
    info2_.isValid = true;
    info2_.fillTruth(photon2, SimTk, SimVtx);
    info2_.fillSC(sclusters[idx2],pcaShowerAnalysis2);
    info2_.eSC = resumEmEnergyTest(sclusters[idx2],clusters);
    getPhotonEnergy(rechitvec,info2_);
  }
  
  //info1_.eSeed = (clusterEmEnergy(sclusters[idx1]->seed(), clusters));
  //info1_.eReco = resumEmEnergy(sclusters[info.matchIndex], clusters);
  //info1_.eRecoPCA = resumEmEnergyTest(sclusters[info.matchIndex], clusters);
  
  //auto 	sc = sclusters[info.matchIndex];
  //	for(unsigned int l =0 ; l<sc->clusters().size(); l++){
  //	}
  
  
  if (info1_.isValid || (!singleGamma_ && info2_.isValid)) tree_->Fill();
  return ;
}

float HGCPhotonReco::resumEmEnergyTest( const edm::Ptr<reco::SuperCluster>& sc,const edm::PtrVector<reco::PFCluster>& clusters){
  
  float total=0;
  for (unsigned int ic =0 ; ic < sc->clusters().size() ; ic++){
    //	//std::cout << "TEST, sc constituent em energies " << (sc->clusters())[ic]->energy() << std::endl;
    //for (unsigned int j =0 ; j < clusters.size() ; j++){
      
      double clus_eta = sc->clusters()[ic]->eta();
      double clus_phi = sc->clusters()[ic]->phi();
      
      if (fabs(clus_eta - sc->seed()->eta())>0.025) continue;
      if (fabs(deltaPhi(clus_phi , sc->seed()->phi()))>0.11) continue;
      
      //if (clusters[j]->position()==(sc->clusters())[ic]->position()) {
	//		//std::cout << "TEST, corresponding cluster " << (clusters[j]->emEnergy()) << std::endl;
      total = total +(sc->clusters()[ic]->energy());
      //break;
      //}
      //}
  
  }
  
  return total;
}

bool HGCPhotonReco::isInFid(const edm::Ptr<reco::GenParticle> & aPhoton){
  return fabs(aPhoton->eta())>1.4 &&  fabs(aPhoton->eta())<3.0;// && ((!singleGamma_ && aPhoton->pt() > 20) || singleGamma_);
}

void HGCPhotonReco::getPhotonEnergy(const edm::PtrVector<HGCRecHit>& rechitvec,
				    info_t & info){

  double phimax = info.phiTrue;//SC;
  double etamax = info.etaTrue;//SC;
  getPhotonEnergy(rechitvec,info,phimax,etamax,false);

  phimax = info.phiPCA;//SC;
  etamax = info.etaPCA;//SC;
  getPhotonEnergy(rechitvec,info,phimax,etamax,true);

}

void HGCPhotonReco::getPhotonEnergy(const edm::PtrVector<HGCRecHit>& rechitvec,info_t & info,const double & phimax, const double & etamax, const bool isPCA){
  
 //get central position for 3x3 or 5x5 arrays
  std::vector<HGCEEDetId> detidmax;
  detidmax.resize(nLayers_,HGCEEDetId());

  //from rechits
  //getMaximumCell(rechitvec,phimax,etamax,detidmax);
  //if (debug_) std::cout << " -- True eta-phi " << etamax << " " << phimax << std::endl;
  //for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
  //std::cout << iL << " Max=(" << xmax[iL] << "," << ymax[iL] << ")" << std::endl;
  //}
  
  //get PU contribution from around the photon
  
  //get energy-weighted position and energy around maximum
  //std::vector<HGCEEDetId> recoPos;
  //recoPos.resize(nLayers_,DetId());
  //std::vector<double> recoE;
  //recoE.resize(nLayers_,0);
  const HGCalTopology& topology = hgcEEGeom_->topology();
  double maxE = 0;
  double theta = 2*atan(exp(-1.*etamax));
  double x0=isPCA? info.xPCA: truthVtx_.x();
  double y0=isPCA? info.yPCA: truthVtx_.y();
  double z0=isPCA? info.zPCA: truthVtx_.z();
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    //get detidmax
    double r = (etamax/fabs(etamax)*zPos_[iL]-z0)/cos(theta);
    GlobalPoint max(r*sin(theta)*cos(phimax)+x0,r*sin(theta)*sin(phimax)+y0,etamax/fabs(etamax)*zPos_[iL]);
    detidmax[iL] = hgcEEGeom_->getClosestCell(max);
    /*const HGCEEDetId & detidcheck = hgcEEGeom_->getClosestCell(max);
    if (debug_ && detidcheck != detidmax[iL]){
      std::cout << " - layer " << iL << std::endl
		<< " -- from rechits: " << detidmax[iL] << " " 
		<< (hgcEEGeom_->getPosition(detidmax[iL])).x() << " "
		<< (hgcEEGeom_->getPosition(detidmax[iL])).y() << " "
		<< (hgcEEGeom_->getPosition(detidmax[iL])).z() << std::endl
		<< " -- from geom: " << detidcheck << " " 
		<< (hgcEEGeom_->getPosition(detidcheck)).x() << " "
		<< (hgcEEGeom_->getPosition(detidcheck)).y() << " "
		<< (hgcEEGeom_->getPosition(detidcheck)).z() << std::endl
		<< " -- diff from truth: " << ((hgcEEGeom_->getPosition(detidcheck)).x()-max.x())*10 << " " << 10*((hgcEEGeom_->getPosition(detidcheck)).y()-max.y()) << " " << 10*((hgcEEGeom_->getPosition(detidcheck)).z()-max.z())
		<< std::endl;
		}*/
    if (!topology.valid(detidmax[iL]) || detidmax[iL].det()!=DetId::Forward || detidmax[iL].subdetId()!=HGCEE) {
      info.invalidDetid=true;
      info.invalidDetidPCA=true;
      continue;
    }
    //loop on signal regions
    std::vector<unsigned> sizeR;
    sizeR.push_back(3);
    sizeR.push_back(5);
    sizeR.push_back(7);
    double etot = 0;
    for (unsigned iS(0);iS<sizeR.size();++iS){
      
      std::vector<double> Exy;
      Exy.resize(sizeR[iS]*sizeR[iS],0);
      if (debug_) std::cout << " -- detidmax[" << iL << "=" << detidmax[iL] << std::endl;
      fillNeighbours(detidmax[iL],Exy,info,isPCA);
      const GlobalPoint & center = hgcEEGeom_->getPosition(detidmax[iL]);
      if (!info.converted) {
	dxMaxTruth_->Fill((center.x()-max.x())*10);
	dyMaxTruth_->Fill((center.y()-max.y())*10);
	dzMaxTruth_->Fill((center.z()-max.z())*10);
	dphiMaxTruthvsLayer_->Fill(iL,DeltaPhi(center.phi(),max.phi()));
      }
      for (unsigned idx(0);idx<Exy.size();++idx){
	if (iS==0) etot += Exy[idx];
	if (!isPCA) {
	  info.eReco[iS] += Exy[idx]*absWeight(iL);
	}
	else info.eRecoPCA[iS] += Exy[idx]*absWeight(iL);
      }
    }
    if (etot>maxE){
      maxE = etot;
      if (!isPCA) info.showerMax = iL;
    }
  }//loop on layers

  GlobalPoint cellPos = hgcEEGeom_->getPosition(detidmax[info.showerMax]);
  if (!isPCA) {
    info.zShowerMax = cellPos.z();
    double theta = 2*atan(exp(-1.*fabs(info.etaTrue)));
    double thetadet =atan((info.zShowerMax-info.zvtxTrue)*tan(theta)/info.zShowerMax); 
    info.etadetTrue = info.zShowerMax/fabs(info.zShowerMax)*-1.*log(tan(thetadet/2.));
    //std::cout << info.zShowerMax 
    //<< " " << info.etaTrue
    //	      << " " << theta 
    //	      << " " << thetadet
    //	      << " " << info.etadetTrue
    //	      << std::endl;
  }
  //getEnergyWeightedPosition(rechitvec,detidmax,recoPos,recoE);
  //getTotalEnergy(rechitvec,recoPos,info);
  info.correctForEtaDep(ROOT::Math::XYZVector(cellPos.x(),cellPos.y(),cellPos.z()),isPCA);
}

void HGCPhotonReco::getMaximumCell(const edm::PtrVector<HGCRecHit>& rechitvec,const double & phimax,const double & etamax,std::vector<HGCEEDetId> & detidmax){

  std::vector<double> dRmin;
  dRmin.resize(nLayers_,10);
  if (debug_) std::cout << " - Processing " << rechitvec.size() << " rechits " << std::endl;
  for (unsigned iH(0); iH<rechitvec.size(); ++iH){//loop on rechits
    const HGCRecHit & lHit = *(rechitvec[iH]);
    const HGCEEDetId & hgcid = lHit.detid();
    unsigned layer = hgcid.layer()-1;
    if (layer >= nLayers_) {
      std::cout << " -- Warning! Wrong layer number: " << layer << " max is set to " << nLayers_ << std::endl;
      continue;
    }
    GlobalPoint cellPos = hgcEEGeom_->getPosition(hgcid);
    double posx = cellPos.x();
    double posy = cellPos.y();
    double posz = cellPos.z();
    //    if (debug_>1) {
    //std::cout << " --  RecoHit " << iH << "/" << rechitvec.size() << " -- layer " << layer << " id " << hgcid << std::endl
    //	      << " --  position x,y,z " << posx << "," << posy << "," << posz << std::endl;
    //}
    //double energy = lHit.energy();
    
    ROOT::Math::XYZVector pos(posx-truthVtx_.x(),posy-truthVtx_.y(),posz-truthVtx_.z());
    double deta = fabs(pos.eta()-etamax);
    double dphi = DeltaPhi(pos.phi(),phimax);
    
    double dR = sqrt(pow(deta,2)+pow(dphi,2));
    if (dR<dRmin[layer]) {
      dRmin[layer] = dR;
      detidmax[layer] = hgcid;
    }
    
    
  }//loop on rechits
  
  for (unsigned iL(0);iL<nLayers_;++iL){
    dRmin_[iL]->Fill(dRmin[iL]);
  }

}

/*void HGCPhotonReco::getEnergyWeightedPosition(
		const edm::PtrVector<HGCRecHit>& rechitvec,
		const std::vector<HGCEEDetId> & detidmax,
		std::vector<HGCEEDetId> & recoPos,
		std::vector<double> & recoE
		//,std::vector<double> & puE,
		//const bool puSubtracted
		)
{
  std::vector<unsigned> nHits;
  nHits.resize(nLayers_,0);

  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    std::vector<double> Exy;
    Exy.resize(9,0);
    fillNeighbours(detidmax[iL],Exy);
    GlobalPoint cellPos = hgcEEGeom_->getPosition(detidmax[iL]);
    double xmax = cellPos.x();
    double ymax = cellPos.y();
    if (nHits[iL]==0) continue;

    double Etot = 0;
    double Ex[3] = {0,0,0};
    double Ey[3] = {0,0,0};
    for (unsigned idx(0);idx<9;++idx){
      Etot += Exy[idx];
    }
    recoE[iL] = Etot;
    hEvsLayer_->Fill(iL,Etot);
    
    Ex[0] = Exy[0]+Exy[3]+Exy[6];
    Ex[1] = Exy[1]+Exy[4]+Exy[7];
    Ex[2] = Exy[2]+Exy[5]+Exy[8];
    Ey[0] = Exy[0]+Exy[1]+Exy[2];
    Ey[1] = Exy[3]+Exy[4]+Exy[5];
    Ey[2] = Exy[6]+Exy[7]+Exy[8];
    double wx[4];
    double wy[4];
    for (unsigned i(0);i<4;++i){
      wx[i] = 0;
      wy[i] = 0;
    }
    double w0 = getW0(iL);
    for (unsigned i(0);i<3;++i){
      wx[i] = std::max(0.,log(Ex[i]/Etot)+w0);
      wy[i] = std::max(0.,log(Ey[i]/Etot)+w0);
      wx[3] += wx[i];
      wy[3] += wy[i];
    }
    double x = xmax;
    //if none pass, discard layer
    if (wx[3]!=0) x += cellSize_*(wx[2]-wx[0])/wx[3];
    else nHits[iL]=0;
    double y = ymax;
    if (wy[3]!=0) y += cellSize_*(wy[2]-wy[0])/wy[3];
    else nHits[iL]=0;
    
    if (nHits[iL]!=0){
      recoPos[iL].SetX(x);
      recoPos[iL].SetY(y);
    }
  }//loop on layers
  
}
*/
void HGCPhotonReco::fillNeighbours(const HGCEEDetId & detidmax,
				   std::vector<double> & Exy,
				   info_t & info,
				   bool isPCA){

  unsigned nCells = sqrt(Exy.size());

  const HGCalTopology& topology = hgcEEGeom_->topology();
  std::vector<HGCEEDetId> neighbours;
  neighbours.resize(Exy.size(),HGCEEDetId());
  neighbours[(Exy.size()-1)/2] = detidmax;
  unsigned nGo = (nCells-1)/2;

  //south lines
  for (unsigned i(0); i<nGo;++i){
    unsigned center = (Exy.size()-1)/2 -nCells*i;
    //std::cout << " Center " << center << std::endl;
    unsigned idxrow = (Exy.size()-1)/2 -nCells*(i+1);
    DetId tmp = topology.goSouth(neighbours[center]);
    if (topology.valid(tmp)) {
      neighbours[idxrow] = HGCEEDetId(tmp);
      unsigned idx=idxrow;
      //std::cout << " South row center " << idx << " " ;
      for (unsigned j(0); j<nGo;++j){
	unsigned nextidx = (Exy.size()-1)/2 -nCells*(i+1) - (j+1);
	//std::cout << nextidx << " " ;
	tmp = topology.goWest(neighbours[idx]);
	if (topology.valid(tmp)) neighbours[nextidx]  = HGCEEDetId(tmp);
	idx=nextidx;
      }
      idx = idxrow;
      for (unsigned j(0); j<nGo;++j){
	unsigned nextidx = (Exy.size()-1)/2 -nCells*(i+1) + (j+1);
	//std::cout << nextidx << " ";
	tmp = topology.goEast(neighbours[idx]);
	if (topology.valid(tmp)) neighbours[nextidx] = HGCEEDetId(tmp);
	idx=nextidx;
      }
      //std::cout << std::endl;
    }
  }

  {
    //center line
    unsigned idx=(Exy.size()-1)/2;
    //std::cout << " Central row center " << idx << " " ;
    for (unsigned j(0); j<nGo;++j){
      unsigned nextidx = (Exy.size()-1)/2 - (j+1);
      //std::cout << nextidx << " " ;
      DetId tmp = topology.goWest(neighbours[idx]);
      if (topology.valid(tmp)) neighbours[nextidx]  = HGCEEDetId(tmp);
      idx=nextidx;
    }
    idx = (Exy.size()-1)/2;
    for (unsigned j(0); j<nGo;++j){
      unsigned nextidx = (Exy.size()-1)/2 + (j+1);
      //std::cout << nextidx << " ";
      DetId tmp = topology.goEast(neighbours[idx]);
      if (topology.valid(tmp)) neighbours[nextidx] = HGCEEDetId(tmp);
      idx=nextidx;
    }
    //std::cout << std::endl;
  }

  //north
  for (unsigned i(0); i<nGo;++i){
    unsigned center = (Exy.size()-1)/2 + nCells*i;
    //std::cout << " Center " << center << std::endl;
    unsigned idxrow = (Exy.size()-1)/2 + nCells*(i+1);
    DetId tmp = topology.goNorth(neighbours[center]);
    if (topology.valid(tmp)) {
      neighbours[idxrow] = HGCEEDetId(tmp);
      unsigned idx=idxrow;
      //std::cout << " North row center " << idx << " " ;
      for (unsigned j(0); j<nGo;++j){
	unsigned nextidx = (Exy.size()-1)/2 + nCells*(i+1) - (j+1);
	//std::cout << nextidx << " " ;
	tmp = topology.goWest(neighbours[idx]);
	if (topology.valid(tmp)) neighbours[nextidx]  = HGCEEDetId(tmp);
	idx=nextidx;
      }
      idx = idxrow;
      for (unsigned j(0); j<nGo;++j){
	unsigned nextidx = (Exy.size()-1)/2 + nCells*(i+1) + (j+1);
	//std::cout << nextidx << " ";
	tmp = topology.goEast(neighbours[idx]);
	if (topology.valid(tmp)) neighbours[nextidx] = HGCEEDetId(tmp);
	idx=nextidx;
      }
      //std::cout << std::endl;
    }
  }
  
  GlobalPoint center = hgcEEGeom_->getPosition(detidmax);

  for (unsigned idx(0);idx<Exy.size();++idx){
    //std::cout << " idx " << idx << std::flush;
    //std::cout << " neighbour detid " << neighbours[idx] << std::endl;
    if (!topology.valid(neighbours[idx]) || neighbours[idx].det()!=DetId::Forward || neighbours[idx].subdetId()!=HGCEE) {
      if (!isPCA) info.invalidNeighbour = true;
      else info.invalidNeighbourPCA = false;
      continue;
    }
    //std::cout << " accessing rechit " << std::endl;
    HGCRecHitCollection::const_iterator theHit = recHits_->find(neighbours[idx]);
    if (theHit==recHits_->end()) continue;
    //std::cout << " rechit found" << std::endl;
    GlobalPoint cellPos = hgcEEGeom_->getPosition(neighbours[idx]);
    if ((neighbours[idx].sector()*neighbours[idx].subsector()) != (detidmax.sector()*detidmax.subsector())) {
      if (!isPCA) info.noPhiCrack = false;
      else info.noPhiCrackPCA = false;
    }
    double posx = cellPos.x();
    double posy = cellPos.y();
    double posz = cellPos.z();
    if (Exy.size()==9){
      dxMaxNeigh_[idx]->Fill((posx-center.x())*10.);
      dyMaxNeigh_[idx]->Fill((posy-center.y())*10.);
      dzMaxNeigh_[idx]->Fill((posz-center.z())*10.);
    }
    if (fabs(posx-center.x())>sqrt(2*pow(nGo,2)) ||
	fabs(posy-center.y())>sqrt(2*pow(nGo,2))){
      if (!isPCA) info.invalidNeighbour = true;
      else info.invalidNeighbourPCA = false;
    }
    //double costheta = fabs(posz)/sqrt(posz*posz+posx*posx+posy*posy);
    double energy = theHit->energy()/mipE_;//*costheta;//in MIP
    Exy[idx] = energy;
    //Exyetacor[idx] = energy*costheta;
    //std::cout << " -end processing neighbour" << std::endl;
  }
}

//void HGCPhotonReco::getTotalEnergy(const edm::PtrVector<HGCRecHit>& rechitvec,
//const std::vector<HGCEEDetId> & eventPos,
//info_t & info)
//{

  //double dx = eventPos[layer].x()-posx;
  //double dy = eventPos[layer].y()-posy;
  //double halfCell = 0.5*cellSize_;
    
  //for (unsigned isr(0); isr<nSR_;++isr){
  //if ( (fabs(dx) <= ((isr+1)*halfCell)) && (fabs(dy) <= ((isr+1)*halfCell))){
	//correction for absorber thickness at 0 angle
//info.eReco[isr] += energy*absWeight(layer);
	//}
	//}
    
  
//}

double HGCPhotonReco::getW0(const unsigned layer){
	if (layer<7) return 4;
	if (layer==7) return 2.55;
	if (layer==8) return 2.9;
	if (layer==9) return 2.45;
	if (layer==10) return 2.75;
	if (layer==11) return 2.35;
	if (layer==12) return 2.55;
	if (layer==13) return 2.2;
	if (layer==14) return 2.35;
	if (layer==15) return 2;
	if (layer==16) return 2.2;
	if (layer==17) return 1.9;
	if (layer==18) return 2.05;
	if (layer==19) return 1.75;
	if (layer==20) return 1.9;
	if (layer==21) return 1.7;
	if (layer==22) return 1.8;
	if (layer==23) return 3;
	if (layer>23) return 4;
	return 0;
}
double HGCPhotonReco::absWeight(const unsigned layer, const bool dedx){
	if (dedx==false){
		if (layer == 0) return 0.08696;
		if (layer == 1) return 1;//0.92
		if (layer == 2) return 0.646989;//88.16/95.4=0.92
		if (layer == 3) return 0.617619;//51.245/95.4=0.537
		if (layer == 4) return 0.646989;
		if (layer == 5) return 0.617619;
		if (layer == 6) return 0.646989;
		if (layer == 7) return 0.617619;
		if (layer == 8) return 0.646989;
		if (layer == 9) return 0.617619;
		if (layer == 10) return 0.646989;
		if (layer == 11) return 0.942829;//74.45/95.4=0.78
		if (layer == 12) return 0.859702;//102.174/95.4=1.071
		if (layer == 13) return 0.942829;
		if (layer == 14) return 0.859702;
		if (layer == 15) return 0.942829;
		if (layer == 16) return 0.859702;
		if (layer == 17) return 0.942829;
		if (layer == 18) return 0.859702;
		if (layer == 19) return 0.942829;
		if (layer == 20) return 0.859702;
		if (layer == 21) return 1.37644;//105.39/95.4=1.1047
		if (layer == 22) return 1.30447;//131.476/95.4=1.378
		if (layer == 23) return 1.37644;
		if (layer == 24) return 1.30447;
		if (layer == 25) return 1.37644;
		if (layer == 26) return 1.30447;
		if (layer == 27) return 1.37644;
		if (layer == 28) return 1.30447;
		if (layer == 29) return 1.37644;//1.79662;//
	}
	else {
		if (layer == 0) return 0.06588;
		if (layer == 1) return 1;//95.4/95.4=1
		if (layer == 2) return 0.92;//88.16/95.4=0.92
		if (layer == 3) return 0.537;//51.245/95.4=0.537
		if (layer == 4) return 0.92;
		if (layer == 5) return 0.537;
		if (layer == 6) return 0.92;
		if (layer == 7) return 0.537;
		if (layer == 8) return 0.92;
		if (layer == 9) return 0.537;
		if (layer == 10) return 0.92;
		if (layer == 11) return 0.78;//74.45/95.4=0.78
		if (layer == 12) return 1.071;//102.174/95.4=1.071
		if (layer == 13) return 0.78;
		if (layer == 14) return 1.071;
		if (layer == 15) return 0.78;
		if (layer == 16) return 1.071;
		if (layer == 17) return 0.78;
		if (layer == 18) return 1.071;
		if (layer == 19) return 0.78;
		if (layer == 20) return 1.071;
		if (layer == 21) return 1.1047;//105.39/95.4=1.1047
		if (layer == 22) return 1.378;//131.476/95.4=1.378
		if (layer == 23) return 1.1047;
		if (layer == 24) return 1.378;
		if (layer == 25) return 1.1047;
		if (layer == 26) return 1.378;
		if (layer == 27) return 1.1047;
		if (layer == 28) return 1.378;
		if (layer == 29) return 1.1047;
	}
	return 1;
}


// ------------ method called once each job just before starting event loop  ------------
	void 
HGCPhotonReco::beginJob()
{
	nNoClusters_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
HGCPhotonReco::endJob() 
{
	std::cout << " -- EndJob:" << std::endl
		<< " -- Number of events with 0 clusters: " << nNoClusters_ << std::endl;
}

// ------------ method called when starting to processes a run  ------------
	void 
HGCPhotonReco::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
  edm::ESHandle<HGCalGeometry> hgcGeo;
  iSetup.get<IdealGeometryRecord>().get(geometrySource_[0],hgcGeo);
  hgcEEGeom_=hgcGeo.product();
}

// ------------ method called when ending the processing of a run  ------------
	void 
HGCPhotonReco::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
HGCPhotonReco::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
HGCPhotonReco::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCPhotonReco::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCPhotonReco);
