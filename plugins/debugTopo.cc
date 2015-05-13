// -*- C++ -*-
//
// Package:    debugTopo
// Class:      debugTopo
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

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"
#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "TTree.h"
#include "TMath.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraphErrors.h"

#include "UserCode/HGCanalysis/interface/HGCAnalysisTools.h"
using namespace std;

//
// class declaration
//
// .h class info
class debugTopo : public edm::EDAnalyzer {
public:
  explicit debugTopo(const edm::ParameterSet&);
  ~debugTopo();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  void fillNeighbours(const HGCEEDetId & detidmax,
		      const bool isP,
		      const unsigned nInEta,
		      const unsigned nInPhi);

private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


  std::vector<std::string> geometrySource_;
  const HGCalGeometry * hgcEEGeom_;

  unsigned nLayers_;
  unsigned nSectors_;
  unsigned nCells_;
  double cellSize_;
  unsigned debug_;
  TH2F *hyvsxn_[10][5];
  TH2F *hyvsxp_[10][5];


};

// constructor
debugTopo::debugTopo(const edm::ParameterSet& iConfig):
  debug_(iConfig.getParameter<unsigned>("debug"))
{ 
  geometrySource_ = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");

  //@TODO
  //get this from geometry !!
  nLayers_ = 10;
  nSectors_ = 5;
  nCells_ = 2386;
  cellSize_ = 1;
  edm::Service<TFileService> fs_;
  
  
  for (unsigned iL(0);iL<nLayers_;++iL){
    for (unsigned iS(0);iS<nSectors_;++iS){
      std::ostringstream label;
      label.str("");
      label << "hyvsxn_" << iL << "_" << iS;
      hyvsxn_[iL][iS] = fs_->make<TH2F>(label.str().c_str(),
				       ";x (cm);y (cm);",
				       1700,0,170,
				       2000,-30,170);
      label.str("");
      label << "hyvsxp_" << iL << "_" << iS;
      hyvsxp_[iL][iS] = fs_->make<TH2F>(label.str().c_str(),
				       ";x (cm);y (cm);",
				       1700,0,170,
				       2000,-30,170);
    }
  }
}

// destructor
debugTopo::~debugTopo()
{
  
}


//
// member functions
//
void
debugTopo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //static bool first = true;

  if (debug_) std::cout << " -- Processing event " << iEvent.run() << " " << iEvent.luminosityBlock() << " " << iEvent.id().event() << std::endl;

  const HGCalTopology& topology = hgcEEGeom_->topology();


  std::cout << " Total number of modules: " << topology.totalModules() << std::endl;
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    for (unsigned iS(0); iS<nSectors_;++iS){//loop on sectors
      //unsigned iC = ;
      for (unsigned iC(0); iC<nCells_;++iC){//loop on cells


	if ((iL!=0 || (iC != 38 && iC != 72 && iC != 122 && iC != 168 && iC != 2117 && iC != 2144)) &&
	    (iL!=9 || (iC != 31 && iC != 65 && iC != 111 && iC != 115 && iC != 170 && iC != 171 && iC != 2194 && iC != 2222)) &&
	    (iC!=0 && iC !=5 && iC != 23 && iC != 109)) continue;
	    

	HGCalTopology::DecodedDetId id;
	id.iCell = iC;
	id.iLay = iL+1;
	id.iSec = iS+1;
	id.iSubSec = -1;
	id.zside = 1;
	id.subdet = HGCEE;
	HGCEEDetId detidmax = topology.encode(id);
	if (!topology.valid(detidmax) || detidmax.det()!=DetId::Forward || detidmax.subdetId()!=HGCEE) {
	  std::cout << iL << " " << iS << " " << iC << " invalid detid " << detidmax << "!" << std::endl;
	  continue;
	}
	GlobalPoint center = hgcEEGeom_->getPosition(detidmax);
	//if (iL==0 && iS == 0 && iC<100){
	  //std::cout << " - " << iL << " " << iS << " " << iC << " xidx " << detidmax. << std::endl;
	hyvsxn_[iL][iS]->Fill(center.x(),center.y(),9);
	
	  //}
	//if (iL==0) std::cout << "Processing cell " << iL << " " << iS << "/-1 " << iC << std::endl;
	fillNeighbours(detidmax,false,3,7);
	//hyvsxn_[iL][iS]->Fill(center.x(),center.y());
	
	id.iSubSec = 1;
	detidmax = topology.encode(id);
	center = hgcEEGeom_->getPosition(detidmax);
	//if (iL==0) std::cout << "Processing cell " << iL << " " << iS << "/1 " << iC << std::endl;
	fillNeighbours(detidmax,true,3,7);
	//if (iL==0 && iS == 0 && iC<100){
	  //std::cout << " + " << iL << " " << iS << " " << iC << " detid " << detidmax << std::endl;
	hyvsxp_[iL][iS]->Fill(center.x(),center.y(),9);
	  //}
	
	}//loop on cells
    }//loop on sections
  }//loop on layers
  
  return ;
}



void debugTopo::fillNeighbours(const HGCEEDetId & detidmax,
			       const bool isPlus,
			       const unsigned nInEta,
			       const unsigned nInPhi){

  //std::cout << " ===================================" << std::endl
  //<< " -- fillneighbours for " << nInEta << " x " << nInPhi << std::endl
  //<< " ===================================" << std::endl;

  const HGCalTopology& topology = hgcEEGeom_->topology();
  std::vector<HGCEEDetId> neighbours;
  neighbours.resize(nInEta*nInPhi,HGCEEDetId());
  unsigned nGoE = 0;
  unsigned nGoW = 0;
  unsigned nGoN = 0;
  unsigned nGoS = 0;
  std::vector<unsigned> nGoMax;
  nGoMax.resize(nInEta*nInPhi,0);
  unsigned centralID = 0;
  //symmetric cells
  if (nInEta==nInPhi){
    centralID = (nInEta*nInPhi-1)/2;
    nGoE = (nInEta-1)/2;
    nGoW = nGoE;
    nGoN = nGoE;
    nGoS = nGoE;
  }
  else if (nInPhi>nInEta) {
    nGoN = (nInEta-1)/2;
    nGoS = nGoN;
    //get direction in which to go from phiSeed-phiSC:
    //and quadrant
    //if (DeltaPhi(info.phiSeed,info.phiSC)>0){
      //center on west side
      //go more east
      centralID = (nInEta-1)/2*nInPhi+(nInEta-1)/2;
      nGoE = (nInEta-1)/2+nInPhi-nInEta;
      nGoW = nGoN;
      //} else {
      //center of east side
      //go more west
      //centralID = (nInEta-1)/2*nInPhi+(nInPhi-1)-(nInEta-1)/2;
      //nGoW = (nInEta-1)/2+nInPhi-nInEta;
      //nGoE = nGoN;
      //}
  }
  else {
    std::cout << " -- Wrong indices nInEta>nInPhi !! Don't know what to do..." << std::endl;
    return;
  }
  neighbours[centralID] = detidmax;
  //south lines
  for (unsigned i(0); i<nGoS;++i){
    unsigned center = centralID -nInPhi*i;
    //std::cout << " Center " << center << std::endl;
    unsigned idxrow = centralID -nInPhi*(i+1);
    DetId tmp = topology.goSouth(neighbours[center]);
    if (topology.valid(tmp)) {
      neighbours[idxrow] = HGCEEDetId(tmp);
      nGoMax[idxrow] = i+1;
      unsigned idx=idxrow;
      //std::cout << " South row center " << idx << " " ;
      for (unsigned j(0); j<nGoW;++j){
	unsigned nextidx = centralID -nInPhi*(i+1) - (j+1);
	//std::cout << nextidx << " " ;
	int subsec = neighbours[idx].subsector();
	if (subsec>0) tmp = topology.goWest(neighbours[idx]);
	else  tmp = topology.goEast(neighbours[idx]);
	if (topology.valid(tmp)) {
	  neighbours[nextidx]  = HGCEEDetId(tmp);
	  nGoMax[nextidx] = i+1+j;
	}
	idx=nextidx;
      }
      idx = idxrow;
      for (unsigned j(0); j<nGoE;++j){
	unsigned nextidx = centralID -nInPhi*(i+1) + (j+1);
	//std::cout << nextidx << " ";
	int subsec = neighbours[idx].subsector();
	if (subsec>0) tmp = topology.goEast(neighbours[idx]);
	else tmp = topology.goWest(neighbours[idx]);
	if (topology.valid(tmp)) {
	  neighbours[nextidx] = HGCEEDetId(tmp);
	  nGoMax[nextidx] = i+1+j;
	}
	idx=nextidx;
      }
      //std::cout << std::endl;
    }
  }

  {
    //center line
    unsigned idx=centralID;
    //std::cout << " Central row center " << idx << " " ;
    for (unsigned j(0); j<nGoW;++j){
      unsigned nextidx = centralID - (j+1);
      //std::cout << nextidx << " " ;
      int subsec = neighbours[idx].subsector();
      DetId tmp;
      if (subsec>0) tmp = topology.goWest(neighbours[idx]);
      else tmp = topology.goEast(neighbours[idx]);
      if (topology.valid(tmp)) {
	neighbours[nextidx]  = HGCEEDetId(tmp);
	nGoMax[nextidx] = 1+j;
      }
      idx=nextidx;
    }
    idx = centralID;
    for (unsigned j(0); j<nGoE;++j){
      unsigned nextidx = centralID + (j+1);
      //std::cout << nextidx << " ";
      int subsec = neighbours[idx].subsector();
      DetId tmp;
      if (subsec>0) tmp = topology.goEast(neighbours[idx]);
      else  tmp = topology.goWest(neighbours[idx]);
      if (topology.valid(tmp)) {
	neighbours[nextidx] = HGCEEDetId(tmp);
	nGoMax[nextidx] = 1+j;
      }
      idx=nextidx;
    }
    //std::cout << std::endl;
  }

  //north
  for (unsigned i(0); i<nGoN;++i){
    unsigned center = centralID + nInPhi*i;
    //std::cout << " Center " << center << std::endl;
    unsigned idxrow = centralID + nInPhi*(i+1);
    DetId tmp = topology.goNorth(neighbours[center]);
    if (topology.valid(tmp)) {
      neighbours[idxrow] = HGCEEDetId(tmp);
      nGoMax[idxrow] = i+1;
      unsigned idx=idxrow;
      //std::cout << " North row center " << idx << " " ;
      for (unsigned j(0); j<nGoW;++j){
	unsigned nextidx = centralID + nInPhi*(i+1) - (j+1);
	//std::cout << nextidx << " " ;
	int subsec = neighbours[idx].subsector();
	if (subsec>0) tmp = topology.goWest(neighbours[idx]);
	else tmp = topology.goEast(neighbours[idx]);
	if (topology.valid(tmp)) {
	  neighbours[nextidx]  = HGCEEDetId(tmp);
	  nGoMax[nextidx] = i+1+j;
	}
	idx=nextidx;
      }
      idx = idxrow;
      for (unsigned j(0); j<nGoE;++j){
	unsigned nextidx = centralID + nInPhi*(i+1) + (j+1);
	//std::cout << nextidx << " ";
	int subsec = neighbours[idx].subsector();
	if (subsec>0) tmp = topology.goEast(neighbours[idx]);
	else tmp = topology.goWest(neighbours[idx]);
	if (topology.valid(tmp)) {
	  neighbours[nextidx] = HGCEEDetId(tmp);
	  nGoMax[nextidx] = i+1+j;
	}
	idx=nextidx;
      }
      //std::cout << std::endl;
    }
  }
  
  //GlobalPoint center = hgcEEGeom_->getPosition(detidmax);

  unsigned iL = detidmax.layer()-1;
  unsigned iS = detidmax.sector()-1;
  for (unsigned idx(0);idx<nInEta*nInPhi;++idx){
    //std::cout << " idx " << idx << std::flush;
    //std::cout << " neighbour detid " << neighbours[idx] << std::endl;
    if (!topology.valid(neighbours[idx]) || neighbours[idx].det()!=DetId::Forward || neighbours[idx].subdetId()!=HGCEE) {
      continue;
    }
    GlobalPoint center = hgcEEGeom_->getPosition(neighbours[idx]);
    //HGCalTopology::DecodedDetId id = topology.decode(neighbours[idx]);
    if (!isPlus && idx!=centralID) hyvsxn_[iL][iS]->Fill(center.x(),center.y(),idx+1);
    if (isPlus && idx!=centralID) hyvsxp_[iL][iS]->Fill(center.x(),center.y(),idx+1);
    //if ((neighbours[idx].sector()*neighbours[idx].subsector()) != (detidmax.sector()*detidmax.subsector())) {
    //}
    //double posx = cellPos.x();
    //double posy = cellPos.y();
    //double posz = cellPos.z();
    //if (fabs(posx-center.x())>sqrt(2*pow(nGoMax[idx],2)) ||
    //fabs(posy-center.y())>sqrt(2*pow(nGoMax[idx],2))){
    //}
  }
}

// ------------ method called once each job just before starting event loop  ------------
	void 
debugTopo::beginJob()
{
  //nNoClusters_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
debugTopo::endJob() 
{
  std::cout << " -- EndJob:" << std::endl;
	  //<< " -- Number of events with 0 clusters: " << nNoClusters_ << std::endl;
}

// ------------ method called when starting to processes a run  ------------
	void 
debugTopo::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
  edm::ESHandle<HGCalGeometry> hgcGeo;
  iSetup.get<IdealGeometryRecord>().get(geometrySource_[0],hgcGeo);
  hgcEEGeom_=hgcGeo.product();
}

// ------------ method called when ending the processing of a run  ------------
	void 
debugTopo::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
debugTopo::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
debugTopo::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
debugTopo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(debugTopo);
