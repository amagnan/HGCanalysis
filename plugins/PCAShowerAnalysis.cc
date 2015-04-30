//#include "RecoEgamma/Examples/interface/PCAShowerAnalysis.h"
#include "PCAShowerAnalysis.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "TMatrixD.h"
#include "TVectorD.h"

PCAShowerAnalysis::PCAShowerAnalysis (const edm::Event& iEvent, const edm::EventSetup& iSetup, 
 bool segmented, bool logweighting, bool debug) : geometry_(0), caloGeomCacheId_(0), principal_(0), logweighting_(logweighting), 
 segmented_(segmented), alreadyfilled_(false), debug_(debug)
{

  iEvent.getByLabel(edm::InputTag("HGCalRecHit:HGCEERecHits"),recHits_);
  principal_ = new TPrincipal(3,"D"); 

  edm::ESHandle<HGCalGeometry> pGeometry;
  unsigned long long newCaloGeomCacheId= iSetup.get<IdealGeometryRecord>().cacheIdentifier() ;
  if (caloGeomCacheId_!=newCaloGeomCacheId) {
    caloGeomCacheId_ = newCaloGeomCacheId ;
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",pGeometry) ;
    geometry_ = pGeometry.product() ;
  }  
  
  // minimal rechit value
  mip_ = 0.000055;//40;
  entryz_ = 320.38;

}

PCAShowerAnalysis::~PCAShowerAnalysis ()
{
  delete principal_;
}

void PCAShowerAnalysis::showerParameters(const reco::SuperCluster* clus, GlobalPoint& barycenter, GlobalVector& axis) 
{ 
  showerParameters(clus->seed(), barycenter, axis);
  return; 
}

GlobalPoint PCAShowerAnalysis::showerBarycenter(const reco::SuperCluster* clus)
{
  return showerBarycenter(clus->seed());
}

GlobalVector PCAShowerAnalysis::showerAxis(const reco::SuperCluster* clus)
{
  return showerAxis(clus->seed());
}


GlobalVector PCAShowerAnalysis::showerEigenValues(const reco::SuperCluster* clus)
{
  return showerEigenValues(clus->seed());
}


GlobalVector PCAShowerAnalysis::showerSigmas(const reco::SuperCluster* clus)
{
  return showerSigmas(clus->seed());
}

void PCAShowerAnalysis::showerParameters( const reco::CaloClusterPtr clus, GlobalPoint& barycenter, GlobalVector& axis)
{
  return showerParameters(&(*clus),barycenter,axis);
}

GlobalPoint PCAShowerAnalysis::showerBarycenter(const reco::CaloClusterPtr clus)
{
  return showerBarycenter(&(*clus));
}

GlobalVector PCAShowerAnalysis::showerAxis(const reco::CaloClusterPtr clus)
{
  return showerAxis(&(*clus));
}
  
GlobalVector PCAShowerAnalysis::showerEigenValues(const reco::CaloClusterPtr clus)
{
  return showerEigenValues(&(*clus));
}
  
GlobalVector PCAShowerAnalysis::showerSigmas(const reco::CaloClusterPtr clus)
{
  return showerSigmas(&(*clus));
}  

void PCAShowerAnalysis::showerParameters(const reco::CaloCluster* clus, GlobalPoint& barycenter, GlobalVector& axis)
{

  if (!alreadyfilled_) {
  
    double variables[3] = {0.,0.,0.};
  
    for (unsigned int ih=0;ih<clus->hitsAndFractions().size();++ih) {
      const DetId & id_ = (clus->hitsAndFractions())[ih].first ;
      HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
      if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	GlobalPoint cellPos = geometry_->getPosition(HGCEEDetId(id_));
	variables[0] = cellPos.x(); variables[1] = cellPos.y(); variables[2] = cellPos.z();
	if (!segmented_) variables[2] = entryz_;
        if (!logweighting_) {
	// energy weighting
	for (int i=0; i<int(theHit->energy()/mip_); i++) principal_->AddRow(variables); 
	} else {
	// a log-weighting, energy not in fraction of total
	double w0 = -log(20.); // threshold, could use here JB's thresholds
	double scale = 250.; // to scale the weight so to get ~same nbr of points as for E-weight 
	                     //  for the highest hit of ~0.1 GeV
	int nhit = int(scale*(w0+log(theHit->energy()/mip_)));
	if (nhit<0) nhit=1;
	for (int i=0; i<nhit; i++) principal_->AddRow(variables);	
	}	     
      }
    }
    
    alreadyfilled_ = true;
    principal_->MakePrincipals();
  
  }
  
  if (debug_) std::cout << "*** Principal component analysis (standalone) ****" << std::endl;
  if (debug_) std::cout << "shower average (x,y,z) = " << "(" << (*principal_->GetMeanValues())[0] << ", " <<
    (*principal_->GetMeanValues())[1] << ", " << (*principal_->GetMeanValues())[2] << ")" << std::endl;
  TMatrixD matrix = *principal_->GetEigenVectors();
  if (debug_) std::cout << "shower main axis (x,y,z) = " << "(" << matrix(0,0) << ", " <<
    matrix(1,0) << ", " << matrix(2,0) << ")" << std::endl;
  
  barycenter = GlobalPoint((*principal_->GetMeanValues())[0],(*principal_->GetMeanValues())[1],(*principal_->GetMeanValues())[2]);	 
  axis = GlobalVector(matrix(0,0),matrix(1,0),matrix(2,0));
   
  // resolve direction ambiguity
  if (axis.z()*barycenter.z()<0) {
    axis = GlobalVector(-matrix(0,0),-matrix(1,0),-matrix(2,0));
    if (debug_) std::cout << "PCA shower dir reverted " << axis << "eta " << axis.eta() << " phi " << axis.phi() << std::endl;
  }
	 
  return;

}

GlobalPoint PCAShowerAnalysis::showerBarycenter(const reco::CaloCluster* clus)
{

  if (!alreadyfilled_) {
  
    double variables[3] = {0.,0.,0.};
  
    for (unsigned int ih=0;ih<clus->hitsAndFractions().size();++ih) {
      const DetId & id_ = (clus->hitsAndFractions())[ih].first ;
      HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
      if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	GlobalPoint cellPos = geometry_->getPosition(HGCEEDetId(id_));
	variables[0] = cellPos.x(); variables[1] = cellPos.y(); variables[2] = cellPos.z();
	if (!segmented_) variables[2] = entryz_;
        if (!logweighting_) {
	// energy weighting
	for (int i=0; i<int(theHit->energy()/mip_); i++) principal_->AddRow(variables); 
	} else {
	// a log-weighting, energy not in fraction of total
	double w0 = -log(20.); // threshold, could yuse here JB's thresholds
	double scale = 250.; // to scale the weight so to get ~same nbr of points as for E-weight 
	                     //  for the highest hit of ~0.1 GeV
	int nhit = int(scale*(w0+log(theHit->energy()/mip_)));
	if (nhit<0) nhit=1;
	for (int i=0; i<nhit; i++) principal_->AddRow(variables);	
	}	     
      }
    }
    
    alreadyfilled_ = true;
    principal_->MakePrincipals();
  
  }
  
  if (debug_) std::cout << "*** Principal component analysis (standalone) ****" << std::endl;
  if (debug_) std::cout << "shower average (x,y,z) = " << "(" << (*principal_->GetMeanValues())[0] << ", " <<
    (*principal_->GetMeanValues())[1] << ", " << (*principal_->GetMeanValues())[2] << ")" << std::endl;
  
  GlobalPoint pcaShowerPos((*principal_->GetMeanValues())[0],(*principal_->GetMeanValues())[1],(*principal_->GetMeanValues())[2]);	 
   
  return pcaShowerPos;
  
}

GlobalVector PCAShowerAnalysis::showerAxis( const reco::CaloCluster* clus)
{
  
  if (!alreadyfilled_) {
  
    double variables[3] = {0.,0.,0.};

    for (unsigned int ih=0;ih<clus->hitsAndFractions().size();++ih) {
      const DetId & id_ = (clus->hitsAndFractions())[ih].first ;
      HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
      if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	GlobalPoint cellPos = geometry_->getPosition(HGCEEDetId(id_));
	variables[0] = cellPos.x(); variables[1] = cellPos.y(); variables[2] = cellPos.z();
	if (!segmented_) variables[2] = entryz_;
        if (!logweighting_) {
	// energy weighting
	for (int i=0; i<int(theHit->energy()/mip_); i++) principal_->AddRow(variables); 
	} else {
	// a log-weighting, energy not in fraction of total
	double w0 = -log(20.); // threshold, could yuse here JB's thresholds
	double scale = 250.; // to scale the weight so to get ~same nbr of points as for E-weight 
	                     //  for the highest hit of ~0.1 GeV
	int nhit = int(scale*(w0+log(theHit->energy()/mip_)));
	if (nhit<0) nhit=1;
	for (int i=0; i<nhit; i++) principal_->AddRow(variables);	
	}	     
      }
    }

    alreadyfilled_ = true;
    principal_->MakePrincipals();
  
  }
  
  if (debug_) std::cout << "*** Principal component analysis (standalone) ****" << std::endl;
  TMatrixD matrix = *principal_->GetEigenVectors();
  if (debug_) std::cout << "shower main axis (x,y,z) = " << "(" << matrix(0,0) << ", " <<
    matrix(1,0) << ", " << matrix(2,0) << ")" << std::endl;
  
  GlobalPoint pcaShowerPos((*principal_->GetMeanValues())[0],(*principal_->GetMeanValues())[1],(*principal_->GetMeanValues())[2]);	 
  GlobalVector pcaShowerDir(matrix(0,0),matrix(1,0),matrix(2,0));
   
  // resolve direction ambiguity
  if (pcaShowerDir.z()*pcaShowerPos.z()<0) {
    pcaShowerDir = GlobalVector(-matrix(0,0),-matrix(1,0),-matrix(2,0));
    if (debug_) std::cout << "PCA shower dir reverted " << pcaShowerDir << "eta " << pcaShowerDir.eta() << " phi " << pcaShowerDir.phi() << std::endl;
  }
	 
  return pcaShowerDir;
  
}

GlobalVector PCAShowerAnalysis::showerEigenValues(const reco::CaloCluster* clus)
{
  
    if (!alreadyfilled_) {
  
    double variables[3] = {0.,0.,0.};

    for (unsigned int ih=0;ih<clus->hitsAndFractions().size();++ih) {
      const DetId & id_ = (clus->hitsAndFractions())[ih].first ;
      HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
      if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	GlobalPoint cellPos = geometry_->getPosition(HGCEEDetId(id_));
	variables[0] = cellPos.x(); variables[1] = cellPos.y(); variables[2] = cellPos.z();
	if (!segmented_) variables[2] = entryz_;
        if (!logweighting_) {
	// energy weighting
	for (int i=0; i<int(theHit->energy()/mip_); i++) principal_->AddRow(variables); 
	} else {
	// a log-weighting, energy not in fraction of total
	double w0 = -log(20.); // threshold, could yuse here JB's thresholds
	double scale = 250.; // to scale the weight so to get ~same nbr of points as for E-weight 
	                     //  for the highest hit of ~0.1 GeV
	int nhit = int(scale*(w0+log(theHit->energy()/mip_)));
	if (nhit<0) nhit=1;
	for (int i=0; i<nhit; i++) principal_->AddRow(variables);	
	}	     
      }
    }

    alreadyfilled_ = true;
    principal_->MakePrincipals();

  }
  
  TVectorD eigenvalues = *principal_->GetEigenValues();
  if (debug_) std::cout << "*** Principal component analysis (standalone) ****" << std::endl;
  if (debug_) std::cout << "shower eigen values = " << "(" << eigenvalues(0) << ", " <<
    eigenvalues(1) << ", " << eigenvalues(2) << ")" << std::endl;
  	
  GlobalVector pcaShowerEigenValues(eigenvalues(0), eigenvalues(1), eigenvalues(2)); 	 
  return pcaShowerEigenValues; 
  
}

GlobalVector PCAShowerAnalysis::showerSigmas( const reco::CaloCluster* clus)
{
  
  if (!alreadyfilled_) {
  
    double variables[3] = {0.,0.,0.};

    for (unsigned int ih=0;ih<clus->hitsAndFractions().size();++ih) {
      const DetId & id_ = (clus->hitsAndFractions())[ih].first ;
      HGCRecHitCollection::const_iterator theHit = recHits_->find(id_);    
      if (id_.det()==DetId::Forward && id_.subdetId()==HGCEE) {
	GlobalPoint cellPos = geometry_->getPosition(HGCEEDetId(id_));
	variables[0] = cellPos.x(); variables[1] = cellPos.y(); variables[2] = cellPos.z();
	if (!segmented_) variables[2] = entryz_;
        if (!logweighting_) {
	// energy weighting
	for (int i=0; i<int(theHit->energy()/mip_); i++) principal_->AddRow(variables); 
	} else {
	// a log-weighting, energy not in fraction of total
	double w0 = -log(20.); // threshold, could yuse here JB's thresholds
	double scale = 250.; // to scale the weight so to get ~same nbr of points as for E-weight 
	                     //  for the highest hit of ~0.1 GeV
	int nhit = int(scale*(w0+log(theHit->energy()/mip_)));
	if (nhit<0) nhit=1;
	for (int i=0; i<nhit; i++) principal_->AddRow(variables);	
	}	     
      }
    }

    alreadyfilled_ = true;
    principal_->MakePrincipals();
  
  }
  
  TVectorD sigmas = *principal_->GetSigmas();
  if (debug_) std::cout << "*** Principal component analysis (standalone) ****" << std::endl;
  if (debug_) std::cout << "shower sigmas = " << "(" << sigmas(0) << ", " <<
    sigmas(1) << ", " << sigmas(2) << ")" << std::endl;
  	
  GlobalVector pcaShowerSigmas(sigmas(0), sigmas(1), sigmas(2)); 	 
  return pcaShowerSigmas; 
  
}
