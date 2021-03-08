#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "RecoEgamma/EgammaTools/interface/EgammaRegressionContainer.h"
#include "RecoEgamma/EgammaTools/interface/EpCombinationTool.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include <vdt/vdtMath.h>

class EGRegressionModifierLPV1 : public ModifyObjectValueBase {
public: 




  struct EleRegs {
    EleRegs(const edm::ParameterSet& iConfig);
    void setEventContent(const edm::EventSetup& iSetup);
    EgammaRegressionContainer ecalOnlyMean;
    EgammaRegressionContainer ecalOnlySigma;
    EpCombinationTool epComb;
  };

  struct PhoRegs {
    PhoRegs(const edm::ParameterSet& iConfig);
    void setEventContent(const edm::EventSetup& iSetup);
    EgammaRegressionContainer ecalOnlyMean;
    EgammaRegressionContainer ecalOnlySigma;    
  };

  EGRegressionModifierLPV1(const edm::ParameterSet& conf);
  ~EGRegressionModifierLPV1() override;
    
  void setEvent(const edm::Event&) final;
  void setEventContent(const edm::EventSetup&) final;
  void setConsumes(edm::ConsumesCollector&) final;
  
  void modifyObject(reco::GsfElectron&) const final;
  void modifyObject(reco::Photon&) const final;
  
  // just calls reco versions
  void modifyObject(pat::Electron&) const final; 
  void modifyObject(pat::Photon&) const final;

private:
  std::array<float,35> getRegData(const reco::GsfElectron& ele)const;
  std::array<float,35> getRegData(const reco::Photon& pho)const;
  void getSeedCrysCoord(const reco::CaloCluster& clus,int& iEtaOrX,int& iPhiOrY)const;

  std::unique_ptr<EleRegs> eleRegs_;
  std::unique_ptr<PhoRegs> phoRegs_;
  
  float rhoValue_;
  edm::InputTag rhoTag_;
  edm::EDGetTokenT<double> rhoToken_;

  bool useClosestToCentreSeedCrysDef_;
  float maxRawEnergyForLowPtEBSigma_;
  float maxRawEnergyForLowPtEESigma_;
  edm::ESHandle<CaloGeometry> caloGeomHandle_;

};

DEFINE_EDM_PLUGIN(ModifyObjectValueFactory,
          EGRegressionModifierLPV1,
          "EGRegressionModifierLPV1");

EGRegressionModifierLPV1::EGRegressionModifierLPV1(const edm::ParameterSet& conf) : 
  ModifyObjectValueBase(conf),
  rhoValue_(0.),
  rhoTag_(conf.getParameter<std::string>("rhoTag")),
  useClosestToCentreSeedCrysDef_(conf.getParameter<bool>("useClosestToCentreSeedCrysDef")),
  maxRawEnergyForLowPtEBSigma_(conf.getParameter<double>("maxRawEnergyForLowPtEBSigma")),
  maxRawEnergyForLowPtEESigma_(conf.getParameter<double>("maxRawEnergyForLowPtEESigma"))
{   
  if(conf.exists("eleRegs")) {
    eleRegs_ = std::make_unique<EleRegs>(conf.getParameter<edm::ParameterSet>("eleRegs"));
  }
  if(conf.exists("phoRegs")) {
    phoRegs_ = std::make_unique<PhoRegs>(conf.getParameter<edm::ParameterSet>("phoRegs"));
  }
}

EGRegressionModifierLPV1::~EGRegressionModifierLPV1(){}


void EGRegressionModifierLPV1::setConsumes(edm::ConsumesCollector& sumes) {

  rhoToken_ = sumes.consumes<double>(rhoTag_);

}



void EGRegressionModifierLPV1::setEvent(const edm::Event& evt)
{

  edm::Handle<double> rhoHandle;
  evt.getByToken(rhoToken_, rhoHandle);
  rhoValue_ = *rhoHandle;
}

void EGRegressionModifierLPV1::setEventContent(const edm::EventSetup& iSetup)
{
  if(eleRegs_) eleRegs_->setEventContent(iSetup);
  if(phoRegs_) phoRegs_->setEventContent(iSetup);
  if(useClosestToCentreSeedCrysDef_) iSetup.get<CaloGeometryRecord>().get(caloGeomHandle_);
}

void EGRegressionModifierLPV1::modifyObject(reco::GsfElectron& ele) const
{
  //check if we have specified an electron regression correction and
  //return the object unmodified if so
  if(!eleRegs_) return;

  const reco::SuperClusterRef& superClus = ele.superCluster();
  
  // skip HGCAL for now
  if( superClus->seed()->seed().det() == DetId::Forward ) return;

  // do not apply corrections in case of missing info (slimmed MiniAOD electrons)
  if(!superClus->clusters().isAvailable()) return; 

  //check if fbrem is filled as its needed for E/p combination so abort if its set to the default value 
  //this will be the case for <5 (or current cuts) for miniAOD electrons
  //  if(ele.fbrem()==reco::GsfElectron::ClassificationVariables::kDefaultValue) return;
  
  if(ele.fbrem()== -1.e30) return;
  
  auto regData = getRegData(ele);
  const float rawEnergy = superClus->rawEnergy(); 
  const float rawESEnergy = superClus->preshowerEnergy();
  //bug here, it should include the ES, kept for backwards compat
  const float rawEt = rawEnergy*superClus->position().rho()/superClus->position().r();
  const bool isSaturated = ele.nSaturatedXtals()!=0;

  const float ecalMean = eleRegs_->ecalOnlyMean(rawEt,ele.isEB(),isSaturated,regData.data());
  const float ecalMeanCorr = ecalMean > 0. ? ecalMean : 1.0;
  //as the sample is trained flat in pt, the regression's only source of very high energy
  //electrons is in the high endcap and therefore it gives a very poor resolution estimate
  //to any electrons with this energy, regardless of their actual eta
  //hence this lovely hack
  if(ele.isEB() && maxRawEnergyForLowPtEBSigma_>=0 && 
     eleRegs_->ecalOnlySigma.useLowEtBin(rawEt,isSaturated)){
    regData[0] = std::min(regData[0],maxRawEnergyForLowPtEBSigma_);
  }
  if(!ele.isEB() && maxRawEnergyForLowPtEESigma_>=0 && 
     eleRegs_->ecalOnlySigma.useLowEtBin(rawEt,isSaturated)){
    regData[0] = std::min(regData[0],maxRawEnergyForLowPtEESigma_);
  }
  const float ecalSigma = eleRegs_->ecalOnlySigma(rawEt,ele.isEB(),isSaturated,regData.data());

  const float corrEnergy = (rawEnergy + rawESEnergy)*ecalMeanCorr;
  const float corrEnergyErr = corrEnergy*ecalSigma;

  ele.setCorrectedEcalEnergy(corrEnergy);
  ele.setCorrectedEcalEnergyError(corrEnergyErr);

  std::pair<float,float> combEnergyAndErr =  eleRegs_->epComb.combine(ele);
  const math::XYZTLorentzVector newP4 = ele.p4()*combEnergyAndErr.first/ele.p4().t();
  ele.correctMomentum(newP4, ele.trackMomentumError(), combEnergyAndErr.second);
}

void EGRegressionModifierLPV1::modifyObject(pat::Electron& ele) const
{
  modifyObject(static_cast<reco::GsfElectron&>(ele));
}

void EGRegressionModifierLPV1::modifyObject(reco::Photon& pho) const
{
  //check if we have specified an photon regression correction and
  //return the object unmodified if so
  if(!phoRegs_) return;
  
  const reco::SuperClusterRef& superClus = pho.superCluster();

  // skip HGCAL for now
  if(superClus->seed()->seed().det() == DetId::Forward ) return;

  // do not apply corrections in case of missing info (happens for some slimmed MiniAOD photons)
  if(!superClus->clusters().isAvailable()) return; 

  auto regData = getRegData(pho); 

  const float rawEnergy = superClus->rawEnergy(); 
  const float rawESEnergy = superClus->preshowerEnergy();
  //bug here, it should include the ES, kept for backwards compat
  const float rawEt = rawEnergy*superClus->position().rho()/superClus->position().r();
  const bool isSaturated = pho.nSaturatedXtals();
  const float ecalMean = phoRegs_->ecalOnlyMean(rawEt,pho.isEB(),isSaturated,regData.data());
  const float ecalMeanCorr = ecalMean > 0. ? ecalMean : 1.0;

  //see the electrons for explaination of this lovely feature
  if(pho.isEB() && maxRawEnergyForLowPtEBSigma_>=0 && 
     phoRegs_->ecalOnlySigma.useLowEtBin(rawEt,isSaturated)){
    regData[0] = std::min(regData[0],maxRawEnergyForLowPtEBSigma_);
  }
  if(!pho.isEB() && maxRawEnergyForLowPtEESigma_>=0 && 
     phoRegs_->ecalOnlySigma.useLowEtBin(rawEt,isSaturated)){
    regData[0] = std::min(regData[0],maxRawEnergyForLowPtEESigma_);
  }
  const float ecalSigma = phoRegs_->ecalOnlySigma(rawEt,pho.isEB(),isSaturated,regData.data());
    
  const double corrEnergy = (rawEnergy + rawESEnergy)*ecalMeanCorr;
  const double corrEnergyErr = corrEnergy*ecalSigma;

  pho.setCorrectedEnergy(reco::Photon::P4type::regression2, corrEnergy, corrEnergyErr, true);     
}

void EGRegressionModifierLPV1::modifyObject(pat::Photon& pho) const {
  modifyObject(static_cast<reco::Photon&>(pho));
}

std::array<float,35> EGRegressionModifierLPV1::getRegData(const reco::GsfElectron& ele)const
{
  std::array<float, 35> data;  

  const reco::SuperClusterRef& superClus = ele.superCluster();
  const edm::Ptr<reco::CaloCluster>& seedClus = superClus->seed();
  
  const bool isEB = ele.isEB();  
  const double rawEnergy = superClus->rawEnergy(); 
  const double rawESEnergy = superClus->preshowerEnergy();
  const int numberOfClusters =  superClus->clusters().size();
  const auto& ssFull5x5 = ele.full5x5_showerShape();

  float e5x5Inverse = ssFull5x5.e5x5 != 0. ? vdt::fast_inv(ssFull5x5.e5x5) : 0.;

  data[0]  = rawEnergy;
  data[1]  = superClus->etaWidth();
  data[2]  = superClus->phiWidth(); 
  data[3]  = superClus->seed()->energy()/rawEnergy;
  data[4]  = ssFull5x5.e5x5/rawEnergy;
  data[5]  = ele.hcalOverEcalBc();
  data[6]  = rhoValue_;
  data[7]  = seedClus->eta() - superClus->position().Eta();
  data[8]  = reco::deltaPhi( seedClus->phi(),superClus->position().Phi());
  data[9]  = ssFull5x5.r9;
  data[10]  = ssFull5x5.sigmaIetaIeta;
  data[11]  = ssFull5x5.sigmaIetaIphi;
  data[12]  = ssFull5x5.sigmaIphiIphi;
  data[13]  = ssFull5x5.eMax*e5x5Inverse;
  data[14]  = ssFull5x5.e2nd*e5x5Inverse;
  data[15]  = ssFull5x5.eTop*e5x5Inverse;
  data[16]  = ssFull5x5.eBottom*e5x5Inverse;
  data[17]  = ssFull5x5.eLeft*e5x5Inverse;
  data[18]  = ssFull5x5.eRight*e5x5Inverse;
  data[19]  = ssFull5x5.e2x5Max*e5x5Inverse;
  data[20]  = ssFull5x5.e2x5Left*e5x5Inverse;
  data[21]  = ssFull5x5.e2x5Right*e5x5Inverse;
  data[22]  = ssFull5x5.e2x5Top*e5x5Inverse;
  data[23]  = ssFull5x5.e2x5Bottom*e5x5Inverse;
  data[24]  = ele.nSaturatedXtals();
  data[25]  = std::max(0,numberOfClusters);

  // extra new variables: nxtals, clus1.nXtals, clus1.energy

  float sc_nxtals=0;
  float maxEne1=0;
  float clus1_nxtals=0;
  
  try{
    if(superClus->clustersSize()>0&& superClus->clustersBegin()!=superClus->clustersEnd()){
      for(auto& cluster : superClus->clusters()) {
	if (cluster->energy() > maxEne1){
	  maxEne1=cluster->energy();
	  clus1_nxtals=(float)cluster->size();
	}
	sc_nxtals=sc_nxtals+(float)cluster->size(); // count nxtals 
      }
    }
  }catch(...){
    // sometimes there are no clusters  
  }
  


  
  if(isEB){
    int iEta,iPhi;
    getSeedCrysCoord(*seedClus,iEta,iPhi);
    int signIEta = iEta > 0 ? +1 : -1;
    data[26] = iEta;
    data[27] = iPhi;
    data[28] = (iEta-signIEta)%5;
    data[29] = (iPhi-1)%2; 
    const int iEtaCorr = iEta - (iEta > 0 ? +1 : -1);
    const int iEtaCorr26 = iEta - (iEta > 0 ? +26 : -26);
    data[30] = std::abs(iEta)<=25 ? iEtaCorr%20 : iEtaCorr26%20;
    data[31] = (iPhi-1)%20;
    data[32] = sc_nxtals;
    data[33] = clus1_nxtals;
    data[34] = maxEne1/rawEnergy;
  }else{
    int iX,iY;
    getSeedCrysCoord(*seedClus,iX,iY);
    data[26] = iX;
    data[27] = iY;
    data[28] = rawESEnergy/rawEnergy;
    data[29] = sc_nxtals;
    data[30] = clus1_nxtals;
    data[31] = maxEne1/rawEnergy;
  }

  return data;
}

std::array<float,35> EGRegressionModifierLPV1::getRegData(const reco::Photon& pho)const
{
  std::array<float, 35> data;  

  const reco::SuperClusterRef& superClus = pho.superCluster();
  const edm::Ptr<reco::CaloCluster>& seedClus = superClus->seed();
  
  const bool isEB = pho.isEB();  
  const double rawEnergy = superClus->rawEnergy(); 
  const double rawESEnergy = superClus->preshowerEnergy(); 
  const int numberOfClusters =  superClus->clusters().size();
  const auto& ssFull5x5 = pho.full5x5_showerShapeVariables();

  float e5x5Inverse = ssFull5x5.e5x5 != 0. ? vdt::fast_inv(ssFull5x5.e5x5) : 0.;

  data[0]  = rawEnergy;
  data[1]  = superClus->etaWidth();
  data[2]  = superClus->phiWidth(); 
  data[3]  = superClus->seed()->energy()/rawEnergy;
  data[4]  = ssFull5x5.e5x5/rawEnergy;
  //interestingly enough this differs from electrons where it uses cone based
  //naively Sam would have thought using cone based is even worse than tower based
  data[5]  = pho.hadronicOverEm();
  data[6]  = rhoValue_;
  data[7]  = seedClus->eta() - superClus->position().Eta();
  data[8]  = reco::deltaPhi( seedClus->phi(),superClus->position().Phi());
  data[9]  = pho.full5x5_r9();
  data[10]  = ssFull5x5.sigmaIetaIeta;
  //interestingly sigmaIEtaIPhi differs in defination here from 
  //electron & sc definations of sigmaIEtaIPhi
  data[11]  = ssFull5x5.sigmaIetaIphi;
  data[12]  = ssFull5x5.sigmaIphiIphi;
  data[13]  = ssFull5x5.maxEnergyXtal*e5x5Inverse;
  data[14]  = ssFull5x5.e2nd*e5x5Inverse;
  data[15]  = ssFull5x5.eTop*e5x5Inverse;
  data[16]  = ssFull5x5.eBottom*e5x5Inverse;
  data[17]  = ssFull5x5.eLeft*e5x5Inverse;
  data[18]  = ssFull5x5.eRight*e5x5Inverse;
  data[19]  = ssFull5x5.e2x5Max*e5x5Inverse;
  data[20]  = ssFull5x5.e2x5Left*e5x5Inverse;
  data[21]  = ssFull5x5.e2x5Right*e5x5Inverse;
  data[22]  = ssFull5x5.e2x5Top*e5x5Inverse;
  data[23]  = ssFull5x5.e2x5Bottom*e5x5Inverse;
  data[24]  = pho.nSaturatedXtals();
  data[25]  = std::max(0,numberOfClusters);
      
  if(isEB){
    int iEta,iPhi;
    getSeedCrysCoord(*seedClus,iEta,iPhi);
    data[26] = iEta;
    data[27] = iPhi;
    int signIEta = iEta > 0 ? +1 : -1;
    data[28] = (iEta-signIEta)%5;
    data[29] = (iPhi-1)%2;
    const int iEtaCorr = iEta - (iEta > 0 ? +1 : -1);
    const int iEtaCorr26 = iEta - (iEta > 0 ? +26 : -26);
    data[30] = std::abs(iEta)<=25 ? iEtaCorr%20 : iEtaCorr26%20;
    data[31] = (iPhi-1)%20;
  }else{
    int iX,iY;
    getSeedCrysCoord(*seedClus,iX,iY);
    data[26] = iX;
    data[27] = iY;
    data[28] = rawESEnergy/rawEnergy;
  }

  return data;
}

void EGRegressionModifierLPV1::getSeedCrysCoord(const reco::CaloCluster& clus,int& iEtaOrX,int& iPhiOrY)const
{
  iEtaOrX = 0;
  iPhiOrY = 0;

  const bool isEB = clus.seed().subdetId()==EcalBarrel;

    if(isEB){
      const EBDetId ebId(clus.seed());
      iEtaOrX = ebId.ieta();
      iPhiOrY = ebId.iphi();
    }else{
      const EEDetId eeId(clus.seed());
      iEtaOrX = eeId.ix();
      iPhiOrY = eeId.iy();
    }

}

EGRegressionModifierLPV1::EleRegs::EleRegs(const edm::ParameterSet& iConfig):
  ecalOnlyMean(iConfig.getParameter<edm::ParameterSet>("ecalOnlyMean")),
  ecalOnlySigma(iConfig.getParameter<edm::ParameterSet>("ecalOnlySigma")),
  epComb(iConfig.getParameter<edm::ParameterSet>("epComb"))
{
  
}

void EGRegressionModifierLPV1::EleRegs::setEventContent(const edm::EventSetup& iSetup)
{
  ecalOnlyMean.setEventContent(iSetup);
  ecalOnlySigma.setEventContent(iSetup);
  epComb.setEventContent(iSetup);
}

EGRegressionModifierLPV1::PhoRegs::PhoRegs(const edm::ParameterSet& iConfig):
  ecalOnlyMean(iConfig.getParameter<edm::ParameterSet>("ecalOnlyMean")),
  ecalOnlySigma(iConfig.getParameter<edm::ParameterSet>("ecalOnlySigma"))
{
  
}
  
void EGRegressionModifierLPV1::PhoRegs::setEventContent(const edm::EventSetup& iSetup)
{
  ecalOnlyMean.setEventContent(iSetup);
  ecalOnlySigma.setEventContent(iSetup);
}
