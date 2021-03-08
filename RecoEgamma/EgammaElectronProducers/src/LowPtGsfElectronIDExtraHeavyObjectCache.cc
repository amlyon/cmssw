#include "CommonTools/MVAUtils/interface/GBRForestTools.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronIDExtraHeavyObjectCache.h"
#include <string>
#include <TVector3.h>

namespace lowptgsfeleidextra {

  ////////////////////////////////////////////////////////////////////////////////
  //
  std::vector<float> Features::get() {
    std::vector<float> output = { 
      eid_rho, 
      eid_sc_eta,           
      eid_shape_full5x5_r9,     
      eid_sc_etaWidth,      
      eid_sc_phiWidth,      
      eid_shape_full5x5_HoverE, 
      eid_trk_nhits,        
      eid_trk_chi2red,      
      eid_gsf_chi2red,      
      eid_brem_frac,        
      eid_gsf_nhits,        
      eid_match_SC_EoverP,  
      eid_match_eclu_EoverP,
      eid_match_SC_dEta,    
      eid_match_SC_dPhi,    
      eid_match_seed_dEta,  
      eid_sc_E,
      eid_trk_p,             
      gsf_mode_p,               
      core_shFracHits,          
      gsf_bdtout1,  
      gsf_dr,                   
      trk_dr,                   
      sc_Nclus,
      sc_clus1_nxtal,
      sc_clus1_dphi,
      sc_clus2_dphi,
      sc_clus1_deta,
      sc_clus2_deta,
      sc_clus1_E,
      sc_clus2_E,
      sc_clus1_E_ov_p,
      sc_clus2_E_ov_p
    };
    return output;
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  void Features::set( const reco::LowPtGsfElectronRef& ele, double rho, float unbiased ) {
    set(edm::refToPtr(ele),rho,unbiased);
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  void Features::set( const reco::GsfElectronPtr& ele, double rho, float unbiased ) {

    // KF tracks                                       
    if ( ele->core().isNonnull() ) {     
      reco::TrackRef trk = ele->closestCtfTrackRef(); 
      if ( trk.isNonnull() ) {
	eid_trk_p = (float)trk->p();                                                   
	eid_trk_nhits = (float)trk->found();                
	eid_trk_chi2red = (float)trk->normalizedChi2();     
	TVector3 trkTV3(0,0,0);
	trkTV3.SetPtEtaPhi(trk->pt(), trk->eta(), trk->phi());  
	TVector3 eleTV3(0,0,0);
	eleTV3.SetPtEtaPhi(ele->pt(), ele->eta(), ele->phi());
	trk_dr = eleTV3.DeltaR(trkTV3);  
      }
    }

    // GSF tracks                                     
    if ( ele->core().isNonnull() ) {
      reco::GsfTrackRef gsf = ele->core()->gsfTrack();
      if ( gsf.isNonnull() ) {
	gsf_mode_p = gsf->pMode();       
	eid_gsf_nhits = (float)gsf->found();       
	eid_gsf_chi2red = gsf->normalizedChi2();       
	TVector3 gsfTV3(0,0,0);
	gsfTV3.SetPtEtaPhi(gsf->ptMode(), gsf->etaMode(), gsf->phiMode());  
	TVector3 eleTV3(0,0,0);
	eleTV3.SetPtEtaPhi(ele->pt(), ele->eta(), ele->phi());
	gsf_dr = eleTV3.DeltaR(gsfTV3);  
      }
    }

    // Super clusters                        
    if ( ele->core().isNonnull() ) {
      reco::SuperClusterRef sc = ele->core()->superCluster();
      if ( sc.isNonnull() ) {
	eid_sc_E = sc->energy();                 
	eid_sc_eta  = sc->eta();                 
	eid_sc_etaWidth = sc->etaWidth();        
	eid_sc_phiWidth = sc->phiWidth();        
	sc_Nclus = sc->clustersSize();          
      }
    }

    // Track-cluster matching              
    if ( ele.isNonnull() ) {
      eid_match_seed_dEta = ele->deltaEtaSeedClusterTrackAtCalo();
      eid_match_eclu_EoverP = (1./ele->ecalEnergy()) - (1./ele->p());
      eid_match_SC_EoverP = ele->eSuperClusterOverP();
      eid_match_SC_dEta = ele->deltaEtaSuperClusterTrackAtVtx();
      eid_match_SC_dPhi = ele->deltaPhiSuperClusterTrackAtVtx();
    }      

    // Shower shape vars        
    if ( ele.isNonnull() ) {
      eid_shape_full5x5_HoverE = ele->full5x5_hcalOverEcal();  
      eid_shape_full5x5_r9     = ele->full5x5_r9();            
    }

    // Misc
    eid_rho = rho;           

    if ( ele.isNonnull() ) {                          
      eid_brem_frac = ele->fbrem();
      core_shFracHits = ele->shFracInnerHits();
    }

    // Unbiased BDT from ElectronSeed   
    gsf_bdtout1 = unbiased;

    // Clusters   
    if ( ele->core().isNonnull() ) {
      reco::GsfTrackRef gsf = ele->core()->gsfTrack();
      if ( gsf.isNonnull() ) {
	reco::SuperClusterRef sc = ele->core()->superCluster();
	if ( sc.isNonnull() ) {	

	  // Propagate electron track to ECAL surface 
	  double mass_ = 0.000511*0.000511; 
	  float p2=pow( gsf->p() ,2 );
	  float energy = sqrt(mass_ + p2);
	  XYZTLorentzVector mom = XYZTLorentzVector(gsf->px(), gsf->py(), gsf->pz(), energy);
	  XYZTLorentzVector pos = XYZTLorentzVector(gsf->vx(), gsf->vy(), gsf->vz(), 0.);
	  float field_z=3.8;
	  BaseParticlePropagator mypart(RawParticle(mom,pos), 0, 0, field_z);
	  mypart.setCharge(gsf->charge());
	  mypart.propagateToEcalEntrance(true);    // true only first half loop , false more than one loop 
	  bool reach_ECAL = mypart.getSuccess();   // 0 does not reach ECAL, 1 yes barrel, 2 yes endcaps    
	  GlobalPoint ecal_pos(mypart.x(), mypart.y(), mypart.z());    // ECAL entry point for track  
	  
	  // Iterate through ECAL clusters and sort in energy
	  int clusNum=0;
	  float maxEne1=-1;
	  float maxEne2=-1;
	  int i1=-1;
	  int i2=-1;
	  try{
	    if(sc->clustersSize()>0 && sc->clustersBegin()!=sc->clustersEnd()){
	      for(auto& cluster : sc->clusters()) {
		if (cluster->energy() > maxEne1){
		  maxEne1=cluster->energy();
		  i1=clusNum;
		}
		clusNum++;
	      }
	      if(sc->clustersSize()>1){
		clusNum=0;
		for(auto& cluster : sc->clusters()) {
		  if (clusNum!=i1) {
		    if (cluster->energy() > maxEne2){
		      maxEne2=cluster->energy();
		      i2=clusNum;
		    }
		  }
		  clusNum++;
		}
	      }
	    } // loop over clusters
	  } catch(...) {
	    std::cout<<"exception caught clusNum="<<clusNum<<" clus size"<<sc->clustersSize()<<" energy="<< sc->energy()<<std::endl;
	  }

	  // Initializations
	  sc_clus1_nxtal  = -999;
	  sc_clus1_dphi   = -999.;
	  sc_clus2_dphi   = -999.;
	  sc_clus1_deta   = -999.;
	  sc_clus2_deta   = -999.;
	  sc_clus1_E      = -999.;
	  sc_clus2_E      = -999.;
	  sc_clus1_E_ov_p = -999.;
	  sc_clus2_E_ov_p = -999.;

	  // track-clusters match
	  clusNum=0;
	  try { 
	    if(sc->clustersSize()>0&& sc->clustersBegin()!=sc->clustersEnd()){
	      for(auto& cluster : sc->clusters()) {
		double pi_=3.1415926535;
		float deta = std::fabs(ecal_pos.eta()-cluster->eta()) ;
		float dphi = std::fabs(ecal_pos.phi()-cluster->phi());
		if (dphi > pi_)  dphi -= 2 * pi_;
		if (ecal_pos.phi()-cluster->phi()<0) dphi=-dphi;
		if (ecal_pos.eta()-cluster->eta()<0) deta=-deta;

		if (clusNum==i1) {
		  sc_clus1_E = cluster->energy();
		  if(gsf->pMode()>0) sc_clus1_E_ov_p = cluster->energy()/gsf->pMode();
		  sc_clus1_nxtal = (int)cluster->size();
		  if (reach_ECAL>0){
		    sc_clus1_deta = deta;
		    sc_clus1_dphi = dphi;
		  } 
		} else if (clusNum==i2) {   
		  sc_clus2_E    = cluster->energy();
		  if(gsf->pMode()>0) sc_clus2_E_ov_p = cluster->energy()/gsf->pMode();
		  if (reach_ECAL>0){
		    sc_clus2_deta = deta;
		    sc_clus2_dphi = dphi;
		  } 
		}
		clusNum++;
	      }
	    }
	  } catch(...) {
	    std::cout<<"caught an exception"<<std::endl;                            
	  }
	}
      }
    }  // clusters

    // Out-of-range
    if (eid_rho<0)   eid_rho=0;
    if (eid_rho>100) eid_rho=100;
    if (eid_sc_eta<-5) eid_sc_eta=-5;
    if (eid_sc_eta>5)  eid_sc_eta=5;
    if (eid_shape_full5x5_r9<0) eid_shape_full5x5_r9=0;
    if (eid_shape_full5x5_r9>2) eid_shape_full5x5_r9=2;
    if (eid_sc_etaWidth<0)    eid_sc_etaWidth=0;
    if (eid_sc_etaWidth>3.14) eid_sc_etaWidth=3.14;
    if (eid_sc_phiWidth<0)    eid_sc_phiWidth=0;
    if (eid_sc_phiWidth>3.14) eid_sc_phiWidth=3.14;
    if (eid_shape_full5x5_HoverE<0) eid_shape_full5x5_HoverE=0;
    if (eid_shape_full5x5_HoverE>50) eid_shape_full5x5_HoverE=50;
    if (eid_trk_nhits<-1) eid_trk_nhits=-1;
    if (eid_trk_nhits>50) eid_trk_nhits=50;
    if (eid_trk_chi2red<-1) eid_trk_chi2red=-1;
    if (eid_trk_chi2red>50) eid_trk_chi2red=50;
    if (eid_gsf_chi2red<-1)  eid_gsf_chi2red=-1;
    if (eid_gsf_chi2red>100) eid_gsf_chi2red=100;
    if (eid_brem_frac<0)  eid_brem_frac=-1;
    if (eid_brem_frac>1)  eid_brem_frac=1;
    if (eid_gsf_nhits<-1) eid_gsf_nhits=-1;
    if (eid_gsf_nhits>50) eid_gsf_nhits=50;
    if (eid_match_SC_EoverP<0)   eid_match_SC_EoverP=0;
    if (eid_match_SC_EoverP>100) eid_match_SC_EoverP=100;
    if (eid_match_eclu_EoverP<-1.) eid_match_eclu_EoverP=-1;
    if (eid_match_eclu_EoverP>1.)  eid_match_eclu_EoverP=1.;
    if (eid_match_SC_dEta<-10) eid_match_SC_dEta=-10;
    if (eid_match_SC_dEta>10)  eid_match_SC_dEta=10;
    if (eid_match_SC_dPhi<-3.14) eid_match_SC_dPhi=-3.14;
    if (eid_match_SC_dPhi>3.14)  eid_match_SC_dPhi=3.14;
    if (eid_match_seed_dEta<-10) eid_match_seed_dEta=-10;
    if (eid_match_seed_dEta>10)  eid_match_seed_dEta=10;    
    if (eid_sc_E<0)    eid_sc_E=0;
    if (eid_sc_E>1000) eid_sc_E=1000;
    if (eid_trk_p<-1)   eid_trk_p=-1;
    if (eid_trk_p>1000) eid_trk_p=1000;
    if (gsf_mode_p<0)    gsf_mode_p=0;
    if (gsf_mode_p>1000) gsf_mode_p=1000;
    if (core_shFracHits<0) core_shFracHits=0;
    if (core_shFracHits>1) core_shFracHits=1;
    if (gsf_bdtout1<-20) gsf_bdtout1=-20;
    if (gsf_bdtout1>20)  gsf_bdtout1=20;
    if (gsf_dr<0) gsf_dr=5;
    if (gsf_dr>5) gsf_dr=5;
    if (trk_dr<0) trk_dr=5;
    if (trk_dr>5) trk_dr=5;
    if (sc_Nclus<0)  sc_Nclus=0;
    if (sc_Nclus>20) sc_Nclus=20;    
    if (sc_clus1_nxtal<0)   sc_clus1_nxtal=0;
    if (sc_clus1_nxtal>100) sc_clus1_nxtal=100;
    if (sc_clus1_dphi<-3.14) sc_clus1_dphi=-5;
    if (sc_clus1_dphi>3.14)  sc_clus1_dphi=5;
    if (sc_clus2_dphi<-3.14) sc_clus2_dphi=-5;
    if (sc_clus2_dphi>3.14)  sc_clus2_dphi=5;
    if (sc_clus1_deta<-5) sc_clus1_deta=-5;
    if (sc_clus1_deta>5)  sc_clus1_deta=5;
    if (sc_clus2_deta<-5) sc_clus2_deta=-5;
    if (sc_clus2_deta>5)  sc_clus2_deta=5;
    if (sc_clus1_E<0)    sc_clus1_E=0;
    if (sc_clus1_E>1000) sc_clus1_E=1000;
    if (sc_clus2_E<0)    sc_clus2_E=0;
    if (sc_clus2_E>1000) sc_clus2_E=1000;
    if (sc_clus1_E_ov_p<0) sc_clus1_E_ov_p=-1;
    if (sc_clus2_E_ov_p<0) sc_clus2_E_ov_p=-1;
  };

  ////////////////////////////////////////////////////////////////////////////////
  //
  HeavyObjectCache::HeavyObjectCache( const edm::ParameterSet& conf ) {

    for ( auto& name : conf.getParameter< std::vector<std::string> >("ModelNames") ) 
      {
	names_.push_back(name);
      }
    for ( auto& weights : conf.getParameter< std::vector<std::string> >("ModelWeights") ) 
      {
	models_.push_back(createGBRForest(edm::FileInPath(weights)));
      }
    for ( auto& thresh : conf.getParameter< std::vector<double> >("ModelThresholds") )
      {
	thresholds_.push_back(thresh);
      }
    if ( names_.size() != models_.size() ) {
      throw cms::Exception("Incorrect configuration")
	<< "'ModelNames' size (" << names_.size()
	<< ") != 'ModelWeights' size (" << models_.size()
	<< ").\n";
    }
    if ( models_.size() != thresholds_.size() ) {
      throw cms::Exception("Incorrect configuration")
	<< "'ModelWeights' size (" << models_.size()
	<< ") != 'ModelThresholds' size (" << thresholds_.size()
	<< ").\n";
    }
    
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  double HeavyObjectCache::eval( const std::string& name,
				 const reco::LowPtGsfElectronRef& ele,
				 double rho,
				 float unbiased ) const
  {
    return eval(name,edm::refToPtr(ele),rho,unbiased);
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  double HeavyObjectCache::eval( const std::string& name,
				 const reco::GsfElectronPtr& ele,
				 double rho,
				 float unbiased ) const
  {
    std::vector<std::string>::const_iterator iter = std::find( names_.begin(), 
							       names_.end(), 
							       name );
    if ( iter != names_.end() ) {
      int index = std::distance(names_.begin(),iter);
      Features features;
      features.set(ele,rho,unbiased);
      std::vector<float> inputs = features.get();
      return models_.at(index)->GetResponse( inputs.data() );
    } else {
      throw cms::Exception("Unknown model name")
	<< "'Name given: '" << name
	<< "'. Check against configuration file.\n";
    }
    return 0.;
  }

} // namespace lowptgsfeleidextra
