#ifndef RecoEgamma_EgammaElectronProducers_LowPtGsfElectronIDExtraHeavyObjectCache_h
#define RecoEgamma_EgammaElectronProducers_LowPtGsfElectronIDExtraHeavyObjectCache_h

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>

namespace reco { 
  class BeamSpot;
  class PreId; 
  typedef edm::Ptr<GsfElectron> GsfElectronPtr;
  typedef edm::Ref< edm::View<GsfElectron> > LowPtGsfElectronRef;
}

namespace lowptgsfeleidextra {

  class Features {
  public:
    float eid_rho = -999.;                  
    float eid_sc_eta = -999.;               
    float eid_shape_full5x5_r9 = -999.;     
    float eid_sc_etaWidth = -999.;          
    float eid_sc_phiWidth = -999.;          
    float eid_shape_full5x5_HoverE = -999.; 
    float eid_trk_nhits = -999.;            
    float eid_trk_chi2red = -999.;          
    float eid_gsf_chi2red = -999.;          
    float eid_brem_frac = -999.;            
    float eid_gsf_nhits = -999.;            
    float eid_match_SC_EoverP = -999.;      
    float eid_match_eclu_EoverP = -999.;    
    float eid_match_SC_dEta   = -999.;      
    float eid_match_SC_dPhi   = -999.;      
    float eid_match_seed_dEta = -999.;      
    float eid_sc_E = -999.;                 
    float eid_trk_p = -999.;                
    float gsf_mode_p = -999.;               
    float core_shFracHits = -999.;          
    float gsf_bdtout1 = -999.;              
    float gsf_dr = -999.;                   
    float trk_dr = -999.; 
    float sc_Nclus = -999.;
    float sc_clus1_nxtal  = -999.;  
    float sc_clus1_dphi = -999.;  
    float sc_clus2_dphi = -999.;  
    float sc_clus1_deta = -999.;  
    float sc_clus2_deta = -999.;  
    float sc_clus1_E = -999.;  
    float sc_clus2_E = -999.;  
    float sc_clus1_E_ov_p = -999.;  
    float sc_clus2_E_ov_p = -999.;    

  public:
    std::vector<float> get();
    void set( const reco::LowPtGsfElectronRef& ele, double rho, float unbiased = -10. );
    void set( const reco::GsfElectronPtr& ele, double rho, float unbiased = -10. );
  };

  class HeavyObjectCache {

  public:

    HeavyObjectCache( const edm::ParameterSet& );

    std::vector<std::string> modelNames() const { return names_; }
    
    double eval( const std::string& name, const reco::LowPtGsfElectronRef&, double rho, float unbiased = -10. ) const;
    double eval( const std::string& name, const reco::GsfElectronPtr&, double rho, float unbiased = -10. ) const;
    
  private:

    std::vector<std::string> names_;
    std::vector< std::unique_ptr<const GBRForest> > models_;
    std::vector<double> thresholds_;

  };
}

#endif // RecoEgamma_EgammaElectronProducers_LowPtGsfElectronIDExtraHeavyObjectCache_h
