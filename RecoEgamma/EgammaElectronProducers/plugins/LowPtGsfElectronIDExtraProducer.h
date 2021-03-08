#ifndef RecoEgamma_EgammaElectronProducers_LowPtGsfElectronIDExtraProducer_h
#define RecoEgamma_EgammaElectronProducers_LowPtGsfElectronIDExtraProducer_h

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronIDExtraHeavyObjectCache.h"
#include <string>
#include <vector>

class LowPtGsfElectronIDExtraProducer final : 
public edm::stream::EDProducer< edm::GlobalCache<lowptgsfeleidextra::HeavyObjectCache> > {
  
 public:
  
  explicit LowPtGsfElectronIDExtraProducer( const edm::ParameterSet&,
					    const lowptgsfeleidextra::HeavyObjectCache* );
  
  ~LowPtGsfElectronIDExtraProducer() override;
  
  static std::unique_ptr<lowptgsfeleidextra::HeavyObjectCache> 
    initializeGlobalCache( const edm::ParameterSet& conf ) {
    return std::make_unique<lowptgsfeleidextra::HeavyObjectCache>(conf);  
  }
  
  static void globalEndJob( lowptgsfeleidextra::HeavyObjectCache const* ) {}

  void produce( edm::Event&, const edm::EventSetup& ) override;
  
  static void fillDescriptions( edm::ConfigurationDescriptions& );

 private:
  
  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > gsfElectrons_; 
  const edm::EDGetTokenT<double> rho_;                                  
  const edm::EDGetTokenT< edm::ValueMap<float> > unbiased_;
  const std::vector<std::string> names_;
  const bool passThrough_;
  const double minPtThreshold_;
  const double maxPtThreshold_;

};

#endif // RecoEgamma_EgammaElectronProducers_LowPtGsfElectronIDExtraProducer_h

