#include "TrackingRecHitProducer.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "FastSimulation/TrackingRecHitProducer/interface/TrackerDetIdSelector.h"
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2DCollection.h"



#include <map>
#include <memory>




void insertRecHits(
    std::auto_ptr<TrackerGSRecHitCollection> targetCollection,
    TrackingRecHitProductPtr product
)
{
    if (product)
    {
        std::vector<SiTrackerGSRecHit2D>& producedRecHits = product->getRecHits();
        edm::OwnVector<SiTrackerGSRecHit2D> ownedRecHits;
        //TODO: this is not very performant; fix!
        for (unsigned int irecHit=0; irecHit<producedRecHits.size(); ++irecHit)
        {
            ownedRecHits.push_back(producedRecHits[irecHit]);
        }
        targetCollection->put(product->getDetId(),ownedRecHits.begin(),ownedRecHits.end());
    }
}

void insertMatchedRecHits(
    std::auto_ptr<TrackerGSMatchedRecHitCollection> targetCollection,
    TrackingRecHitProductPtr product
)
{
    if (product)
    {
        std::vector<SiTrackerGSMatchedRecHit2D>& producedRecHits = product->getMatchedRecHits();
        edm::OwnVector<SiTrackerGSMatchedRecHit2D> ownedRecHits;
        //TODO: this is not very performant; fix!
        for (unsigned int irecHit=0; irecHit<producedRecHits.size(); ++irecHit)
        {
            ownedRecHits.push_back(producedRecHits[irecHit]);
        }
        targetCollection->put(product->getDetId(),ownedRecHits.begin(),ownedRecHits.end());
    }
}


TrackingRecHitProducer::TrackingRecHitProducer(const edm::ParameterSet& config):
    _trackerGeometry(nullptr),
    _trackerTopology(nullptr)
{
    edm::ConsumesCollector consumeCollector = consumesCollector();
    const edm::ParameterSet& pluginConfigs = config.getParameter<edm::ParameterSet>("plugins");
    const std::vector<std::string> psetNames = pluginConfigs.getParameterNamesForType<edm::ParameterSet>();

    for (unsigned int iplugin = 0; iplugin<psetNames.size(); ++iplugin)
    {
        const edm::ParameterSet& pluginConfig = pluginConfigs.getParameter<edm::ParameterSet>(psetNames[iplugin]);
        const std::string pluginName = pluginConfig.getParameter<std::string>("type");
        TrackingRecHitAlgorithm* recHitAlgorithm = TrackingRecHitAlgorithmFactory::get()->tryToCreate(pluginName,psetNames[iplugin],pluginConfig,consumeCollector);
        if (recHitAlgorithm)
        {
            _recHitAlgorithms.push_back(recHitAlgorithm);
        }
        else
        {
            edm::LogWarning("TrackingRecHitAlgorithm plugin not found: ") << "plugin name = "<<pluginName<<"\nconfiguration=\n"<<pluginConfig.dump();
        }
    }

    edm::InputTag simHitTag = config.getParameter<edm::InputTag>("simHits");
    _simHitToken = consumes<std::vector<PSimHit>>(simHitTag);

    produces<TrackerGSRecHitCollection>("TrackerGSRecHits");
    produces<TrackerGSMatchedRecHitCollection>("TrackerGSMatchedRecHits");

}

void TrackingRecHitProducer::beginRun(edm::Run const&, const edm::EventSetup& eventSetup)
{
    edm::ESHandle<TrackerGeometry> trackerGeometryHandle;
    edm::ESHandle<TrackerTopology> trackerTopologyHandle;
    eventSetup.get<TrackerDigiGeometryRecord>().get(trackerGeometryHandle);
    eventSetup.get<IdealGeometryRecord>().get(trackerTopologyHandle);
    _trackerGeometry = trackerGeometryHandle.product();
    _trackerTopology = trackerTopologyHandle.product();
    
    for (TrackingRecHitAlgorithm* algo: _recHitAlgorithms)
    {
        algo->setupTrackerTopology(_trackerTopology);
    }
    
    
    //build pipes for all detIds
    const std::vector<DetId>& detIds = _trackerGeometry->detIds();
    
    for (const DetId& detId: detIds)
    {
        TrackerDetIdSelector selector(detId,*_trackerTopology);
       
        TrackingRecHitPipe pipe;
        for (TrackingRecHitAlgorithm* algo: _recHitAlgorithms)
        {
            if (selector.passSelection(algo->getSelectionString()))
            {
                pipe.addAlgorithm(algo);
            }
        }
         _detIdPipes[detId.rawId()]=pipe;
    }
}

void TrackingRecHitProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
    //build DetId -> PSimHit map
    edm::Handle<std::vector<PSimHit>> simHits;
    event.getByToken(_simHitToken,simHits);
    std::map<unsigned int,std::vector<const PSimHit*>> hitsPerDetId;
    for (unsigned int ihit = 0; ihit < simHits->size(); ++ihit)
    {
        const PSimHit* simHit = &(*simHits)[ihit];
        hitsPerDetId[simHit->detUnitId()].push_back(simHit);
    }
    
    std::auto_ptr<TrackerGSRecHitCollection> recHits(new TrackerGSRecHitCollection());
    std::auto_ptr<TrackerGSMatchedRecHitCollection> matchedRecHits(new TrackerGSMatchedRecHitCollection());
    //run pipes
    for (std::map<unsigned int,std::vector<const PSimHit*>>::iterator simHitsIt = hitsPerDetId.begin(); simHitsIt != hitsPerDetId.end(); ++simHitsIt)
    {
        const DetId& detId = simHitsIt->first;
        std::map<unsigned int, TrackingRecHitPipe>::const_iterator pipeIt = _detIdPipes.find(detId);
        if (pipeIt!=_detIdPipes.cend())
        {
            std::vector<const PSimHit*>& simHits = simHitsIt->second;
            const TrackingRecHitPipe& pipe = pipeIt->second;

            TrackingRecHitProductPtr product = std::make_shared<TrackingRecHitProduct>(detId,simHits);

            product = pipe.produce(product);

            insertRecHits(recHits,product);
            insertMatchedRecHits(matchedRecHits,product);

        }
        else
        {
            //there should be at least an empty pipe for each DetId
            throw cms::Exception("FastSimulation/TrackingRecHitProducer","A PSimHit carries a DetId which does not belong to the TrackerGeometry: "+_trackerTopology->print(simHitsIt->first));
        }
    }

    event.put(recHits,"TrackerGSRecHits");
    event.put(matchedRecHits,"TrackerGSMatchedRecHits");
}

TrackingRecHitProducer::~TrackingRecHitProducer()
{
}


DEFINE_FWK_MODULE(TrackingRecHitProducer);
