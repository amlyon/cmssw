#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "DataFormats/TrajectoryState/interface/TrackCharge.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include <iostream>

using namespace reco;
using namespace std;



GsfTransientTrack::GsfTransientTrack() : 
  GsfTrack(), tkr_(), hasTime(false), timeExt_(0.), dtErrorExt_(0.),
  theField(nullptr), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false), theTIPExtrapolator()
{}

GsfTransientTrack::GsfTransientTrack( const GsfTrack & tk , const MagneticField* field) : 
  GsfTrack(tk),
  tkr_(), hasTime(false), timeExt_(0.), dtErrorExt_(0.),
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false)
{
  
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

GsfTransientTrack::GsfTransientTrack( const GsfTrack & tk , const double time,
                                      const double dtime,
                                      const MagneticField* field) : 
  GsfTrack(tk),
  tkr_(), hasTime(true), timeExt_(time), dtErrorExt_(dtime), 
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false)
{
  
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}


GsfTransientTrack::GsfTransientTrack( const GsfTrackRef & tk , const MagneticField* field) : 
  GsfTrack(*tk), 
  tkr_(tk), hasTime(false), timeExt_(0.), dtErrorExt_(0.), 
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false),
  theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum))
{
  
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

GsfTransientTrack::GsfTransientTrack( const GsfTrackRef & tk ,
                                      const double time,
                                      const double dtime, 
                                      const MagneticField* field) : 
  GsfTrack(*tk), 
  tkr_(tk), hasTime(true), timeExt_(time), dtErrorExt_(dtime), 
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false),
  theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum))
{
  
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

GsfTransientTrack::GsfTransientTrack( const GsfTrack & tk , const MagneticField* field,
				      const edm::ESHandle<GlobalTrackingGeometry>& tg) :
  GsfTrack(tk),
  tkr_(), hasTime(false), timeExt_(0.), dtErrorExt_(0.), 
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false), theTrackingGeometry(tg),
  theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum))
{
  
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

GsfTransientTrack::GsfTransientTrack( const GsfTrack & tk , const double time,
                                      const double dtime,
                                      const MagneticField* field,
				      const edm::ESHandle<GlobalTrackingGeometry>& tg) :
  GsfTrack(tk),
  tkr_(), 
  hasTime(true), timeExt_(time), dtErrorExt_(dtime),
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false), theTrackingGeometry(tg),
  theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum))
{
  
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

GsfTransientTrack::GsfTransientTrack( const GsfTrackRef & tk , const MagneticField* field, 
				      const edm::ESHandle<GlobalTrackingGeometry>& tg) :
  GsfTrack(*tk),
  tkr_(tk), hasTime(false), timeExt_(0.), dtErrorExt_(0.),
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false), theTrackingGeometry(tg),
  theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum))
{
  
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

GsfTransientTrack::GsfTransientTrack( const GsfTrackRef & tk , const double time,
                                      const double dtime,
                                      const MagneticField* field, 
				      const edm::ESHandle<GlobalTrackingGeometry>& tg) :
  GsfTrack(*tk),
  tkr_(tk), hasTime(true), timeExt_(time), dtErrorExt_(dtime),
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false), theTrackingGeometry(tg),
  theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum))
{
  
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

//from gsfMode
GsfTransientTrack::GsfTransientTrack( const GsfTrackRef & tk , const MagneticField* field,
				      const edm::ESHandle<GlobalTrackingGeometry>& tg, 
				      const math::XYZVector& momVal, const int chVal):
  GsfTrack(*tk),
  tkr_(), hasTime(false), timeExt_(0.), dtErrorExt_(0.),
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false), theTrackingGeometry(tg),
  theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum)){

  Basic3DVector<float> pos (tk->referencePoint());
  GlobalPoint gpos(pos);
  Basic3DVector<float> mom (momVal);
  GlobalVector gmom( mom);
  TrackCharge tkCh(chVal);
  GlobalTrajectoryParameters par( gpos, gmom, TrackCharge(tkCh), field);

  //from RecoTracker/TrackProducer/src/GsfTrackProducerBase.cc
  reco::GsfTrack::CovarianceMatrixMode covMode = tk->covarianceMode();
  reco::TrackBase::CovarianceMatrix covMean = tk->covariance();

  AlgebraicSymMatrix55 covErr;
  for (unsigned int iv1 = 0; iv1 < 5; ++iv1) {
    if(iv1 < reco::GsfTrack::dimensionMode) covErr(iv1, iv1) = covMode(iv1, iv1);
    else covErr(iv1, iv1) = covMean(iv1, iv1);
  }
  for (unsigned int iv1 = 0; iv1 < 5; ++iv1) {
    for (unsigned int iv2 = 0; iv2 < iv1; ++iv2) {
      double cov12 = covMean(iv1, iv2) * sqrt(covErr(iv1, iv1) / covMean(iv1, iv1) *
					      covErr(iv2, iv2) / covMean(iv2, iv2));
      covErr(iv1, iv2) = covErr(iv2, iv1) = cov12;
    }
  }

  CurvilinearTrajectoryError err(covErr);
  initialFTS = FreeTrajectoryState( par, err);
}

//from regression
GsfTransientTrack::GsfTransientTrack( const GsfTrackRef & tk , const MagneticField* field,
				      const edm::ESHandle<GlobalTrackingGeometry>& tg,
				      const math::XYZVector& momVal, const int chVal, const float regERatio):
  GsfTrack(*tk),
  tkr_(), hasTime(false), timeExt_(0.), dtErrorExt_(0.),
  theField(field), initialTSOSAvailable(false),
  initialTSCPAvailable(false), blStateAvailable(false), theTrackingGeometry(tg),
  theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum)){

  Basic3DVector<float> pos (tk->referencePoint());
  GlobalPoint gpos(pos);
  Basic3DVector<float> mom (momVal);
  GlobalVector gmom( mom);
  TrackCharge tkCh(chVal);
  GlobalTrajectoryParameters par( gpos, gmom, TrackCharge(tkCh), field);

  //from RecoTracker/TrackProducer/src/GsfTrackProducerBase.cc
  reco::GsfTrack::CovarianceMatrixMode covMode = tk->covarianceMode();
  reco::TrackBase::CovarianceMatrix covMean = tk->covariance();

  AlgebraicSymMatrix55 covErr;
  for (unsigned int iv1 = 0; iv1 < 5; ++iv1) {
    if(iv1 < reco::GsfTrack::dimensionMode) covErr(iv1, iv1) = covMode(iv1, iv1);
    else covErr(iv1, iv1) = covMean(iv1, iv1);
    if(iv1 == 0) covErr(iv1, iv1) = covErr(iv1, iv1) * regERatio * regERatio;
  }
  for (unsigned int iv1 = 0; iv1 < 5; ++iv1) {
    for (unsigned int iv2 = 0; iv2 < iv1; ++iv2) {
      double cov12 = covMean(iv1, iv2) * sqrt(covErr(iv1, iv1) / covMean(iv1, iv1) *
					      covErr(iv2, iv2) / covMean(iv2, iv2));
      covErr(iv1, iv2) = covErr(iv2, iv1) = cov12;
    }
  }

  CurvilinearTrajectoryError err(covErr);
  initialFTS = FreeTrajectoryState( par, err);
}


GsfTransientTrack::GsfTransientTrack( const GsfTransientTrack & tt ) :
  GsfTrack(tt),
  tkr_(tt.persistentTrackRef()), 
  hasTime(tt.hasTime), timeExt_(tt.timeExt_), dtErrorExt_(tt.dtErrorExt_),
  theField(tt.field()), 
  initialFTS(tt.initialFreeState()), initialTSOSAvailable(false),
  initialTSCPAvailable(false),
  theTIPExtrapolator(AnalyticalPropagator(tt.field(), alongMomentum))
{
  if (tt.initialTSOSAvailable) {
    initialTSOS= tt.impactPointState();
    initialTSOSAvailable = true;
  }
  if (tt.initialTSCPAvailable) {
    initialTSCP= tt.impactPointTSCP();
    initialTSCPAvailable = true;
  }
}



void GsfTransientTrack::setES(const edm::EventSetup& setup) {

  setup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry); 

}

void GsfTransientTrack::setTrackingGeometry(const edm::ESHandle<GlobalTrackingGeometry>& tg) {

  theTrackingGeometry = tg;

}

void GsfTransientTrack::setBeamSpot(const BeamSpot& beamSpot)
{
  theBeamSpot = beamSpot;
  blStateAvailable = false;
}


TrajectoryStateOnSurface GsfTransientTrack::impactPointState() const
{
  if (!initialTSOSAvailable) calculateTSOSAtVertex();
  return initialTSOS;
}

TrajectoryStateClosestToPoint GsfTransientTrack::impactPointTSCP() const
{
  if (!initialTSCPAvailable) {
    initialTSCP = builder(initialFTS, initialFTS.position());
    initialTSCPAvailable = true;
  }
  return initialTSCP;
}

TrajectoryStateOnSurface GsfTransientTrack::outermostMeasurementState() const
{
    MultiTrajectoryStateTransform theMTransform;
    return theMTransform.outerStateOnSurface((*this),*theTrackingGeometry,theField);
}

TrajectoryStateOnSurface GsfTransientTrack::innermostMeasurementState() const
{
    MultiTrajectoryStateTransform theMTransform;   
    return theMTransform.innerStateOnSurface((*this),*theTrackingGeometry,theField);
}

void GsfTransientTrack::calculateTSOSAtVertex() const
{
  TransverseImpactPointExtrapolator tipe(theField);
  initialTSOS = tipe.extrapolate(initialFTS, initialFTS.position());
  initialTSOSAvailable = true;
}

TrajectoryStateOnSurface 
GsfTransientTrack::stateOnSurface(const GlobalPoint & point) const
{
  return theTIPExtrapolator.extrapolate(innermostMeasurementState(), point);
}


TrajectoryStateClosestToPoint 
GsfTransientTrack::trajectoryStateClosestToPoint( const GlobalPoint & point ) const
{
  return builder(stateOnSurface(point), point);
}

TrajectoryStateClosestToBeamLine GsfTransientTrack::stateAtBeamLine() const
{
  if (!blStateAvailable) {
    TSCBLBuilderNoMaterial blsBuilder;
    trajectoryStateClosestToBeamLine = blsBuilder(initialFTS, theBeamSpot);
    blStateAvailable = true;
  }
  return trajectoryStateClosestToBeamLine;
}

