#include "TauAnalysis/SVfit/interface/SVfitEventBuilderBase.h"

#include "DataFormats/Common/interface/Handle.h"

#include "TauAnalysis/SVfit/interface/SVfitParameter.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitTauDecayHypothesis.h"

using namespace svFit_namespace;

SVfitEventBuilderBase::SVfitEventBuilderBase(const edm::ParameterSet& cfg)
  : SVfitBuilderBase(cfg),
    numResonanceBuilders_(0),
    eventVertexRefitAlgorithm_(0),
    idxFitParameter_pvShiftX_(-1),
    idxFitParameter_pvShiftY_(-1),
    idxFitParameter_pvShiftZ_(-1)    
{
  srcBeamSpot_ = cfg.getParameter<edm::InputTag>("srcBeamSpot");

  edm::ParameterSet cfg_resonances = cfg.getParameter<edm::ParameterSet>("resonances");
  typedef std::vector<std::string> vstring;
  vstring resonanceNames = cfg_resonances.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator resonanceName = resonanceNames.begin();
	resonanceName != resonanceNames.end(); ++resonanceName ) {
    edm::ParameterSet cfg_resonance = cfg_resonances.getParameter<edm::ParameterSet>(*resonanceName);
    edm::ParameterSet cfg_builder = cfg_resonance.getParameter<edm::ParameterSet>("builder");
    cfg_builder.addParameter<edm::ParameterSet>("daughters", cfg_resonance.getParameter<edm::ParameterSet>("daughters"));
    cfg_builder.addParameter<std::string>("prodResonanceLabel", *resonanceName);
    std::string pluginType = cfg_builder.getParameter<std::string>("pluginType");
    SVfitResonanceBuilderBase* resonanceBuilder =
      SVfitResonanceBuilderPluginFactory::get()->create(pluginType, cfg_builder);
    resonanceBuilders_.push_back(resonanceBuilder);
    ++numResonanceBuilders_;
  }

  eventVertexRefitAlgorithm_ = new SVfitEventVertexRefitter(cfg);
}

SVfitEventBuilderBase::~SVfitEventBuilderBase()
{
  for ( std::vector<SVfitResonanceBuilderBase*>::iterator it = resonanceBuilders_.begin();
	it != resonanceBuilders_.end(); ++it ) {
    delete (*it);
  }

  delete eventVertexRefitAlgorithm_;
}

void SVfitEventBuilderBase::beginJob(SVfitAlgorithmBase* algorithm)
{
  for ( std::vector<SVfitResonanceBuilderBase*>::iterator resonanceBuilder = resonanceBuilders_.begin();
	resonanceBuilder != resonanceBuilders_.end(); ++resonanceBuilder ) {
    (*resonanceBuilder)->beginJob(algorithm);
  }

  idxFitParameter_pvShiftX_ = getFitParameterIdx(algorithm, "*", svFit_namespace::kPV_shiftX, true); // optional parameter
  idxFitParameter_pvShiftY_ = getFitParameterIdx(algorithm, "*", svFit_namespace::kPV_shiftY, true);
  idxFitParameter_pvShiftZ_ = getFitParameterIdx(algorithm, "*", svFit_namespace::kPV_shiftZ, true);
  doFitParameter_pvShift_ = (idxFitParameter_pvShiftX_ != -1 && idxFitParameter_pvShiftY_ != -1 && idxFitParameter_pvShiftZ_ != -1);

  // CV: always perform event vertex refit
  doEventVertexRefit_ = true;
}

void SVfitEventBuilderBase::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  SVfitBuilderBase::beginEvent(evt, es);
  for ( std::vector<SVfitResonanceBuilderBase*>::iterator resonanceBuilder = resonanceBuilders_.begin();
	resonanceBuilder != resonanceBuilders_.end(); ++resonanceBuilder ) {
    (*resonanceBuilder)->beginEvent(evt, es);
  }

  if ( doEventVertexRefit_ ) {
    eventVertexRefitAlgorithm_->beginEvent(evt, es);
  }

  edm::Handle<reco::BeamSpot> beamSpot;
  evt.getByLabel(srcBeamSpot_, beamSpot);

  beamSpot_ = (*beamSpot);
}

SVfitEventHypothesis* SVfitEventBuilderBase::build(const inputParticleMap& inputParticles, const reco::Vertex* eventVertex) const
{
  inputParticleMap::const_iterator metPtr = inputParticles.find("met");
  assert(metPtr != inputParticles.end());

  SVfitEventHypothesis* event = new SVfitEventHypothesis(metPtr->second);

  reco::Candidate::LorentzVector p4(0,0,0,0);

  for ( std::vector<SVfitResonanceBuilderBase*>::const_iterator resonanceBuilder = resonanceBuilders_.begin();
	resonanceBuilder != resonanceBuilders_.end(); ++resonanceBuilder ) {
    SVfitResonanceHypothesis* resonance = (*resonanceBuilder)->build(inputParticles);    
    resonance->setEventHypothesis(event);

    p4 += resonance->p4();

    lastBuiltResonances_[*resonanceBuilder] = event->numResonances();

    event->resonances_.push_back(resonance);
  }

  event->p4_ = p4;

  event->eventVertexIsValid_ = false;

//--- refit primary event vertex, excluding tracks of tau decay products
  std::vector<const reco::Track*> tauTracks;
  size_t numResonances = event->numResonances();
  for ( size_t iResonance = 0; iResonance < numResonances; ++iResonance ) {
    SVfitResonanceHypothesis* resonance = event->resonance(iResonance);
    size_t numDaughters = resonance->numDaughters();
    for ( size_t iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
      const SVfitSingleParticleHypothesis* daughter = resonance->daughter(iDaughter);
      const SVfitTauDecayHypothesis* tauDecayProducts = dynamic_cast<const SVfitTauDecayHypothesis*>(daughter);
      if ( tauDecayProducts ) {
	tauTracks.insert(tauTracks.begin(), tauDecayProducts->tracks().begin(), tauDecayProducts->tracks().end());
      }
    }
  }
  
  if ( doEventVertexRefit_ ) {
    TransientVertex eventVertex_refitted = eventVertexRefitAlgorithm_->refit(eventVertex, &tauTracks);
    if ( eventVertex_refitted.isValid() ) {      
      event->eventVertexPos_(0) = eventVertex_refitted.position().x();
      event->eventVertexPos_(1) = eventVertex_refitted.position().y();
      event->eventVertexPos_(2) = eventVertex_refitted.position().z();

      event->eventVertexCov_ = eventVertex_refitted.positionError().matrix_new();

      event->eventVertexIsValid_ = true;
    } else {
      // take estimate of event vertex x/y position from beam-spot;
      // take z position from original vertex and assign large uncertainty
      double zVertex, sigmaZ;
      if ( eventVertex && eventVertex->isValid() ) {
	zVertex = eventVertex->position().z();
	sigmaZ = square(2.); // account for bias of tauTrack to reconstructed vertex position by assigning sigmaZ = 2 cm
      } else {
	zVertex = 0.;
	sigmaZ = square(10.); // no information about z position of vertex available, 
	                      // assign width of zVertex distribution as uncertainty
      }
      
      event->eventVertexPos_(0) = beamSpot_.x(zVertex);
      event->eventVertexPos_(1) = beamSpot_.y(zVertex);
      event->eventVertexPos_(2) = zVertex;
      
      for ( int iRow = 0; iRow < 3; ++iRow ) {
	for ( int iColumn = 0; iColumn < 3; ++iColumn ) {
	  if      ( iRow <= 1 && iColumn <= 1 ) event->eventVertexCov_(iRow, iColumn) = beamSpot_.covariance(iRow, iColumn);
	  else if ( iRow == 2 && iColumn == 2 ) event->eventVertexCov_(iRow, iColumn) = sigmaZ;
	  else                                  event->eventVertexCov_(iRow, iColumn) = 0.;
	}
      }

      if ( eventVertex && eventVertex->isValid() ) event->eventVertexIsValid_ = true;
      else event->eventVertexIsValid_ = false;
    }
  } else {
    if ( eventVertex && eventVertex->isValid() ) {
      event->eventVertexPos_(0) = eventVertex->position().x();
      event->eventVertexPos_(1) = eventVertex->position().y();
      event->eventVertexPos_(2) = eventVertex->position().z();

      event->eventVertexCov_ = eventVertex->covariance();

      event->eventVertexIsValid_ = true;
    } else {
      event->eventVertexPos_(0) = beamSpot_.x(0.);
      event->eventVertexPos_(1) = beamSpot_.y(0.);
      event->eventVertexPos_(2) = 0.;

      for ( int iRow = 0; iRow < 3; ++iRow ) {
	for ( int iColumn = 0; iColumn < 3; ++iColumn ) {
	  if      ( iRow <= 1 && iColumn <= 1 ) event->eventVertexCov_(iRow, iColumn) = beamSpot_.covariance(iRow, iColumn);
	  else if ( iRow == 2 && iColumn == 2 ) event->eventVertexCov_(iRow, iColumn) = square(10.); // no information about z position of vertex available, 
                                                                                                     // assign width of zVertex distribution as uncertainty
	  else                                  event->eventVertexCov_(iRow, iColumn) = 0.;
	}
      }

      event->eventVertexIsValid_ = false;
    }
  }
 
  if ( verbosity_ ) {
    std::cout << "<SVfitEventBuilderBase::build>:" << std::endl;
    std::cout << " pluginName = " << pluginName_ << std::endl;
    std::cout << " recEventVertexIsValid = " << event->eventVertexIsValid_ << std::endl;
    printVector(" recEventVertexPos", event->eventVertexPos());
    printMatrix(" recEventVertexCov", event->eventVertexCov_);
  }

  event->barcode_ = barcodeCounter_;
  ++barcodeCounter_;

  return event;
}

bool SVfitEventBuilderBase::applyFitParameter(SVfitEventHypothesis* event, const double* param) const
{
  bool isValidSolution = true;

  if ( doFitParameter_pvShift_ ) {
    event->eventVertexShift_(0) = param[idxFitParameter_pvShiftX_];
    event->eventVertexShift_(1) = param[idxFitParameter_pvShiftY_];
    event->eventVertexShift_(2) = param[idxFitParameter_pvShiftZ_];
  }

  reco::Candidate::LorentzVector dp4(0,0,0,0);

  for ( unsigned iResonanceBuilder = 0; iResonanceBuilder < numResonanceBuilders_; ++iResonanceBuilder ) {
    SVfitResonanceHypothesis* resonance = event->resonance(iResonanceBuilder);
    isValidSolution &= resonanceBuilders_[iResonanceBuilder]->applyFitParameter(resonance, param);

    dp4 += resonance->dp4_fitted();
  }

  event->dp4_ = dp4;

  return isValidSolution;
}

void SVfitEventBuilderBase::print(std::ostream& stream) const
{
  stream << "<SVfitEventBuilderBase::print>:" << std::endl;
  stream << " pluginName = " << pluginName_ << std::endl;
  stream << " pluginType = " << pluginType_ << std::endl;
  for ( std::vector<SVfitResonanceBuilderBase*>::const_iterator resonanceBuilder = resonanceBuilders_.begin();
	resonanceBuilder != resonanceBuilders_.end(); ++resonanceBuilder ) {
    (*resonanceBuilder)->print(stream);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(SVfitEventBuilderPluginFactory, "SVfitEventBuilderPluginFactory");
