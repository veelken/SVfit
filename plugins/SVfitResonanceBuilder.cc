#include "TauAnalysis/SVfit/plugins/SVfitResonanceBuilder.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <vector>
#include <string>

SVfitResonanceBuilder::SVfitResonanceBuilder(const edm::ParameterSet& cfg)
  : SVfitResonanceBuilderBase(cfg)
{
  edm::ParameterSet cfg_daughters = cfg.getParameter<edm::ParameterSet>("daughters");
  typedef std::vector<std::string> vstring;
  vstring daughterNames = cfg_daughters.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator daughterName = daughterNames.begin();
	daughterName != daughterNames.end(); ++daughterName ) {
    edm::ParameterSet cfg_daughter = cfg_daughters.getParameter<edm::ParameterSet>(*daughterName);
    edm::ParameterSet cfg_builder = cfg_daughter.getParameter<edm::ParameterSet>("builder");
    cfg_builder.addParameter<std::string>("prodParticleLabel", *daughterName);
    std::string pluginType = cfg_builder.getParameter<std::string>("pluginType");
    SVfitSingleParticleBuilderBase* daughterBuilder =
      SVfitSingleParticleBuilderPluginFactory::get()->create(pluginType, cfg_builder);
    daughterBuilders_.push_back(daughterBuilder);
    ++numDaughterBuilders_;
  }

  if ( cfg.exists("polStates") ) {
    typedef std::vector<std::string> vstring;
    vstring polStates_string = cfg.getParameter<vstring>("polStates");
    for ( vstring::const_iterator polState_string = polStates_string.begin();
	  polState_string != polStates_string.end(); ++polState_string ) {
      int polHandedness = -1;
      if      ( (*polState_string) == "undefined" ) polHandedness = SVfitResonanceHypothesis::kPolUndefined;
      else if ( (*polState_string) == "LR"        ) polHandedness = SVfitResonanceHypothesis::kPolLR;
      else if ( (*polState_string) == "RL"        ) polHandedness = SVfitResonanceHypothesis::kPolRL;
      else if ( (*polState_string) == "LL"        ) polHandedness = SVfitResonanceHypothesis::kPolLL;
      else if ( (*polState_string) == "RR"        ) polHandedness = SVfitResonanceHypothesis::kPolRR;
      else throw cms::Exception("SVfitResonanceBuilder")
	<< " Invalid Configuration Parameter 'polState' = " << (*polState_string) << " !!\n";
      polHandedness_.push_back(polHandedness);
    }
  } else {
    polHandedness_.push_back(SVfitResonanceHypothesis::kPolUndefined);
  }
  numPolStates_ = polHandedness_.size();
}

SVfitResonanceHypothesis* SVfitResonanceBuilder::build(const inputParticleMap& inputParticles) const
{
  SVfitResonanceHypothesis* resonance = SVfitResonanceBuilderBase::build(inputParticles);

//--- set polarization status for resonance
  resonance->polHandedness_ = polHandedness_;
  resonance->numPolStates_ = numPolStates_;

//--- set polarization status for daughters 
  if ( resonance->numDaughters() == 2 ) {
    for ( size_t iDaughter = 0; iDaughter < resonance->numDaughters(); ++iDaughter ) {
      SVfitSingleParticleHypothesis* daughter = resonance->daughter(iDaughter);

      daughter->polHandedness_.resize(numPolStates_);
      daughter->polSign_.resize(numPolStates_);
      daughter->numPolStates_ = numPolStates_;

      for ( unsigned iPolState = 0; iPolState < numPolStates_; ++iPolState ) {
	int resonance_polHandedness = polHandedness_[iPolState];
	int daughter_polHandedness = -1;
	if ( resonance_polHandedness == SVfitResonanceHypothesis::kPolUndefined ) {
	  daughter_polHandedness = SVfitSingleParticleHypothesis::kPolUndefined;
	} else if ( resonance_polHandedness == SVfitResonanceHypothesis::kPolLR ) {
	  if      ( iDaughter == 0 ) daughter_polHandedness = SVfitSingleParticleHypothesis::kPolL;
	  else if ( iDaughter == 1 ) daughter_polHandedness = SVfitSingleParticleHypothesis::kPolR;
	} else if ( resonance_polHandedness == SVfitResonanceHypothesis::kPolRL ) {
	  if      ( iDaughter == 0 ) daughter_polHandedness = SVfitSingleParticleHypothesis::kPolL;
	  else if ( iDaughter == 1 ) daughter_polHandedness = SVfitSingleParticleHypothesis::kPolR;
	} else if ( resonance_polHandedness == SVfitResonanceHypothesis::kPolLL ) {
	  daughter_polHandedness = SVfitSingleParticleHypothesis::kPolL;
	} else if ( resonance_polHandedness == SVfitResonanceHypothesis::kPolRR ) {
	  daughter_polHandedness = SVfitSingleParticleHypothesis::kPolR;
	} 
	assert(daughter_polHandedness != -1);
	daughter->polHandedness_[iPolState] = daughter_polHandedness;

	int daughter_polSign = 0;
	// CV: left-handed  tau- and right-handed tau+ are assigned polarization -1,
	//     right-handed tau- and left-handed  tau+ are assigned polarization +1
	if        ( daughter_polHandedness == SVfitSingleParticleHypothesis::kPolL ) {
	  if      ( daughter->particle()->charge() < -0.5 ) daughter_polSign = -1;
	  else if ( daughter->particle()->charge() > +0.5 ) daughter_polSign = +1;
	} else if ( daughter_polHandedness == SVfitSingleParticleHypothesis::kPolR ) {
	  if      ( daughter->particle()->charge() < -0.5 ) daughter_polSign = +1;
	  else if ( daughter->particle()->charge() > +0.5 ) daughter_polSign = -1;
	} 
	daughter->polSign_[iPolState] = daughter_polSign;
      }
    }
  } else if ( numPolStates_ == 1 && polHandedness_[0] == SVfitResonanceHypothesis::kPolUndefined ) {
    for ( size_t iDaughter = 0; iDaughter < resonance->numDaughters(); ++iDaughter ) {
      SVfitSingleParticleHypothesis* daughter = resonance->daughter(iDaughter);
      daughter->polHandedness_.resize(numPolStates_);
      daughter->polHandedness_[0] = SVfitSingleParticleHypothesis::kPolUndefined;
      daughter->polSign_.resize(numPolStates_);
      daughter->polSign_[0] = 0;
      daughter->numPolStates_ = numPolStates_;
    }
  } else throw cms::Exception("SVfitResonanceBuilder")
      << " Support for Polarization not implemented for case of " << resonance->numDaughters() << " daughters yet !!\n";

  return resonance;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitResonanceBuilderPluginFactory, SVfitResonanceBuilder, "SVfitResonanceBuilder");
