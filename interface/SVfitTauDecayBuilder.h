#ifndef TauAnalysis_SVfit_SVfitTauDecayBuilder_h
#define TauAnalysis_SVfit_SVfitTauDecayBuilder_h

/** \class SVfitTauDecayBuilder
 *
 * Base-class for building objects that come from tau decays.
 *
 * \author Evan K. Friis, Christian Veelken, UC Davis
 *
 */

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/View.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TauAnalysis/SVfit/interface/SVfitSingleParticleBuilderBase.h"
#include "TauAnalysis/SVfit/interface/SVfitTrackService.h"
#include "TauAnalysis/SVfit/interface/SVfitDecayVertexFitter.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitTauDecayHypothesis.h"

class SVfitSingleParticleHypothesis;
class SVfitAlgorithmBase;

class SVfitTauDecayBuilder : public SVfitSingleParticleBuilderBase
{
  public:
    SVfitTauDecayBuilder(const edm::ParameterSet&);
    virtual ~SVfitTauDecayBuilder();

    // Setup the parameters of the fit.
    virtual void beginJob(SVfitAlgorithmBase*);

    // Access TrackBuilder
    virtual void beginEvent(const edm::Event&, const edm::EventSetup&);

    // Add Tracking information
    // 
    // NOTE: this function needs to be called after the event has been build,
    //       as the tracking information depends on the primary event vertex
    //      (which in turn depends on the tracks of daughter particles)
    //
    virtual void finalize(SVfitSingleParticleHypothesis*) const;

    // Build the tau decay hypothesis from the fit parameters
    virtual bool applyFitParameter(SVfitSingleParticleHypothesis*, const double*) const;

    /* Abstract functions overridden by the different decay type builders */
    // Overridden to allocate the specific decay type.
    virtual bool nuSystemIsMassless() const = 0;
    // The decay mode
    virtual int getDecayMode(const reco::Candidate*) const = 0;
    // Get the track(s) associated to a given Candidate
    virtual std::vector<const reco::Track*> extractTracks(const reco::Candidate*) const = 0;

    virtual void print(std::ostream&) const;

  protected:
    // Initialize data-members common to tau --> e/mu and tau --> had decays
    void initialize(SVfitTauDecayHypothesis*, const reco::Candidate*) const;

    SVfitAlgorithmBase* algorithm_;

    edm::Service<SVfitTrackService> trackService_;
    const TransientTrackBuilder* trackBuilder_;

    mutable std::vector<reco::TransientTrack> selectedTracks_;

    SVfitDecayVertexFitter* decayVertexFitAlgorithm_;
    bool fitDecayVertex_;

    unsigned trackMinNumHits_;
    unsigned trackMinNumPixelHits_;
    double   trackMaxChi2DoF_;
    double   trackMaxDeltaPoverP_;
    double   trackMinPt_;

    int idxFitParameter_visEnFracX_;
    int idxFitParameter_phi_lab_;
    int idxFitParameter_visMass_;   // used for hadronic decays only
    int idxFitParameter_nuInvMass_; // used for leptonic decays only
    int idxFitParameter_decayDistance_shift_;
    
    /// optional parameters for setting reconstructed to Monte Carlo truth values
    edm::InputTag srcGenTaus_;
    typedef edm::View<reco::GenParticle> GenParticleView;
    edm::Handle<GenParticleView> genParticles_;
    double dRmatch_;

    bool fixToGenVisEnFracX_;
    bool fixToGenPhiLab_;
    bool fixToGenVisMass_;
    bool fixToGenNuInvMass_;
    bool fixRecToGenVertex_;
    bool fixToGenDecayDistance_;
    bool fixToGenVisP4_;
    bool initializeToGen_;

    mutable double genVisEnFracX_;
    mutable double genPhiLab_;
    mutable double genVisMass_;
    mutable double genNuInvMass_;
    mutable reco::Candidate::Point genDecayVertex_;
    mutable double genDecayDistance_;
    mutable reco::Candidate::LorentzVector genVisP4_;
};

void applyOptionalFitParameter(const double*, int, double&);

#endif /* end of include guard: TauAnalysis_SVfit_SVfitTauDecayBuilder_h */
