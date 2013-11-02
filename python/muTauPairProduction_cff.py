import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objProdConfigurator import *
from TauAnalysis.CandidateTools.resolutions_cfi import *
from TauAnalysis.CandidateTools.svFitAlgorithmDiTau_cfi import *
from TauAnalysis.CandidateTools.svFitAlgorithmTauDecayKineMC_cfi import *
from RecoMET.METProducers.METSigParams_cfi import *

#--------------------------------------------------------------------------------
# produce combinations of muon + tau-jet pairs
#--------------------------------------------------------------------------------

allMuTauPairs = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatTausForMuTauElectronVetoCumulative'),
    dRmin12 = cms.double(0.3),
    srcMET = cms.InputTag('patMETs'),
    srcPrimaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),
    srcGenParticles = cms.InputTag('genParticles'),
    recoMode = cms.string(""),
    doSVreco = cms.bool(True),
    svFit = cms.PSet(),
    scaleFuncImprovedCollinearApprox = cms.string('1'),
    doPFMEtSign = cms.bool(True),
    pfMEtSign = cms.PSet(
        srcPFJets = cms.InputTag('ak5PFJets'),
        srcPFCandidates = cms.InputTag('particleFlow'),
        resolution = METSignificance_params,
        dRoverlapPFJet = cms.double(0.3),
        dRoverlapPFCandidate = cms.double(0.1)
    ),
    doMtautauMin = cms.bool(True),                     
    verbosity = cms.untracked.int32(0)
)

#--------------------------------------------------------------------------------
# configure (new) SVfit algorithm
# (using combination of PS + MET likelihoods + logM regularization term
#  to reconstruct mass of tau lepton pair, as described in CMS AN-11-165)
allMuTauPairs.svFit.psKine_MEt_logM_fit = cms.PSet()
allMuTauPairs.svFit.psKine_MEt_logM_fit.config = copy.deepcopy(svFitConfig_template)
allMuTauPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg1 = cms.PSet(
    src = allMuTauPairs.srcLeg1,
    likelihoodFunctions = cms.VPSet(svFitMuonLikelihoodMatrixElement.clone(
        applySinThetaFactor = cms.bool(True))),
    builder = svFitTauToMuBuilder
)
allMuTauPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg2 = cms.PSet(
    src = allMuTauPairs.srcLeg2,
    likelihoodFunctions = cms.VPSet(svFitTauLikelihoodPhaseSpace.clone(
        applySinThetaFactor = cms.bool(True))),
    builder = svFitTauToHadBuilder
)
allMuTauPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.likelihoodFunctions = cms.VPSet(svFitResonanceLikelihoodLogM)
allMuTauPairs.svFit.psKine_MEt_logM_fit.algorithm = cms.PSet(
    pluginName = cms.string("svFitAlgorithmByLikelihoodMaximization"),
    pluginType = cms.string("SVfitAlgorithmByLikelihoodMaximization"),
    minimizer  = cms.vstring("Minuit2", "Migrad"),
    maxObjFunctionCalls = cms.uint32(5000),
    verbosity = cms.int32(0)
)

allMuTauPairs.svFit.psKine_MEt_int = cms.PSet()
allMuTauPairs.svFit.psKine_MEt_int.config = copy.deepcopy(svFitConfig_template)
allMuTauPairs.svFit.psKine_MEt_int.config.event.resonances.A.daughters.leg1 = cms.PSet(
    src = allMuTauPairs.srcLeg1,
    likelihoodFunctions = cms.VPSet(svFitMuonLikelihoodMatrixElement.clone(
        applySinThetaFactor = cms.bool(False))),
    builder = svFitTauToMuBuilder
)
allMuTauPairs.svFit.psKine_MEt_int.config.event.resonances.A.daughters.leg2 = cms.PSet(
    src = allMuTauPairs.srcLeg2,
    likelihoodFunctions = cms.VPSet(svFitTauLikelihoodPhaseSpace.clone(
        applySinThetaFactor = cms.bool(False))),
    builder = svFitTauToHadBuilder
)
allMuTauPairs.svFit.psKine_MEt_int.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
allMuTauPairs.svFit.psKine_MEt_int.algorithm = cms.PSet(
    pluginName    = cms.string("svFitAlgorithmByIntegration"),
    pluginType    = cms.string("SVfitAlgorithmByIntegration"),
    parameters    = svFitProducerByIntegration.algorithm.parameters,
    vegasOptions  = svFitProducerByIntegration.algorithm.vegasOptions,
    max_or_median = cms.string("max")
)
#--------------------------------------------------------------------------------

muTauPairProdConfigurator = objProdConfigurator(
    allMuTauPairs,
    pyModuleName = __name__
)

produceMuTauPairs = muTauPairProdConfigurator.configure(pyNameSpace = locals())

allMuTauPairsPFtype1MET = copy.deepcopy(allMuTauPairs)
allMuTauPairsPFtype1MET.srcMET = cms.InputTag('patPFtype1METs')
produceMuTauPairs += allMuTauPairsPFtype1MET

# define additional collections of muon + tau-jet candidates
# with loose track and ECAL isolation applied on muon leg
# (NOTE: to be used for the purpose of factorizing efficiencies
#        of muon isolation from other event selection criteria,
#        in order to avoid problems with limited Monte Carlo statistics)

allMuTauPairsLooseMuonIsolation = allMuTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatMuonsTrkIPlooseIsolationCumulative'),
)

muTauPairProdConfiguratorLooseMuonIsolation = objProdConfigurator(
    allMuTauPairsLooseMuonIsolation,
    pyModuleName = __name__
)

produceMuTauPairsLooseMuonIsolation = muTauPairProdConfiguratorLooseMuonIsolation.configure(pyNameSpace = locals())

produceMuTauPairsAll = cms.Sequence(produceMuTauPairs * produceMuTauPairsLooseMuonIsolation)
