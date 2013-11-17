import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.SVfit.tools.objProdConfigurator import *
from TauAnalysis.SVfit.resolutions_cfi import *
from TauAnalysis.SVfit.svFitAlgorithmDiTau_cfi import *
from RecoMET.METProducers.METSigParams_cfi import *

#--------------------------------------------------------------------------------
# produce combinations of electron + muons pairs
#--------------------------------------------------------------------------------

allElecMuPairs = cms.EDProducer("PATElecMuPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatElectronsForElecMuTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    dRmin12 = cms.double(-1.),
    srcMET = cms.InputTag('patMETs'),
    srcPrimaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),                               
    srcGenParticles = cms.InputTag('genParticles'),                                  
    recoMode = cms.string(""),
    doSVreco = cms.bool(True),                          
    svFit = cms.PSet(),                            
    scaleFuncImprovedCollinearApprox = cms.string('1'),                            
    verbosity = cms.untracked.int32(0)
)

#--------------------------------------------------------------------------------
# configure (new) SVfit algorithm
# (using combination of PS + MET likelihoods + logM regularization term
#  to reconstruct mass of tau lepton pair, as described in CMS AN-11-165)
allElecMuPairs.svFit.psKine_MEt_logM_fit = cms.PSet()
allElecMuPairs.svFit.psKine_MEt_logM_fit.config = copy.deepcopy(svFitConfig_template)
allElecMuPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg1 = cms.PSet(
    src = allElecMuPairs.srcLeg1,
    likelihoodFunctions = cms.VPSet(svFitElectronLikelihoodMatrixElement.clone(
        applySinThetaFactor = cms.bool(True))),
    builder = svFitTauToElecBuilder
)
allElecMuPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg2 = cms.PSet(
    src = allElecMuPairs.srcLeg2,
    likelihoodFunctions = cms.VPSet(svFitMuonLikelihoodMatrixElement.clone(
        applySinThetaFactor = cms.bool(True))),
    builder = svFitTauToMuBuilder
)
allElecMuPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.likelihoodFunctions = cms.VPSet(svFitResonanceLikelihoodLogM)
allElecMuPairs.svFit.psKine_MEt_logM_fit.algorithm = cms.PSet(
    pluginName = cms.string("svFitAlgorithmByLikelihoodMaximization"),
    pluginType = cms.string("SVfitAlgorithmByLikelihoodMaximization"),                                    
    minimizer  = cms.vstring("Minuit2", "Migrad"),
    maxObjFunctionCalls = cms.uint32(5000),  
    verbosity = cms.int32(0)
)

allElecMuPairs.svFit.psKine_MEt_int = cms.PSet()
allElecMuPairs.svFit.psKine_MEt_int.config = copy.deepcopy(svFitConfig_template)
allElecMuPairs.svFit.psKine_MEt_int.config.event.resonances.A.daughters.leg1 = cms.PSet(
    src = allElecMuPairs.srcLeg1,
    likelihoodFunctions = cms.VPSet(svFitElectronLikelihoodMatrixElement.clone(
        applySinThetaFactor = cms.bool(False))),
    builder = svFitTauToElecBuilder
)
allElecMuPairs.svFit.psKine_MEt_int.config.event.resonances.A.daughters.leg2 = cms.PSet(
    src = allElecMuPairs.srcLeg2,
    likelihoodFunctions = cms.VPSet(svFitMuonLikelihoodMatrixElement.clone(
        applySinThetaFactor = cms.bool(False))),
    builder = svFitTauToMuBuilder
)
allElecMuPairs.svFit.psKine_MEt_int.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
allElecMuPairs.svFit.psKine_MEt_int.algorithm = cms.PSet(
    pluginName    = cms.string("svFitAlgorithmByIntegration"),
    pluginType    = cms.string("SVfitAlgorithmByIntegration"),
    parameters    = svFitProducerByIntegration.algorithm.parameters,
    vegasOptions  = svFitProducerByIntegration.algorithm.vegasOptions,
    max_or_median = cms.string("max")
)
#--------------------------------------------------------------------------------

produceElecMuPairs = cms.Sequence(allElecMuPairs)

# define additional collections of electron + muon candidates
# with loose track and ECAL isolation applied on electron leg
# (NOTE: to be used for the purpose of factorizing efficiencies
#        of electron isolation from other event selection criteria,
#        in order to avoid problems with limited Monte Carlo statistics)

allElecMuPairsLooseElectronIsolation = cms.EDProducer("PATElecMuPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatElectronsForElecMuTrkIPlooseIsolationCumulative'),
    srcLeg2 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    dRmin12 = cms.double(-1.),
    srcMET = cms.InputTag('patMETs'),
    srcGenParticles = cms.InputTag('genParticles'),                                                        
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                                  
    verbosity = cms.untracked.int32(0)
)

produceElecMuPairsLooseElectronIsolation = cms.Sequence(allElecMuPairsLooseElectronIsolation)

produceElecMuPairsAll = cms.Sequence(produceElecMuPairs * produceElecMuPairsLooseElectronIsolation)


