import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.SVfit.tools.objProdConfigurator import *
from TauAnalysis.SVfit.resolutions_cfi import *
from TauAnalysis.SVfit.svFitAlgorithmDiTau_cfi import *
from RecoMET.METProducers.METSigParams_cfi import *

#--------------------------------------------------------------------------------
# produce combinations of electron + tau-jet pairs
#--------------------------------------------------------------------------------

allElecTauPairs = cms.EDProducer("PATElecTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatElectronsForElecTauTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatTausForElecTauMuonVetoCumulative'),
    dRmin12 = cms.double(0.3),
    srcMET = cms.InputTag('patMETs'),
    srcPrimaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),                             
    srcGenParticles = cms.InputTag('genParticles'),                                   
    recoMode = cms.string(""),
    doSVreco = cms.bool(True),                          
    svFit = cms.PSet(),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                            
    doPFMEtSign = cms.bool(False),
    pfMEtSign = cms.PSet(
        srcPFJets = cms.InputTag('ak5PFJets'),
        srcPFCandidates = cms.InputTag('particleFlow'),
        resolution = METSignificance_params,
        dRoverlapPFJet = cms.double(0.3),
        dRoverlapPFCandidate = cms.double(0.1)
    ),
    doMtautauMin = cms.bool(False),                             
    verbosity = cms.untracked.int32(0)
)

#--------------------------------------------------------------------------------
# configure (new) SVfit algorithm
# (using combination of PS + MET likelihoods + logM regularization term
#  to reconstruct mass of tau lepton pair, as described in CMS AN-11-165)
allElecTauPairs.svFit.psKine_MEt_logM_fit = cms.PSet()
allElecTauPairs.svFit.psKine_MEt_logM_fit.config = copy.deepcopy(svFitConfig_template)
allElecTauPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg1 = cms.PSet(
    src = allElecTauPairs.srcLeg1,
    likelihoodFunctions = cms.VPSet(svFitElectronLikelihoodMatrixElement.clone(
        applySinThetaFactor = cms.bool(True))),
    builder = svFitTauToElecBuilder
)
allElecTauPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg2 = cms.PSet(
    src = allElecTauPairs.srcLeg2,
    likelihoodFunctions = cms.VPSet(svFitTauLikelihoodPhaseSpace.clone(
        applySinThetaFactor = cms.bool(True))),
    builder = svFitTauToHadBuilder
)
allElecTauPairs.svFit.psKine_MEt_logM_fit.config.event.resonances.A.likelihoodFunctions = cms.VPSet(svFitResonanceLikelihoodLogM)
allElecTauPairs.svFit.psKine_MEt_logM_fit.algorithm = cms.PSet(
    pluginName = cms.string("svFitAlgorithmByLikelihoodMaximization"),
    pluginType = cms.string("SVfitAlgorithmByLikelihoodMaximization"),                                    
    minimizer  = cms.vstring("Minuit2", "Migrad"),
    maxObjFunctionCalls = cms.uint32(5000),  
    verbosity = cms.int32(0)
)

allElecTauPairs.svFit.psKine_MEt_int = cms.PSet()
allElecTauPairs.svFit.psKine_MEt_int.config = copy.deepcopy(svFitConfig_template)
allElecTauPairs.svFit.psKine_MEt_int.config.event.resonances.A.daughters.leg1 = cms.PSet(
    src = allElecTauPairs.srcLeg1,
    likelihoodFunctions = cms.VPSet(svFitElectronLikelihoodMatrixElement.clone(
        applySinThetaFactor = cms.bool(False))),
    builder = svFitTauToElecBuilder
)
allElecTauPairs.svFit.psKine_MEt_int.config.event.resonances.A.daughters.leg2 = cms.PSet(
    src = allElecTauPairs.srcLeg2,
    likelihoodFunctions = cms.VPSet(svFitTauLikelihoodPhaseSpace.clone(
        applySinThetaFactor = cms.bool(False))),
    builder = svFitTauToHadBuilder
)
allElecTauPairs.svFit.psKine_MEt_int.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
allElecTauPairs.svFit.psKine_MEt_int.algorithm = cms.PSet(
    pluginName    = cms.string("svFitAlgorithmByIntegration"),
    pluginType    = cms.string("SVfitAlgorithmByIntegration"),
    parameters    = svFitProducerByIntegration.algorithm.parameters,
    vegasOptions  = svFitProducerByIntegration.algorithm.vegasOptions,
    max_or_median = cms.string("max")
)
#--------------------------------------------------------------------------------

elecTauPairProdConfigurator = objProdConfigurator(
	allElecTauPairs,
	pyModuleName = __name__
)
                
produceElecTauPairs = elecTauPairProdConfigurator.configure(pyNameSpace = locals())

# define additional collections of electron + tau-jet candidates
# with loose track and ECAL isolation (and conversion rejection) applied on electron leg
# (NOTE: to be used for the purpose of factorizing efficiencies
#        of electron isolation from other event selection criteria,
#        in order to avoid problems with limited Monte Carlo statistics)

allElecTauPairsLooseElectronIsolation = allElecTauPairs.clone()
allElecTauPairsLooseElectronIsolation.srcLeg1 = cms.InputTag('selectedPatElectronsForElecTauTrkIPlooseIsolationCumulative')

elecTauPairProdConfiguratorLooseElectronIsolation = objProdConfigurator(
	allElecTauPairsLooseElectronIsolation,
	pyModuleName = __name__
)               
                
produceElecTauPairsLooseElectronIsolation = elecTauPairProdConfiguratorLooseElectronIsolation.configure(pyNameSpace = locals())

produceElecTauPairsAll = cms.Sequence( produceElecTauPairs )
        #* produceElecTauPairsLooseElectronIsolation )
