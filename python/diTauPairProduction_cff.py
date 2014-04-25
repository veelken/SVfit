import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.SVfit.tools.objProdConfigurator import *
from TauAnalysis.SVfit.resolutions_cfi import *
from TauAnalysis.SVfit.svFitAlgorithmDiTau_cfi import *
from RecoMET.METProducers.METSigParams_cfi import *

#--------------------------------------------------------------------------------
# produce combinations of tau-jet + tau-jet pairs
#--------------------------------------------------------------------------------

allDiTauPairs = cms.EDProducer("PATDiTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('cleanPatTaus'),
    srcLeg2 = cms.InputTag('cleanPatTaus'),
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

produceDiTauPairs = cms.Sequence(allDiTauPairs)

produceDiTauPairsAll = cms.Sequence(produceDiTauPairs)

