import FWCore.ParameterSet.Config as cms

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
import RecoMET.METProducers.METSigParams_cfi as met_config
from TauAnalysis.SVfit.svFitAlgorithmTrackQualityCuts_cfi import trackQualityCuts

svFitTrackService = cms.Service("SVfitTrackService")

svFitElectronLikelihoodPhaseSpace = cms.PSet(
    pluginName = cms.string("svFitTauToElecLikelihoodPhaseSpace"),
    pluginType = cms.string("SVfitTauToElecLikelihoodPhaseSpace"),
    applySinThetaFactor = cms.bool(False),
    verbosity = cms.int32(0)
)

svFitElectronLikelihoodMatrixElement = cms.PSet(
    pluginName = cms.string("svFitTauToElecLikelihoodMatrixElement"),
    pluginType = cms.string("SVfitTauToElecLikelihoodMatrixElement"),
    applySinThetaFactor = cms.bool(False),
    verbosity = cms.int32(0)
)

svFitElectronLikelihoodDummy = cms.PSet(
    pluginName = cms.string("svFitTauToElecLikelihoodDummy"),
    pluginType = cms.string("SVfitTauToElecLikelihoodDummy"),
    verbosity = cms.int32(0)
)

svFitTauToElecBuilder = cms.PSet(
    pluginName = cms.string("svFitTauToElecBuilder"),
    pluginType = cms.string("SVfitTauToElecBuilder"),
    trackQualityCuts = trackQualityCuts,
    fitDecayVertex = cms.bool(False),
    verbosity = cms.int32(0)
)

svFitMuonLikelihoodPhaseSpace = cms.PSet(
    pluginName = cms.string("svFitTauToMuLikelihoodPhaseSpace"),
    pluginType = cms.string("SVfitTauToMuLikelihoodPhaseSpace"),
    applySinThetaFactor = cms.bool(False),
    verbosity = cms.int32(0)
)

svFitMuonLikelihoodMatrixElement = cms.PSet(
    pluginName = cms.string("svFitTauToMuLikelihoodMatrixElement"),
    pluginType = cms.string("SVfitTauToMuLikelihoodMatrixElement"),
    applySinThetaFactor = cms.bool(False),
    verbosity = cms.int32(0)
)

svFitMuonLikelihoodDummy = cms.PSet(
    pluginName = cms.string("svFitTauToMuLikelihoodDummy"),
    pluginType = cms.string("SVfitTauToMuLikelihoodDummy"),
    verbosity = cms.int32(0)
)

svFitTauToMuBuilder = cms.PSet(
    pluginName = cms.string("svFitTauToMuBuilder"),
    pluginType = cms.string("SVfitTauToMuBuilder"),
    trackQualityCuts = trackQualityCuts,
    fitDecayVertex = cms.bool(False),
    verbosity = cms.int32(0)
)

svFitTauLikelihoodPhaseSpace = cms.PSet(
    pluginName = cms.string("svFitTauToHadLikelihoodPhaseSpace"),
    pluginType = cms.string("SVfitTauToHadLikelihoodPhaseSpace"),
    applySinThetaFactor = cms.bool(False),
    varyVisMass = cms.bool(False),
    inputFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitGenTauHadMassPDF.root"),
    histogramName = cms.string('DQMData/genTauMassAnalyzer/genTauJetMass'),
    shiftVisMass = cms.bool(False),
    visMassResolution = cms.VPSet(
        cms.PSet(
            tauDecayModes = cms.vint32(0),
            inputFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitVisMassAndPtResolutionPDF.root"),
            histogramName = cms.string('recMinusGenTauMass_recDecayModeEq0'),
        ),
        cms.PSet(
            tauDecayModes = cms.vint32(1, 2),
            inputFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitVisMassAndPtResolutionPDF.root"),
            histogramName = cms.string('recMinusGenTauMass_recDecayModeEq1'),
        ),
        cms.PSet(
            tauDecayModes = cms.vint32(10),
            inputFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitVisMassAndPtResolutionPDF.root"),
            histogramName = cms.string('recMinusGenTauMass_recDecayModeEq10'),
        )
    ),
    shiftVisPt = cms.bool(False),
    visPtResolution = cms.VPSet(
        cms.PSet(
            tauDecayModes = cms.vint32(0),
            inputFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitVisMassAndPtResolutionPDF.root"),
            histogramName = cms.string('recTauPtDivGenTauPt_recDecayModeEq0'),
        ),
        cms.PSet(
            tauDecayModes = cms.vint32(1, 2),
            inputFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitVisMassAndPtResolutionPDF.root"),
            histogramName = cms.string('recTauPtDivGenTauPt_recDecayModeEq1'),
        ),
        cms.PSet(
            tauDecayModes = cms.vint32(10),
            inputFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitVisMassAndPtResolutionPDF.root"),
            histogramName = cms.string('recTauPtDivGenTauPt_recDecayModeEq10'),
        )
    ),
    verbosity = cms.int32(0)
)

svFitTauLikelihoodMatrixElement = cms.PSet(
    pluginName = cms.string("svFitTauToHadLikelihoodMatrixElement"),
    pluginType = cms.string("SVfitTauToHadLikelihoodMatrixElement"),
    VMshapeFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitVMpdfs.root"),
    recToGenTauDecayModeMapFileName = cms.FileInPath("TauAnalysis/SVfit/data/svFitRecToGenTauDecayModeMap.root"),
    applySinThetaFactor = cms.bool(False),
    verbosity = cms.int32(0)
)

svFitTauLikelihoodDummy = cms.PSet(
    pluginName = cms.string("svFitTauToHadLikelihoodDummy"),
    pluginType = cms.string("SVfitTauToHadLikelihoodDummy"),
    verbosity = cms.int32(0)
)

svFitTauToHadBuilder = cms.PSet(
    pluginName = cms.string("svFitTauToHadBuilder"),
    pluginType = cms.string("SVfitTauToHadBuilder"),
    trackQualityCuts = trackQualityCuts,
    fitDecayVertex = cms.bool(True),
    verbosity = cms.int32(0)
)

svFitResonanceLikelihoodPhaseSpace = cms.PSet(
    pluginName = cms.string("svFitResonanceLikelihoodPhaseSpace"),
    pluginType = cms.string("SVfitResonanceLikelihoodPhaseSpace"),
    power = cms.double(1.0),
    verbosity = cms.int32(0)
)

svFitResonanceLikelihoodLogM = cms.PSet(
    pluginName = cms.string("svFitResonanceLikelihoodLogM"),
    pluginType = cms.string("SVfitResonanceLikelihoodRegularization"),
    nll = cms.string("TMath::Log(mass)"),
    power = cms.double(1.0)
)

svFitResonanceLikelihoodLogEff = cms.PSet(
    pluginName = cms.string("svFitResonanceLikelihoodEff_power100"),
    pluginType = cms.string("SVfitResonanceLikelihoodRegularization"),
    nll = cms.string("TMath::Log(TMath::Max(5.00e-3, 4.21e-2*(2.52e-2 + TMath::Erf((x - 4.40e+1)*6.90e-3))))"),
    power = cms.double(1.0)
)

svFitResonanceLikelihoodPolarization = cms.PSet(
    pluginName = cms.string("svFitResonanceLikelihoodPolarization"),
    pluginType = cms.string("SVfitResonanceLikelihoodPolarization"),
    LR = cms.PSet(
        formula = cms.string("[0]*([1] + [2] - x)/[2]"), # CV: linearly decrease weight of Z specific polarization
                                                         #     between mZ and mZ + 30 GeV
        xMin = cms.double(91.188),
        xMax = cms.double(121.188),
        parameter = cms.PSet(
            par0 = cms.double(0.576),
            par1 = cms.double(91.188), # mZ / GeV
            par2 = cms.double(30.)     # width of transition region between Z and Higgs specific polarization
        )
    ),
    RL = cms.PSet(
        formula = cms.string("[0]*([1] + [2] - x)/[2]"), # CV: linearly decrease weight of Z specific polarization
                                                         #     between mZ and mZ + 30 GeV
        xMin = cms.double(91.188),
        xMax = cms.double(121.188),
        parameter = cms.PSet(
            par0 = cms.double(0.424),
            par1 = cms.double(91.188), # mZ / GeV
            par2 = cms.double(30.)     # width of transition region between Z and Higgs specific polarization
        )
    ),
    LL = cms.PSet(
        formula = cms.string("[0]*(x - [1])/[2]"), # CV: linearly increase weight of Higgs specific polarization
                                                   #     between mZ and mZ + 30 GeV
        xMin = cms.double(91.188), 
        xMax = cms.double(121.188),
        parameter = cms.PSet(
            par0 = cms.double(0.50),
            par1 = cms.double(91.188), # mZ / GeV
            par2 = cms.double(30.),    # width of transition region between Z and Higgs specific polarization
        )
    ),
    RR = cms.PSet(
        formula = cms.string("[0]*(x - [1])/[2]"), # CV: linearly increase weight of Higgs specific polarization
                                                   #     between mZ and mZ + 30 GeV
        xMin = cms.double(91.188), 
        xMax = cms.double(121.188),
        parameter = cms.PSet(
            par0 = cms.double(0.50),
            par1 = cms.double(91.188), # mZ / GeV
            par2 = cms.double(30.),    # width of transition region between Z and Higgs specific polarization
        )
    ),
    power = cms.double(1.0),
    verbosity = cms.int32(0)
)

svFitResonanceLikelihoodPrior = cms.PSet(
    pluginName = cms.string("svFitResonanceLikelihoodPrior"),
    pluginType = cms.string("SVfitResonanceLikelihoodPrior"),
    formula = cms.string("1. + [0]*TMath::Gaus(x, [1], [2])"),
    xMin = cms.double(91.188), # do not apply prior probability correction below mZ
    xMax = cms.double(1.e+3),
    parameter = cms.PSet(
        par0 = cms.double(3.),
        par1 = cms.double(91.188), # mZ / GeV
        par2 = cms.double(2.495)   # GammaZ / GeV
    ),
    power = cms.double(1.0)
)

svFitResonanceBuilder = cms.PSet(
    pluginName = cms.string("svFitResonanceBuilder"),
    pluginType = cms.string("SVfitResonanceBuilder"),
    polStates = cms.vstring( # polarization states to be considered when evaluating likelihoods
        ##"LR", "RL", # Z case
        ##"LL", "RR"  # Higgs case
        "undefined"
    )
)

tailProbCorr_Data_2011 = cms.PSet( # use for 2011 Data
    formula = cms.string("[0] + [1]*x + [2]*0.5*(3.*x*x - 1.) + [3]*0.2*(5.*x*x*x -3.*x) + [4]*0.125*(35.*x*x*x*x - 30.*x*x + 3.)"),
    xMin = cms.double(0.),
    xMax = cms.double(5.),
    parameter = cms.PSet(
        par0 = cms.double(9.94505e-01),
        par1 = cms.double(5.85923e-04),
        par2 = cms.double(1.41680e-03),
        par3 = cms.double(-6.20178e-05),
        par4 = cms.double(5.55139e-05)
    )
)

tailProbCorr_MC_2011 = cms.PSet( # use for 2011 Monte Carlo
    formula = cms.string("[0] + [1]*x + [2]*0.5*(3.*x*x - 1.) + [3]*0.2*(5.*x*x*x -3.*x) + [4]*0.125*(35.*x*x*x*x - 30.*x*x + 3.)"),
    xMin = cms.double(0.),
    xMax = cms.double(5.),
    parameter = cms.PSet(
        par0 = cms.double(9.94448e-01),
        par1 = cms.double(-6.74415e-04),
        par2 = cms.double(1.52430e-03),
        par3 = cms.double(1.07572e-04),
        par4 = cms.double(5.27405e-05)
    )
)

svFitEventLikelihoodMEt2 = cms.PSet(
    pluginName = cms.string("svFitEventLikelihoodMEt2"),
    pluginType = cms.string("SVfitEventLikelihoodMEt2"),
    srcPFJets = cms.InputTag('ak5PFJets'),
    srcPFCandidates = cms.InputTag('particleFlow'),
    resolution = met_config.METSignificance_params,
    dRoverlapPFJet = cms.double(0.3),
    dRoverlapPFCandidate = cms.double(0.1),
    #tailProbCorr = tailProbCorr_MC_2011,
    sfMEtCov = cms.double(1.0), # CV: use 1.0 for Type-1 corrected PFMET, 0.70 for No-PU MET
    power = cms.double(1.0),
    verbosity = cms.int32(0)
)

response_Data_2012 = "8.73728e-01*0.5*(1.0 - TMath::Erf(-9.22702e-03*TMath::Power(x, 1.29731e+00)))" # response of MVA MET algorithm (non-unity response training) in 2012 data
response_MC_2012 = "9.02974e-01*0.5*(1.0 - TMath::Erf(-1.46769e-02*TMath::Power(x, 1.18952e+00)))" # response of MVA MET algorithm (non-unity response training) in Summer'12 MC

svFitEventLikelihoodMEt2a = cms.PSet(
    pluginName = cms.string("svFitEventLikelihoodMEt2a"),
    pluginType = cms.string("SVfitEventLikelihoodMEt2a"),
    response = cms.string(response_MC_2012),
    #tailProbCorr = tailProbCorr_MC_2011,
    power = cms.double(1.0),
    verbosity = cms.int32(0)
)

svFitEventLikelihoodMEt3 = cms.PSet(
    pluginName = cms.string("svFitEventLikelihoodMEt3"),
    pluginType = cms.string("SVfitEventLikelihoodMEt3"),
    srcPFJets = cms.InputTag('ak5PFJets'),
    srcPFCandidates = cms.InputTag('particleFlow'),
    resolution = met_config.METSignificance_params,
    dRoverlapPFJet = cms.double(0.3),
    dRoverlapPFCandidate = cms.double(0.1),
    pfJetPtThreshold = cms.double(20.),
    pfCandPtThreshold = cms.double(5.),    
    numToys = cms.uint32(10000000),    
    power = cms.double(1.0),
    monitorMEtUncertainty = cms.bool(False),
    monitorFilePath = cms.string('/data1/veelken/tmp/svFitStudies/'),
    ##monitorFilePath = cms.string('/tmp/veelken/'),
    verbosity = cms.int32(0)
)

svFitEventBuilder = cms.PSet(
    pluginName = cms.string("svFitEventBuilder"),
    pluginType = cms.string("SVfitEventBuilder"),
    srcBeamSpot = cms.InputTag('offlineBeamSpot'),
    algorithm = cms.string("AdaptiveVertexFitter"),
    applyBeamSpotConstraint = cms.bool(True)
)

svFitConfig_template = cms.PSet(
    event = cms.PSet(
        resonances = cms.PSet(
            A = cms.PSet(
                daughters = cms.PSet(
                    leg1 = cms.PSet(
                        src = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
                        likelihoodFunctions = cms.VPSet(svFitMuonLikelihoodMatrixElement),
                        builder = svFitTauToMuBuilder
                    ),
                    leg2 = cms.PSet(
                        src = cms.InputTag('selectedPatTausForMuTauElectronVetoCumulative'),
                        likelihoodFunctions = cms.VPSet(svFitTauLikelihoodPhaseSpace),
                        builder = svFitTauToHadBuilder
                    )
                ),
                likelihoodFunctions = cms.VPSet(),
                builder = svFitResonanceBuilder
            )
        ),
        srcMEt = cms.InputTag('patPFMETs'),
        srcPrimaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
        likelihoodFunctions = cms.VPSet(svFitEventLikelihoodMEt2),
        builder = svFitEventBuilder
    )
)

svFitProducerByIntegration = cms.EDProducer("SVfitProducerByIntegration",
    config = svFitConfig_template.clone(),
    algorithm = cms.PSet(
        pluginName = cms.string("svFitAlgorithmByIntegration"),
        pluginType = cms.string("SVfitAlgorithmByIntegration"),
        parameters = cms.PSet(
            mass_A = cms.PSet(
                min = cms.double(5.),
                max = cms.double(2000.),                                         
                stepSizeFactor = cms.double(1.025), # nextM = max(stepSizeFactor*currentM, minStepSize)
                minStepSize = cms.double(0.),      
                replace = cms.string("leg1.x"),
                by = cms.string("(A.p4.mass/mass_A)*(A.p4.mass/mass_A)/leg2.x"),
                deltaFuncDerrivative = cms.string("2.*leg1.x/mass_A")                                                   
            )
        ),
        vegasOptions = cms.PSet(
            numCallsGridOpt = cms.uint32(1000),
            numCallsIntEval = cms.uint32(10000),
            maxChi2 = cms.double(2.),
            maxIntEvalIter = cms.uint32(5),                                          
            precision = cms.double(0.00001)
        ),
        max_or_median = cms.string("max"),
        applyJacobiFactors = cms.bool(False),
        verbosity = cms.int32(0)
    ),
    dRmin = cms.double(0.3),
    instanceLabel = cms.string("")
)
svFitProducerByIntegration.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = \
  cms.bool(False)
svFitProducerByIntegration.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = \
  cms.bool(False)

svFitProducerByIntegration2 = cms.EDProducer("SVfitProducerByIntegration2",
    config = svFitConfig_template.clone(),
    algorithm = cms.PSet(
        pluginName = cms.string("svFitAlgorithmByIntegration2"),
        pluginType = cms.string("SVfitAlgorithmByIntegration2"),
        markovChainOptions = cms.PSet(
            mode = cms.string("Metropolis"),
            initMode = cms.string("Gaus"),
            numIterBurnin = cms.uint32(10000),
            numIterSampling = cms.uint32(100000),
            numIterSimAnnealingPhase1 = cms.uint32(2000),
            numIterSimAnnealingPhase2 = cms.uint32(6000),
            T0 = cms.double(15.),
            alpha = cms.double(0.999),
            numChains = cms.uint32(1),
            numBatches = cms.uint32(1),
            L = cms.uint32(1),
            epsilon0 = cms.double(1.e-2),
            nu = cms.double(0.71)
        ),
        max_or_median = cms.string("max"),
        applyJacobiFactors = cms.bool(False),                                          
        verbosity = cms.int32(0)
    ),
    dRmin = cms.double(0.3),
    instanceLabel = cms.string("")
)
svFitProducerByIntegration2.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = \
  cms.bool(False)
svFitProducerByIntegration2.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = \
  cms.bool(False)

svFitProducerByLikelihoodMaximization = cms.EDProducer("SVfitProducer",
    config = svFitConfig_template.clone(),
    algorithm = cms.PSet(
        pluginName = cms.string("svFitAlgorithmByLikelihoodMaximization"),
        pluginType = cms.string("SVfitAlgorithmByLikelihoodMaximization"),
        minimizer  = cms.vstring("Minuit2", "Migrad"),
        maxObjFunctionCalls = cms.uint32(5000),
        verbosity = cms.int32(0)
    ),
    dRmin = cms.double(0.3),
    instanceLabel = cms.string("")
)
svFitProducerByLikelihoodMaximization.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = \
  cms.bool(True)
svFitProducerByLikelihoodMaximization.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = \
  cms.bool(True)
svFitProducerByLikelihoodMaximization.config.event.resonances.A.likelihoodFunctions = cms.VPSet(svFitResonanceLikelihoodLogM)
