import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer
## for 14_0_X
from PhysicsTools.NanoAOD.simplePATJetFlatTableProducer_cfi import simplePATJetFlatTableProducer
from PhysicsTools.NanoAOD.simplePATMuonFlatTableProducer_cfi import simplePATMuonFlatTableProducer

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import *
from PhysicsTools.NanoAOD.jetsAK4_Puppi_cff import *
#from PhysicsTools.NanoAOD.custom_muon_cff import *
#from PhysicsTools.NanoAOD.jets_cff import *


def customize_process_and_associate(process, isMC, disTauTagOutputOpt = 1, useCHSJets = True) :
    # Lost tracks
    process.lostTrackTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("lostTracks"),
        #cut = cms.string(""),
        cut = cms.string("pt > 1"),
        name= cms.string("LostTrack"),
        doc = cms.string("Lost tracks"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            CandVars,
        )
    )

    process.disMuonTable = simplePATMuonFlatTableProducer.clone(
        src = cms.InputTag("slimmedDisplacedMuons"),
        name = cms.string("DisMuon"),
        doc = cms.string("Displaced Muon Collection"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(CandVars,
            ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
            dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
            dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
            dxybs = Var("dB('BS2D')",float,doc="dxy (with sign) wrt the beam spot, in cm",precision=10),
            dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
            dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
            trkChi2 = Var("? globalTrack().isNonnull() ? globalTrack().normalizedChi2() : ? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from either globalTrack or innerTrack "),
            muonHits = Var("? globalTrack().isNonnull() ? globalTrack().hitPattern().numberOfValidMuonHits() : ?  innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidMuonHits() :-99",float,doc="Number of valid Muon Hits from either globalTrack or innerTrack"),
            pixelHits = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidPixelHits() : -99", float, doc="Numbr of valid pixel hits"),
            validFraction = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().validFraction() : -99", float, doc="Inner Track Valid Fraction"),
            positionChi2 = Var("combinedQuality().chi2LocalPosition", float, doc="chi2 Local Position"),
            trkKink = Var("combinedQuality().trkKink", float, doc="Track Kink"),
            ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
            sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
            segmentComp   = Var("segmentCompatibility()", float, doc = "muon segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
            nStations = Var("numberOfMatchedStations", "uint8", doc = "number of matched stations with default arbitration (segment & track)"),
            nTrackerLayers = Var("?track.isNonnull?innerTrack().hitPattern().trackerLayersWithMeasurement():0", "uint8", doc = "number of layers in the tracker"),
            highPurity = Var("?track.isNonnull?innerTrack().quality('highPurity'):0", bool, doc = "inner track is high purity"),
            jetIdx = Var("?hasUserCand('jet')?userCand('jet').key():-1", "int16", doc="index of the associated jet (-1 if none)"),
            svIdx = Var("?hasUserCand('vertex')?userCand('vertex').key():-1", "int16", doc="index of matching secondary vertex"),
            ),
    )

    process.disMuonsMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
        src         = process.disMuonTable.src,                         # final reco collection
        matched     = cms.InputTag("finalGenParticles"),     # final mc-truth particle collection
        mcPdgId     = cms.vint32(13),               # one or more PDG ID (13 = mu); absolute values (see below)
        checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
        mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
        maxDeltaR   = cms.double(0.3),              # Minimum deltaR for the match
        maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
        resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
        resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
    )

    process.disMuonMCTable = cms.EDProducer("CandMCMatchTableProducer",
        src     = process.disMuonTable.src,
        mcMap   = cms.InputTag("disMuonsMCMatchForTable"),
        objName = process.disMuonTable.name,
        objType = cms.string("Muon"), #cms.string("Muon"),
        branchName = cms.string("genPart"),
        docString = cms.string("MC matching to status==1 muons"),
    )
    
    process.disMuonMCTask = cms.Task(process.disMuonsMCMatchForTable, process.disMuonMCTable)  
    process.disMuonTablesTask = cms.Task(process.disMuonTable)

    muonTableForID = muonTable.clone()
    muonTableForID.variables.muonHits = Var("? globalTrack().isNonnull() ? globalTrack().hitPattern().numberOfValidMuonHits() : ?  innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidMuonHits() :-99",float,doc="Number of valid Muon Hits from either globalTrack or innerTrack")
    muonTableForID.variables.pixelHits = Var("? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().hitPattern().numberOfValidPixelHits() : -99", float, doc="Numbr of valid pixel hits")
    muonTableForID.variables.trkChi2 = Var("? globalTrack().isNonnull() ? globalTrack().normalizedChi2() : ? innerTrack().isNonnull() && innerTrack().isAvailable() ? innerTrack().normalizedChi2() : -99",float,doc="Normalized Chi Square from either globalTrack or innerTrack ")
    muonTableForID.variables.positionChi2 = Var("combinedQuality().chi2LocalPosition", float, doc="chi2 Local Position")
    muonTableForID.variables.trkKink = Var("combinedQuality().trkKink", float, doc="Track Kink")
    process.globalReplace("muonTable", muonTableForID)


    # PF candidates
    process.isFromTauForPfCand = cms.EDProducer("IsFromPatTauMapProducer",
        packedPFCandidates = cms.InputTag("packedPFCandidates"),
        #patTaus = cms.InputTag("slimmedTaus"),
        patTaus = cms.InputTag("linkedObjects", "taus"),
        #patTaus = cms.InputTag("selectedPatTaus"),
    )
    
    trk_cond = "hasTrackDetails"
    
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Packed_ParticleFlow_Candidates
    # https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html
    # lostInnerHits: https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html#ab9ef9a12f92e02fa61653ba77ee34274
    # fromPV: https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html#a1e86b4e893738b7cbae410b7f106f339
    process.pfCandTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("packedPFCandidates"),
        #cut = cms.string(""),
        cut = cms.string("pt > 1"),
        name= cms.string("PFCandidate"),
        doc = cms.string("PF candidates"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            CandVars,
            fromPV                  = Var("fromPV"                              , int       , doc = "isolated track comes from PV"),
            lostInnerHits           = Var("lostInnerHits"                       , int       , doc = "Lost inner hits"),
            hasTrackDetails         = Var("hasTrackDetails"                     , bool      , doc = "True if a bestTrack can be extracted from this Candidate"),
            phiAtVtx                = Var("phiAtVtx"                            , float     , doc = "Phi of the candidate's track at the vertex; this is identical to phi() for the vast majority of the particles, but the two might differ for some of them if the calorimeters had contributed significantly in defining the 4-vector of the particle"),
            dxy                     = Var(f"?{trk_cond}?dxy:-999"                , float     , doc = "dxy w.r.t. associated PV"),
            dxyError                = Var(f"?{trk_cond}?dxyError:-999"           , float     , doc = "Error on dxy"),
            dz                      = Var(f"?{trk_cond}?dzAssociatedPV:-999"     , float     , doc = "dz w.r.t. associated PV"),
            dzError                 = Var(f"?{trk_cond}?dzError:-999"            , float     , doc = "Error on dz"),
            vx                      = Var("vx"                                  , float     , doc = "Vertex x"),
            vy                      = Var("vx"                                  , float     , doc = "Vertex y"),
            vz                      = Var("vz"                                  , float     , doc = "Vertex z"),
        ),
        externalVariables = cms.PSet(
            isTauIdxSignalCand     = ExtVar("isFromTauForPfCand:isTauIdxSignalCand"       , int, doc = "Index of the tau if it belongs to pat::Tau::signalCands(); else -1"),
            isTauIdxIsoCand        = ExtVar("isFromTauForPfCand:isTauIdxIsoCand"          , int, doc = "Index of the tau if it belongs to pat::Tau::isolationCands(); else -1"),
            isTauIdxLeadChHadCand  = ExtVar("isFromTauForPfCand:isTauIdxLeadChHadCand"    , int, doc = "Index of the tau if it is pat::Tau::leadChargedHadrCand(); else -1"),
        )
    )
    
    
    # Unfiltered taus
    process.finalTaus.cut = cms.string("pt > 18")
    process.tauTable.doc = cms.string("slimmedTaus after basic selection (" + process.finalTaus.cut.value()+")")
    
    
    # CaloJets // does not work
    process.caloJetTable = cms.EDProducer("SimplePATJetFlatTableProducer",
        src = cms.InputTag("slimmedJets"),
        #cut = cms.string(""),
        cut = cms.string("pt > 15"),
        name= cms.string("Jet"),
        doc = cms.string("slimmedJets"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            P4Vars,
            emEnergyFraction                = Var("emEnergyFraction"                , float),
            #emEnergyInEB                    = Var("emEnergyInEB"                    , float),
            #emEnergyInEE                    = Var("emEnergyInEE"                    , float),
            #emEnergyInHF                    = Var("emEnergyInHF"                    , float),
            energyFractionHadronic          = Var("energyFractionHadronic"          , float),
            #hadEnergyInHB                   = Var("hadEnergyInHB"                   , float),
            #hadEnergyInHE                   = Var("hadEnergyInHE"                   , float),
            #hadEnergyInHF                   = Var("hadEnergyInHF"                   , float),
            #hadEnergyInHO                   = Var("hadEnergyInHO"                   , float),
            #maxEInEmTowers                  = Var("maxEInEmTowers"                  , float),
            #maxEInHadTowers                 = Var("maxEInHadTowers"                 , float),
            #towersArea                      = Var("towersArea"                      , float),
            ## following to be added back
#             detectorP4pt                    = Var("detectorP4.Pt"                   , float),
#             detectorP4eta                   = Var("detectorP4.Eta"                  , float),
#             detectorP4phi                   = Var("detectorP4.Phi"                  , float),
#             detectorP4mass                  = Var("detectorP4.M"                    , float),
#             detectorP4energy                = Var("detectorP4.E"                    , float),
        ),
    )
    
    if isMC:     
        # GenParticles
        myGenParticleTable = genParticleTable.clone()
        myGenParticleTable.variables.vertexX        = Var("vertex.X"      , float)
        myGenParticleTable.variables.vertexY        = Var("vertex.Y"      , float)
        myGenParticleTable.variables.vertexZ        = Var("vertex.Z"      , float)
        myGenParticleTable.variables.vertexRho      = Var("vertex.Rho"    , float)
        myGenParticleTable.variables.vertexR        = Var("vertex.R"      , float)
    
        process.globalReplace("genParticleTable", myGenParticleTable) ## was commented out before JUn 19)
        
    
    if (disTauTagOutputOpt > 0) :
        print("DisTauTagOutputOpt is ", disTauTagOutputOpt)
        process.disTauTag = cms.EDProducer(
            "DisTauTag",
#             graphPath = cms.string("data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"),
            graphPath = cms.string("/afs/cern.ch/work/f/fiorendi/private/displacedTaus/desy/LLStaus_Run2/Production/data/models/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"),
            #jets = cms.InputTag("finalJets"),
            jets = process.jetTable.src,
            pfCandidates = cms.InputTag('packedPFCandidates'),
            save_inputs  = cms.bool(False)
        )
        
        d_disTauTagVars = {
            "disTauTag_score0":     ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
            "disTauTag_score1":     ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
        }
    
    ##process.jetTable.externalVariables = process.jetTable.externalVariables.clone(
    ##    #disTauTag_score0         = ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
    ##    #disTauTag_score1         = ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
    ##    **d_disTauTagVars
    ##)
    

    # Create the task
    if (disTauTagOutputOpt == 0) :
         
        process.custom_nanoaod_task = cms.Task(
            process.lostTrackTable,
            process.isFromTauForPfCand,
#             process.disMuonTablesTask,
#             process.disMuonMCTask,
            process.pfCandTable,
            process.caloJetTable,
        )
    
    elif (disTauTagOutputOpt == 1) :
        
        print ('adding disTau edproducer')
        if useCHSJets:
          process.jetTable.externalVariables = process.jetTable.externalVariables.clone(**d_disTauTagVars)
        ## for puppi jets, use this!
        else:
          process.jetPuppiTable.externalVariables = process.jetPuppiTable.externalVariables.clone(**d_disTauTagVars)
        
        process.custom_nanoaod_task = cms.Task(
#             process.lostTrackTable,
            process.isFromTauForPfCand,
#             process.pfCandTable,
#             process.caloJetTable,
            process.disMuonTablesTask,
            process.disMuonMCTask,
            process.disTauTag,
        )
    
    elif (disTauTagOutputOpt == 2) :
        
        process.jetTable.variables = cms.PSet()
        process.jetTable.externalVariables = cms.PSet(**d_disTauTagVars)
        
        process.custom_nanoaod_task = cms.Task(process.disTauTag)
    
    
    
    # Associate the task to the associate
    process.schedule.associate(process.custom_nanoaod_task)