# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --mc --eventcontent NANOEDMAODSIM --datatier NANOAODSIM --conditions 126X_mcRun3_2022_realistic_postEE_v1 --step NANO --nThreads 4 --scenario pp --era Run3,run3_nanoAOD_124 --filein /store/group/lpcdisptau/Staus_M_100_100mm_13p6TeV_Run3Summer22EE/MINIAODSIM/230503_144828/0000/TSG-Run3Summer22EEMiniAOD_inMINIAODSIM_1.root --fileout file:nanoAODv11.root
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.Eras.Modifier_run3_common_cff import run3_common

process = cms.Process('NANO',Run3, run3_common)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from LLStaus_Run2.Production.arg_config import *
args = get_args()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(args.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
#     lumisToProcess = cms.untracked.VLuminosityBlockRange('1:30406-1:30406', '1:30408-1:30408')
)
    
if len(args.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = args.lumiFile).getVLuminosityBlockRange()

if args.eventRange != '':
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(',', args.eventRange))

if args.maxEvents > 0:
    process.maxEvents.input = args.maxEvents

process.options = cms.untracked.PSet()


process.options = cms.untracked.PSet(
#     FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
#     SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)



# process.finalJetsPuppi.cut = cms.string('pt > 15')


## Replace AK4 Puppi with AK4 CHS for Run-3
## for CMSSW_14_0_7
use_CHS_jets = True

if use_CHS_jets:
    run3_common.toModify(
        process.linkedObjects, jets="finalJets" ## run 2
    )
    _nanoTableTaskCommonRun3 = process.nanoTableTaskCommon.copy()
    _nanoTableTaskCommonRun3.add(process.jetTask)
    _nanoTableTaskCommonRun3.add(process.jetForMETTask)
    ## remove puppi table otherwise it tries to save the jet table twice
    process.jetPuppiForMETTask.remove(process.corrT1METJetPuppiTable)
    
    process.jetTablesTask.remove(process.bjetNN)
    process.jetTablesTask.remove(process.cjetNN)
    _nanoTableTaskCommonRun3.replace(process.jetPuppiTablesTask, process.jetTablesTask)
    
    process.jetTable.externalVariables = cms.PSet()
    
    run3_common.toReplaceWith(
        process.nanoTableTaskCommon, _nanoTableTaskCommonRun3
    )
    run3_common.toModify(
        process.ptRatioRelForEle, srcJet="updatedJets"
    )
    run3_common.toModify(
        process.ptRatioRelForMu, srcJet="updatedJets"
    )

outputCommands = process.NANOAODSIMEventContent.outputCommands

#### add dis tau
# Output definition
assert(args.disTauTagOutputOpt in [0, 1, 2])

disTauTaggerOnly = False

if args.disTauTagOutputOpt == 0 :
    
    outputCommands = process.NANOAODSIMEventContent.outputCommands

elif args.disTauTagOutputOpt == 1 :
    outputCommands = process.NANOAODSIMEventContent.outputCommands
#     outputCommands = cms.untracked.vstring( 'keep *')
    args.outFile = args.outFile.replace(".root", "_with-disTauTagScore.root")

elif args.disTauTagOutputOpt == 2 :
    
    disTauTaggerOnly = True
    
    outputCommands = cms.untracked.vstring(
        "drop *",
        #"keep *_*_*disTauTag*_*",
        "keep nanoaodFlatTable_jetPuppiTable_*_*",
    )
    
    args.outFile = args.outFile.replace(".root", "_only-disTauTagScore.root")



process.NANOEDMAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(args.outFile),#'file:nanoAODv12.root'),
#     outputCommands = process.NANOAODSIMEventContent.outputCommands,
    outputCommands = outputCommands,
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v5', '')

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOEDMAODSIMoutput_step = cms.EndPath(process.NANOEDMAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOEDMAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
# process.options.numberOfThreads = 4
# process.options.numberOfStreams = 0

from PhysicsTools.NanoAOD.custom_btv_cff import *
def BTVCustomNanoAODStaus(process):
    addPFCands(process,btvNano_switch.btvNano_addallPF_switch,btvNano_switch.btvNano_addAK4_switch,btvNano_switch.btvNano_addAK8_switch)
    
    ### for MC
    process.load("PhysicsTools.NanoAOD.btvMC_cff")
    process.nanoSequenceMC+=ak4onlyPFCandsMCSequence
    return process

process = BTVCustomNanoAODStaus(process)


if use_CHS_jets:
  process.finalJetsAK4Constituents.src = src = cms.InputTag("finalJets")
  process.customAK4ConstituentsTable.jets = cms.InputTag("finalJets")

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon 
#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeCommon(process)

## comment for now
## process.linkedObjects.jets = cms.InputTag("finalJetsPuppi")
from LLStaus_Run2.Production.customize_nanoaod_eventcontent_cff import *
customize_process_and_associate(process, 1, disTauTagOutputOpt = args.disTauTagOutputOpt, useCHSJets = use_CHS_jets)

# if args.disTauTagOutputOpt > 0 :
#     process.schedule.insert(0,process.distau_path)
# End of customisation functions




# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

# if (args.debugEDM) :
#     
#     process.out = cms.OutputModule("PoolOutputModule",
#         fileName = cms.untracked.string("debugEDM.root")
#     )
#     
#     process.output_step = cms.EndPath(process.out)
#     process.schedule.extend([process.output_step])
