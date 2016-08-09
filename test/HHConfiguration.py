
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework
from cp3_llbb.Framework import METProducer

runOnData = True

globalTag_ = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'
processName_ = 'PAT'
if runOnData :
    globalTag_ = '76X_dataRun2_16Dec2015_v0'
    processName_ = 'RECO'

framework = Framework.Framework(runOnData, eras.Run2_25ns, globalTag=globalTag_, processName=processName_)

framework.addAnalyzer('hh_analyzer', cms.PSet(
        type = cms.string('hh_analyzer'),
        prefix = cms.string('hh_'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            # Producers
            electronsProducer = cms.string('electrons'),
            muonsProducer = cms.string('muons'),
            jetsProducer = cms.string('jets'),
            metProducer = cms.string('met'),
            nohfMETProducer = cms.string('nohf_met'),
            # Here are the default value (just to show what is configurable)
            # NB : Isolation and IDs did not change from 74 to 76
            electronIsoCut_EB_Loose = cms.untracked.double(0.0893), # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            electronIsoCut_EE_Loose = cms.untracked.double(0.121),
            electronIsoCut_EB_Tight = cms.untracked.double(0.0354),
            electronIsoCut_EE_Tight = cms.untracked.double(0.0646),
            leadingElectronPtCut = cms.untracked.double(20),
            subleadingElectronPtCut = cms.untracked.double(15),
            electronEtaCut = cms.untracked.double(2.5),
            muonLooseIsoCut = cms.untracked.double(.25), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonTightIsoCut = cms.untracked.double(.15), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            leadingMuonPtCut = cms.untracked.double(20),
            subleadingMuonPtCut = cms.untracked.double(10),
            muonEtaCut = cms.untracked.double(2.4),
            electrons_loose_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
            electrons_medium_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
            electrons_tight_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
            jetEtaCut = cms.untracked.double(2.4),
            jetPtCut = cms.untracked.double(20),
            # BTAG INFO
            discr_name =  cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
            discr_cut_loose =  cms.untracked.double(0.460),
            discr_cut_medium =  cms.untracked.double(0.800),
            discr_cut_tight =  cms.untracked.double(0.935),
            minDR_l_j_Cut = cms.untracked.double(0.3),
            hltDRCut = cms.untracked.double(0.1),
            hltDPtCut = cms.untracked.double(0.5),  # cut will be DPt/Pt < hltDPtCut
            applyBJetRegression = cms.untracked.bool(False), # BE SURE TO ACTIVATE computeRegression FLAG BELOW
            hlt_efficiencies = cms.untracked.PSet(
                Ele17_12Leg1 = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Electron_HLT_Ele17_12Leg1_76X_Tight_HWW.json'),
                Ele17_12Leg2 = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Electron_HLT_Ele17_12Leg1_76X_Tight_HWW.json'),
                DoubleIsoMu17Mu8_IsoMu17leg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Muon_DoubleMu_IsoMu17leg_Run2015D_25ns_PTvsETA_HWW_76X.json'),
                DoubleIsoMu17Mu8_IsoMu8orIsoTkMu8leg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Muon_DoubleMu_IsoMu8orIsoTkMu8leg_Run2015D_25ns_PTvsETA_HWW_76X.json'),
                DoubleIsoMu17Mu8_IsoMu8leg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Muon_DoubleMu_IsoMu8leg_Run2015D_25ns_PTvsETA_HWW_76X.json'),
                #DoubleIsoMu17Mu8_TkMu8leg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Muon_TnP_DoubleIsoMu17Mu8_TkMu8leg_Run2015D_25ns_PTvsETA_binBig_HWW_ID_M_ISO_T.json')
            )
        )
    )
)
# Add PUPPI MET
puppiCfg = cms.PSet(METProducer.default_configuration.clone())
puppiCfg.prefix = cms.string('puppimet_')
puppiCfg.parameters.met = cms.untracked.InputTag('slimmedMETsPuppi')
framework.addProducer('puppimet', puppiCfg)

# Remove fat jets
framework.removeProducer('fat_jets')

#framework.getProducer('jets').parameters.cut = cms.untracked.string("pt > 20")
#framework.getProducer('jets').parameters.computeRegression = cms.untracked.bool(True)

framework.redoJEC()
framework.smearJets()
framework.doSystematics(['jec', 'jer'])
process = framework.create()

if runOnData : 
    process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'
        #'/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/000298CD-87A6-E511-9E56-002590593878.root'
        )
else : 
    process.source.fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall15MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/50FF8034-BEB9-E511-A09C-001EC9ADDD58.root'
        )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
