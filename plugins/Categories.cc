#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/HHAnalysis/interface/Categories.h>

// ***** ***** *****
// Dilepton categories
// ***** ***** *****

const std::vector<HH::Lepton>& DileptonCategory::getLeptons(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>(m_analyzer_name);
    return hh_analyzer.leptons;
}

const std::vector<HH::Dilepton>& DileptonCategory::getDileptons(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>(m_analyzer_name);
    return hh_analyzer.ll;
}

const std::vector<HH::DileptonMetDijet>& DileptonCategory::getDileptonMetDijets(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>(m_analyzer_name);
    return hh_analyzer.llmetjj;
}

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool MuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    return true;
};

bool MuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    return true;
};

void MuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_Mu17_Mu8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
    manager.new_cut("fire_trigger_Mu17_TkMu8", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    manager.new_cut("fire_trigger_IsoMu27", "HLT_IsoMu27_v*");
};

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_Mu8");
        if (path.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_TkMu8");
        if (path.find("HLT_IsoMu27_v") != std::string::npos) manager.pass_cut("fire_trigger_IsoMu27");
    }
}

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    return true;
};

bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    return true;
};

void ElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_Ele17_Ele12", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ*");
    //manager.new_cut("fire_trigger_Ele23_WPLoose", "HLT_Ele23_WPLoose_Gsf_v*");
};

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Ele17_Ele12");
        //if (path.find("HLT_Ele23_WPLoose_Gsf_v") != std::string::npos) manager.pass_cut("fire_trigger_Ele23_WPLoose");
    }
}

// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****
bool ElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    return true;
};

bool ElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    return true;
};

void ElMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_Mu8_Ele17", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_*");
};

void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu8_Ele17");
    }
}

// ***** ***** *****
// Dilepton Mu-El category
// ***** ***** *****
bool MuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    return true;
};

bool MuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    return true;
};

void MuElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_Mu17_Ele12", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_*");
};

void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_Ele12");
    }
}
