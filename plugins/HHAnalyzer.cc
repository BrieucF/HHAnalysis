#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>
#include <cp3_llbb/Framework/interface/BTagsAnalyzer.h>
#include <cp3_llbb/HHAnalysis/interface/Categories.h>
#include <cp3_llbb/HHAnalysis/interface/GenStatusFlags.h>

#include <cp3_llbb/Framework/interface/EventProducer.h>
#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/LeptonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cmath>

#define HHANADEBUG 0
#define HH_GEN_DEBUG (false)
#define TT_GEN_DEBUG (false)

void HHAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
    edm::ParameterSet newconfig = edm::ParameterSet(config);
    newconfig.addUntrackedParameter("m_analyzer_name", this->m_name);
    manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons", newconfig);
    manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons", newconfig);
    manager.new_category<ElMuCategory>("elmu", "Category with leading leptons as electron, subleading as muon", newconfig);
    manager.new_category<MuElCategory>("muel", "Category with leading leptons as muon, subleading as electron", newconfig);
}


void HHAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&, const ProducersManager& producers, const AnalyzersManager&, const CategoryManager&) {

    // GEN MATCHING
    const GenParticlesProducer* gp = nullptr;
    if (!event.isRealData())
    {
        // TTBAR MC TRUTH
        gp = &producers.get<GenParticlesProducer>("gen_particles");
        std::function<bool(size_t, size_t)> pruned_decays_from = [&pruned_decays_from, &gp](size_t particle_index, size_t mother_index) -> bool {
            // Iterator over all pruned particles to find if the particle `particle_index` has `mother_index` in its decay history
            if (gp->pruned_mothers_index[particle_index].empty())
                return false;

            size_t index = gp->pruned_mothers_index[particle_index][0];

            if (index == mother_index) {
                return true;
            }

            if (pruned_decays_from(index, mother_index))
                return true;

            return false;
        };

        // 'Pruned' particles are from the hard process
        // 'Packed' particles are stable particles

#if TT_GEN_DEBUG
        std::function<void(size_t)> print_mother_chain = [&gp, &print_mother_chain](size_t p) {

            if (gp->pruned_mothers_index[p].empty()) {
                std::cout << std::endl;
                return;
            }

            size_t index = gp->pruned_mothers_index[p][0];
                std::cout << " <- #" << index << "(" << gp->pruned_pdg_id[index] << ")";
                print_mother_chain(index);
        };
#endif

#define ASSIGN_INDEX( X ) \
        if (flags.isLastCopy()) { \
            gen_##X = i; \
        }\
        if (flags.isFirstCopy()) { \
            gen_##X##_beforeFSR = i; \
        }

    // Assign index to X if it's empty, or Y if not
#define ASSIGN_INDEX2(X, Y, ERROR) \
        if (flags.isLastCopy()) { \
            if (gen_##X == 0) \
                gen_##X = i; \
            else if (gen_##Y == 0)\
                gen_##Y = i; \
            else \
                std::cout << ERROR << std::endl; \
        } \
        if (flags.isFirstCopy()) { \
            if (gen_##X##_beforeFSR == 0) \
                gen_##X##_beforeFSR = i; \
            else if (gen_##Y##_beforeFSR == 0)\
                gen_##Y##_beforeFSR = i; \
            else \
                std::cout << ERROR << std::endl; \
        }

        gen_t = 0; // Index of the top quark
        gen_t_beforeFSR = 0; // Index of the top quark, before any FSR
        gen_tbar = 0; // Index of the anti-top quark
        gen_tbar_beforeFSR = 0; // Index of the anti-top quark, before any FSR

        gen_b = 0; // Index of the b quark coming from the top decay
        gen_b_beforeFSR = 0; // Index of the b quark coming from the top decay, before any FSR
        gen_bbar = 0; // Index of the anti-b quark coming from the anti-top decay
        gen_bbar_beforeFSR = 0; // Index of the anti-b quark coming from the anti-top decay, before any FSR

        gen_jet1_t = 0; // Index of the first jet from the top decay chain
        gen_jet1_t_beforeFSR = 0; // Index of the first jet from the top decay chain, before any FSR
        gen_jet2_t = 0; // Index of the second jet from the top decay chain
        gen_jet2_t_beforeFSR = 0; // Index of the second jet from the top decay chain, before any FSR

        gen_jet1_tbar = 0; // Index of the first jet from the anti-top decay chain
        gen_jet1_tbar_beforeFSR = 0; // Index of the first jet from the anti-top decay chain, before any FSR
        gen_jet2_tbar = 0; // Index of the second jet from the anti-top decay chain
        gen_jet2_tbar_beforeFSR = 0; // Index of the second jet from the anti-top decay chain, before any FSR

        gen_lepton_t = 0; // Index of the lepton from the top decay chain
        gen_lepton_t_beforeFSR = 0; // Index of the lepton from the top decay chain, before any FSR
        gen_neutrino_t = 0; // Index of the neutrino from the top decay chain
        gen_neutrino_t_beforeFSR = 0; // Index of the neutrino from the top decay chain, before any FSR

        gen_lepton_tbar = 0; // Index of the lepton from the anti-top decay chain
        gen_lepton_tbar_beforeFSR = 0; // Index of the lepton from the anti-top decay chain, before any FSR
        gen_neutrino_tbar = 0; // Index of the neutrino from the anti-top decay chain
        gen_neutrino_tbar_beforeFSR = 0; // Index of the neutrino from the anti-top decay chain, before any FSR
        for (size_t i = 0; i < gp->pruned_pdg_id.size(); i++) {

            int16_t pdg_id = gp->pruned_pdg_id[i];
            uint16_t a_pdg_id = std::abs(pdg_id);

            // We only care of particles with PDG id <= 16 (16 is neutrino tau)
            if (a_pdg_id > 16)
                continue;

            GenStatusFlags flags(gp->pruned_status_flags[i]);

            if (! flags.isLastCopy() && ! flags.isFirstCopy())
                continue;

            if (! flags.fromHardProcess())
                continue;

#if TT_GEN_DEBUG
            std::cout << "---" << std::endl;
            std::cout << "Gen particle #" << i << ": PDG id: " << gp->pruned_pdg_id[i];
            print_mother_chain(i);
            flags.dump();
#endif

            if (pdg_id == 6) {
                ASSIGN_INDEX(t);
                continue;
            } else if (pdg_id == -6) {
                ASSIGN_INDEX(tbar);
                continue;
            }

            if (gen_t == 0 || gen_tbar == 0) {
                // Don't bother if we don't have found the tops
                continue;
            }

            bool from_t_decay = pruned_decays_from(i, gen_t);
            bool from_tbar_decay = pruned_decays_from(i, gen_tbar);

            // Only keep particles coming from the tops decay
            if (! from_t_decay && ! from_tbar_decay)
                continue;

            if (pdg_id == 5) {
                // Maybe it's a b coming from the W decay
                if (!flags.isFirstCopy() && flags.isLastCopy() && gen_b == 0) {

                    // This can be a B decaying from a W
                    // However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
                    // Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
                    // If yes, then it's not the B coming directly from the top decay
                    if ((gen_jet1_t_beforeFSR != 0 && std::abs(gp->pruned_pdg_id[gen_jet1_t_beforeFSR]) == 5) ||
                        (gen_jet2_t_beforeFSR != 0 && std::abs(gp->pruned_pdg_id[gen_jet2_t_beforeFSR]) == 5) ||
                        (gen_jet1_tbar_beforeFSR != 0 && std::abs(gp->pruned_pdg_id[gen_jet1_tbar_beforeFSR]) == 5) ||
                        (gen_jet2_tbar_beforeFSR != 0 && std::abs(gp->pruned_pdg_id[gen_jet2_tbar_beforeFSR]) == 5)) {

#if TT_GEN_DEBUG
                        std::cout << "A quark coming from W decay is a b" << std::endl;
#endif

                        if (! (gen_jet1_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_tbar_beforeFSR)) &&
                            ! (gen_jet2_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_tbar_beforeFSR)) &&
                            ! (gen_jet1_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_t_beforeFSR)) &&
                            ! (gen_jet2_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_t_beforeFSR))) {
#if TT_GEN_DEBUG
                            std::cout << "This after-FSR b quark is not coming from a W decay" << std::endl;
#endif
                            gen_b = i;
                            continue;
                        }
#if TT_GEN_DEBUG
                        else {
                            std::cout << "This after-FSR b quark comes from a W decay" << std::endl;
                        }
#endif
                    } else {
#if TT_GEN_DEBUG
                        std::cout << "Assigning gen_b" << std::endl;
#endif
                        gen_b = i;
                        continue;
                    }
                } else if (flags.isFirstCopy() && gen_b_beforeFSR == 0) {
                    gen_b_beforeFSR = i;
                    continue;
                } else {
#if TT_GEN_DEBUG
                    std::cout << "This should not happen!" << std::endl;
#endif
                }
            } else if (pdg_id == -5) {
                if (!flags.isFirstCopy() && flags.isLastCopy() && gen_bbar == 0) {

                    // This can be a B decaying from a W
                    // However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
                    // Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
                    // If yes, then it's not the B coming directly from the top decay
                    if ((gen_jet1_t_beforeFSR != 0 && std::abs(gp->pruned_pdg_id[gen_jet1_t_beforeFSR]) == 5) ||
                        (gen_jet2_t_beforeFSR != 0 && std::abs(gp->pruned_pdg_id[gen_jet2_t_beforeFSR]) == 5) ||
                        (gen_jet1_tbar_beforeFSR != 0 && std::abs(gp->pruned_pdg_id[gen_jet1_tbar_beforeFSR]) == 5) ||
                        (gen_jet2_tbar_beforeFSR != 0 && std::abs(gp->pruned_pdg_id[gen_jet2_tbar_beforeFSR]) == 5)) {

#if TT_GEN_DEBUG
                        std::cout << "A quark coming from W decay is a bbar" << std::endl;
#endif

                        if (! (gen_jet1_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_tbar_beforeFSR)) &&
                            ! (gen_jet2_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_tbar_beforeFSR)) &&
                            ! (gen_jet1_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_t_beforeFSR)) &&
                            ! (gen_jet2_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_t_beforeFSR))) {
#if TT_GEN_DEBUG
                            std::cout << "This after-fsr b anti-quark is not coming from a W decay" << std::endl;
#endif
                            gen_bbar = i;
                            continue;
                        }
#if TT_GEN_DEBUG
                        else {
                            std::cout << "This after-fsr b anti-quark comes from a W decay" << std::endl;
                        }
#endif
                    } else {
#if TT_GEN_DEBUG
                        std::cout << "Assigning gen_bbar" << std::endl;
#endif
                        gen_bbar = i;
                        continue;
                    }
                } else if (flags.isFirstCopy() && gen_bbar_beforeFSR == 0) {
                    gen_bbar_beforeFSR = i;
                    continue;
                }
            }

            if ((gen_tbar == 0) || (gen_t == 0))
                continue;

            if (gen_t != 0 && from_t_decay) {
#if TT_GEN_DEBUG
            std::cout << "Coming from the top chain decay" << std::endl;
#endif
                if (a_pdg_id >= 1 && a_pdg_id <= 5) {
                    ASSIGN_INDEX2(jet1_t, jet2_t, "Error: more than two quarks coming from top decay");
                } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
                    ASSIGN_INDEX(lepton_t);
                } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
                    ASSIGN_INDEX(neutrino_t);
                } else {
                    std::cout << "Error: unknown particle coming from top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
                }
            } else if (gen_tbar != 0 && from_tbar_decay) {
#if TT_GEN_DEBUG
            std::cout << "Coming from the anti-top chain decay" << std::endl;
#endif
                if (a_pdg_id >= 1 && a_pdg_id <= 5) {
                    ASSIGN_INDEX2(jet1_tbar, jet2_tbar, "Error: more than two quarks coming from anti-top decay");
                } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
                    ASSIGN_INDEX(lepton_tbar);
                } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
                    ASSIGN_INDEX(neutrino_tbar);
                } else {
                    std::cout << "Error: unknown particle coming from anti-top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
                }
            }
        }

        if (!gen_t || !gen_tbar) {
#if TT_GEN_DEBUG
            std::cout << "This is not a ttbar event" << std::endl;
#endif
            gen_ttbar_decay_type = NotTT;
        }
        if (gen_ttbar_decay_type != NotTT) {

        if ((gen_jet1_t != 0) && (gen_jet2_t != 0) && (gen_jet1_tbar != 0) && (gen_jet2_tbar != 0)) {
#if TT_GEN_DEBUG
            std::cout << "Hadronic ttbar decay" << std::endl;
#endif
            gen_ttbar_decay_type = Hadronic;
        } else if (
                ((gen_lepton_t != 0) && (gen_lepton_tbar == 0)) ||
                ((gen_lepton_t == 0) && (gen_lepton_tbar != 0))
                ) {

#if TT_GEN_DEBUG
            std::cout << "Semileptonic ttbar decay" << std::endl;
#endif

            uint16_t lepton_pdg_id;
            if (gen_lepton_t != 0)
                lepton_pdg_id = std::abs(gp->pruned_pdg_id[gen_lepton_t]);
            else
                lepton_pdg_id = std::abs(gp->pruned_pdg_id[gen_lepton_tbar]);

            if (lepton_pdg_id == 11)
                gen_ttbar_decay_type = Semileptonic_e;
            else if (lepton_pdg_id == 13)
                gen_ttbar_decay_type = Semileptonic_mu;
            else
                gen_ttbar_decay_type = Semileptonic_tau;
        } else if (gen_lepton_t != 0 && gen_lepton_tbar != 0) {
            uint16_t lepton_t_pdg_id = std::abs(gp->pruned_pdg_id[gen_lepton_t]);
            uint16_t lepton_tbar_pdg_id = std::abs(gp->pruned_pdg_id[gen_lepton_tbar]);

#if TT_GEN_DEBUG
            std::cout << "Dileptonic ttbar decay" << std::endl;
#endif

            if (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 11)
                gen_ttbar_decay_type = Dileptonic_ee;
            else if (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 13)
                gen_ttbar_decay_type = Dileptonic_mumu;
            else if (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 15)
                gen_ttbar_decay_type = Dileptonic_tautau;
            else if (
                    (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 13) ||
                    (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 11)
                    ) {
                gen_ttbar_decay_type = Dileptonic_mue;
            }
            else if (
                    (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 15) ||
                    (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 11)
                    ) {
                gen_ttbar_decay_type = Dileptonic_etau;
            }
            else if (
                    (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 15) ||
                    (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 13)
                    ) {
                gen_ttbar_decay_type = Dileptonic_mutau;
            } else {
                std::cout << "Error: unknown dileptonic ttbar decay." << std::endl;
                gen_ttbar_decay_type = NotTT;
            }
        } else {
            std::cout << "Error: unknown ttbar decay." << std::endl;
            gen_ttbar_decay_type = UnknownTT;
        }
        // Retrieve the 4-vector of the hard process particles for transfer function
        if (gen_b_beforeFSR != 0)
            gen_ttbar_b_beforeFSR_p4 = gp->pruned_p4[gen_b_beforeFSR];
        if (gen_bbar_beforeFSR != 0)
            gen_ttbar_bbar_beforeFSR_p4 = gp->pruned_p4[gen_bbar_beforeFSR]; 
        if (gen_lepton_t_beforeFSR != 0)
            gen_ttbar_lepton_t_beforeFSR_p4 = gp->pruned_p4[gen_lepton_t_beforeFSR];
        if (gen_lepton_tbar_beforeFSR != 0)
            gen_ttbar_lepton_tbar_beforeFSR_p4 = gp->pruned_p4[gen_lepton_tbar_beforeFSR];
        }
        } // end of if isMC

    //float mh = event.isRealData() ? 125.02 : 125.0;
    LorentzVector null_p4(0., 0., 0., 0.);
    const EventProducer& fwevent = producers.get<EventProducer>("event");
    float event_weight = fwevent.weight;
    float tmp_count_has2leptons = 0.;
    float tmp_count_has2leptons_elel = 0.;
    float tmp_count_has2leptons_elmu = 0.;
    float tmp_count_has2leptons_muel = 0.;
    float tmp_count_has2leptons_mumu = 0.;
    float tmp_count_has2leptons_1llmetjj = 0.;
    float tmp_count_has2leptons_elel_1llmetjj = 0.;
    float tmp_count_has2leptons_elmu_1llmetjj = 0.;
    float tmp_count_has2leptons_muel_1llmetjj = 0.;
    float tmp_count_has2leptons_mumu_1llmetjj = 0.;
    float tmp_count_has2leptons_1llmetjj_2btagM = 0.;
    float tmp_count_has2leptons_elel_1llmetjj_2btagM = 0.;
    float tmp_count_has2leptons_elmu_1llmetjj_2btagM = 0.;
    float tmp_count_has2leptons_muel_1llmetjj_2btagM = 0.;
    float tmp_count_has2leptons_mumu_1llmetjj_2btagM = 0.;

    // ***** ***** *****
    // Trigger Matching
    // ***** ***** *****

    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    //Function that tries to match `lepton` with an online object, using a deltaR and a deltaPt cut   
    //Returns the index inside the HLTProducer collection, or -1 if no match is found.
    //(Taken from https://github.com/blinkseb/TTAnalysis/blob/c2a2d5de3e4281943c19c582afb452b8ef6457f1/plugins/TTAnalyzer.cc#L533)
    auto matchOfflineLepton = [&](HH::Lepton& lepton) {

        if (lepton.hlt_already_tried_matching)
            return lepton.hlt_idx;
        float min_dr = std::numeric_limits<float>::max();
        float final_dpt_over_pt = std::numeric_limits<float>::max();

        int8_t index = -1;
        for (size_t hlt_object = 0; hlt_object < hlt.object_p4.size(); hlt_object++) {

            float dr = ROOT::Math::VectorUtil::DeltaR(lepton.p4, hlt.object_p4[hlt_object]);
            float dpt_over_pt = std::abs(lepton.p4.Pt() - hlt.object_p4[hlt_object].Pt()) / lepton.p4.Pt();
            if (dr < min_dr && dpt_over_pt < final_dpt_over_pt) {
                min_dr = dr;
                final_dpt_over_pt = dpt_over_pt;
                if (dr < m_hltDRCut && dpt_over_pt < m_hltDPtCut)
                    index = hlt_object;
            }
        }
        lepton.hlt_idx = index;
        lepton.hlt_already_tried_matching = true;
        lepton.hlt_DR_matchedObject = min_dr;
        lepton.hlt_DPtOverPt_matchedObject = final_dpt_over_pt;
        return index;
    };

    // ********** 
    // Leptons and dileptons
    // ********** 
    const ElectronsProducer& allelectrons = producers.get<ElectronsProducer>(m_electrons_producer);
    const MuonsProducer& allmuons = producers.get<MuonsProducer>(m_muons_producer);
    float tt_genLepMatching_maxDR = 0.1;
    float tt_genJetMatching_maxDR = 0.2;

    nMuonsL = 0;
    nMuonsT = 0;
    nMuonsHWW = 0;
    nElectronsL = 0;
    nElectronsT = 0;
    nElectronsHWW = 0;
    nLeptonsL = 0;
    nLeptonsT = 0;
    nLeptonsHWW = 0;

    leptons.clear();
    ll.clear();

    // Electrons
    for (unsigned int ielectron = 0; ielectron < allelectrons.p4.size(); ielectron++)
    {
        if ( !(fabs(allelectrons.p4[ielectron].Eta()) < m_electronEtaCut) )
           continue;

        HH::Lepton ele;
        ele.passLep1PtCut = allelectrons.p4[ielectron].Pt() > m_leadingElectronPtCut;
        ele.passLep2PtCut = allelectrons.p4[ielectron].Pt() > m_subleadingElectronPtCut;
        ele.id_L = allelectrons.ids[ielectron][m_electron_loose_wp_name];
        ele.id_M = allelectrons.ids[ielectron][m_electron_medium_wp_name];
        ele.id_T = allelectrons.ids[ielectron][m_electron_tight_wp_name];
        if (allelectrons.isEB[ielectron])
            ele.id_HWW = allelectrons.dz[ielectron] < 0.373 && allelectrons.dxy[ielectron] < 0.1;
        else
            ele.id_HWW = allelectrons.dz[ielectron] < 0.602 && allelectrons.dxy[ielectron] < 0.2;
        ele.iso_L = allelectrons.isEB[ielectron] ? (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EB_Loose) : (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EE_Loose);
        ele.iso_T = allelectrons.isEB[ielectron] ? (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EB_Tight) : (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EE_Tight);
        ele.iso_HWW = ele.iso_T && allelectrons.relEcalPFClusterIso[ielectron] < 0.45 && allelectrons.relHcalPFClusterIso[ielectron] < 0.25 && allelectrons.relTrackIso[ielectron] < 0.2;
        if (ele.id_L && ele.iso_L){
            nElectronsL++;
            nLeptonsL++;
        }
        if (ele.id_T && ele.iso_T){
            nElectronsT++;
            nLeptonsT++;
        }
        if (ele.id_HWW && ele.iso_HWW){
            nElectronsHWW++;
            nLeptonsHWW++;
        }
        ele.p4 = allelectrons.p4[ielectron];
        if (!hlt.paths.empty())
           matchOfflineLepton(ele);
        ele.charge = allelectrons.charge[ielectron];
        ele.idx = ielectron;
        ele.isMu = false;
        ele.isEl = true;
        if (!(ele.iso_HWW && ele.id_HWW))
            continue;
        // GEN MATCHING
        ele.tt_matched = false;
        ele.tt_gen_p4 = null_p4;
        ele.tt_gen_DR = std::numeric_limits<float>::max();
        ele.tt_parton_p4 = null_p4;
        ele.tt_parton_DR = std::numeric_limits<float>::max();
        if (!event.isRealData() && gen_lepton_t != 0 && gen_lepton_t_beforeFSR != 0 && gen_neutrino_tbar !=0 && gen_lepton_tbar_beforeFSR != 0) {
            float temp_DR_lep_t = std::numeric_limits<float>::max(); 
            float temp_DR_lep_antit = std::numeric_limits<float>::max(); 
            float DR_lepParton = std::numeric_limits<float>::max(); 
            uint16_t parton_id = 0;
            uint16_t parton_id_beforeFSR = 0;
            bool goodID = false;
            if (std::abs(gp->pruned_pdg_id[gen_lepton_t]) == 11) { 
                temp_DR_lep_t = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[gen_lepton_t], ele.p4);
                goodID = true;
            }
            if (std::abs(gp->pruned_pdg_id[gen_lepton_tbar]) == 11) { 
                temp_DR_lep_antit = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[gen_lepton_tbar], ele.p4);
                goodID = true;
            }
            if (goodID) {
                parton_id = (temp_DR_lep_t < temp_DR_lep_antit) ? gen_lepton_t : gen_lepton_tbar;
                DR_lepParton = (temp_DR_lep_t < temp_DR_lep_antit) ? temp_DR_lep_t : temp_DR_lep_antit;
                parton_id_beforeFSR = (temp_DR_lep_t < temp_DR_lep_antit) ? gen_lepton_t_beforeFSR : gen_lepton_tbar_beforeFSR;
                ele.tt_gen_DR = DR_lepParton;
                ele.tt_gen_p4 = gp->pruned_p4[parton_id];
                ele.tt_parton_p4 = gp->pruned_p4[parton_id_beforeFSR];
                ele.tt_parton_DR = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[parton_id_beforeFSR], ele.p4);
                if (ele.tt_gen_DR < tt_genLepMatching_maxDR)
                    ele.tt_matched = true;
            }
        }
        ele.gen_matched = allelectrons.matched[ielectron];
        ele.gen_p4 = ele.gen_matched ? allelectrons.gen_p4[ielectron] : null_p4;
        ele.gen_DR = ele.gen_matched ? ROOT::Math::VectorUtil::DeltaR(ele.p4, ele.gen_p4): -1.;
        ele.gen_DPtOverPt = ele.gen_matched ? (ele.p4.Pt() - ele.gen_p4.Pt()) / ele.p4.Pt() : -10.;
        electrons.push_back(ielectron);
        leptons.push_back(ele);
    }//end of loop on electrons
    nElectrons = electrons.size();

    // Muons
    for (unsigned int imuon = 0; imuon < allmuons.p4.size(); imuon++)
    {
        if (!(fabs(allmuons.p4[imuon].Eta()) < m_muonEtaCut))
            continue;
        HH::Lepton mu;
        mu.passLep1PtCut = allmuons.p4[imuon].Pt() > m_leadingMuonPtCut; 
        mu.passLep2PtCut = allmuons.p4[imuon].Pt() > m_subleadingMuonPtCut; 
        mu.id_L = allmuons.isLoose[imuon];
        mu.id_M = allmuons.isMedium[imuon];
        mu.id_T = allmuons.isTight[imuon];
        mu.id_HWW = mu.id_M && (mu.p4.Pt() < 20. ? fabs(allmuons.dxy[imuon]) < 0.01 : fabs(allmuons.dxy[imuon]) < 0.02) && (fabs(allmuons.dz[imuon]) < 0.1);
        mu.iso_L = allmuons.relativeIsoR04_deltaBeta[imuon] < m_muonLooseIsoCut;
        mu.iso_T = allmuons.relativeIsoR04_deltaBeta[imuon] < m_muonTightIsoCut;
        mu.iso_HWW = mu.iso_T; // For the isolation use relative PF isolation (cone size = 0.4) with deltaBeta PU corrections (a.l.a. Run I) and WP < 0.15
        if (mu.id_L && mu.iso_L){
            nMuonsL++;
            nLeptonsL++;
        }
        if (mu.id_T && mu.iso_T){
            nMuonsT++;
            nLeptonsT++;
        }
        if (mu.id_HWW && mu.iso_HWW){
            nMuonsHWW++;
            nLeptonsHWW++;
        }
        if (!(mu.iso_HWW && mu.id_HWW))
            continue;
        mu.p4 = allmuons.p4[imuon];
        if (!hlt.paths.empty())
           matchOfflineLepton(mu);
        mu.charge = allmuons.charge[imuon];
        mu.idx = imuon;
        mu.isMu = true;
        mu.isEl = false;
        // GEN MATCHING
        mu.tt_matched = false;
        mu.tt_gen_p4 = null_p4;
        mu.tt_gen_DR = std::numeric_limits<float>::max();
        mu.tt_parton_p4 = null_p4;
        mu.tt_parton_DR = std::numeric_limits<float>::max();
        if (!event.isRealData() && gen_lepton_t != 0 && gen_lepton_t_beforeFSR != 0 && gen_neutrino_tbar !=0 && gen_lepton_tbar_beforeFSR != 0) {
            float temp_DR_lep_t = std::numeric_limits<float>::max(); 
            float temp_DR_lep_antit = std::numeric_limits<float>::max(); 
            float DR_lepParton = std::numeric_limits<float>::max(); 
            uint16_t parton_id = 0;
            uint16_t parton_id_beforeFSR = 0;
            bool goodID = false;
            if (std::abs(gp->pruned_pdg_id[gen_lepton_t]) == 13) { 
                temp_DR_lep_t = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[gen_lepton_t], mu.p4);
                goodID = true;
            }
            if (std::abs(gp->pruned_pdg_id[gen_lepton_tbar]) == 13) { 
                temp_DR_lep_antit = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[gen_lepton_tbar], mu.p4);
                goodID = true;
            }
            if (goodID) {
                parton_id = (temp_DR_lep_t < temp_DR_lep_antit) ? gen_lepton_t : gen_lepton_tbar;
                DR_lepParton = (temp_DR_lep_t < temp_DR_lep_antit) ? temp_DR_lep_t : temp_DR_lep_antit;
                parton_id_beforeFSR = (temp_DR_lep_t < temp_DR_lep_antit) ? gen_lepton_t_beforeFSR : gen_lepton_tbar_beforeFSR;
                mu.tt_gen_DR = DR_lepParton;
                mu.tt_gen_p4 = gp->pruned_p4[parton_id];
                mu.tt_parton_p4 = gp->pruned_p4[parton_id_beforeFSR];
                mu.tt_parton_DR = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[parton_id_beforeFSR], mu.p4);
                if (mu.tt_gen_DR < tt_genLepMatching_maxDR)
                    mu.tt_matched = true;
            }
        }
        mu.gen_matched = allmuons.matched[imuon];
        mu.gen_p4 = mu.gen_matched ? allmuons.gen_p4[imuon] : null_p4;
        mu.gen_DR = mu.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mu.p4, mu.gen_p4): -1.;
        mu.gen_DPtOverPt = mu.gen_matched ? (mu.p4.Pt() - mu.gen_p4.Pt()) / mu.p4.Pt() : -10.;
        muons.push_back(imuon);
        leptons.push_back(mu);
    }//end of loop on muons

    nMuons = muons.size();
    nLeptons = leptons.size();

    // sort leptons by pt (ignoring flavour, id and iso)
    std::sort(leptons.begin(), leptons.end(), [](const HH::Lepton& lep1, const HH::Lepton& lep2) { return lep1.p4.Pt() > lep2.p4.Pt(); });     

    // Dileptons
    for (unsigned int ilep1 = 0; ilep1 < leptons.size(); ilep1++)
    {
        for (unsigned int ilep2 = ilep1+1; ilep2 < leptons.size(); ilep2++)
        {
            HH::Dilepton dilep;
            dilep.p4 = leptons[ilep1].p4 + leptons[ilep2].p4;
            dilep.idxs = std::make_pair(leptons[ilep1].idx, leptons[ilep2].idx);
            dilep.ilep1 = ilep1;
            dilep.ilep2 = ilep2;
            dilep.isOS = leptons[ilep1].charge * leptons[ilep2].charge < 0;
            dilep.isMuMu = leptons[ilep1].isMu && leptons[ilep2].isMu;
            dilep.isElEl = leptons[ilep1].isEl && leptons[ilep2].isEl;
            dilep.isElMu = leptons[ilep1].isEl && leptons[ilep2].isMu;
            dilep.isMuEl = leptons[ilep1].isMu && leptons[ilep2].isEl;
            dilep.isSF = dilep.isMuMu || dilep.isElEl;
            dilep.id_LL = leptons[ilep1].id_L && leptons[ilep2].id_L;
            dilep.id_LM = (leptons[ilep1].id_L && leptons[ilep2].id_M) || (leptons[ilep2].id_L && leptons[ilep1].id_M);
            dilep.id_LT = (leptons[ilep1].id_L && leptons[ilep2].id_T) || (leptons[ilep2].id_L && leptons[ilep1].id_T);
            dilep.id_LHWW = (leptons[ilep1].id_L && leptons[ilep2].id_HWW) || (leptons[ilep2].id_L && leptons[ilep1].id_HWW);
            dilep.id_ML = (leptons[ilep1].id_M && leptons[ilep2].id_L) || (leptons[ilep2].id_M && leptons[ilep1].id_L);
            dilep.id_MM = leptons[ilep1].id_M && leptons[ilep2].id_M;
            dilep.id_MT = (leptons[ilep1].id_T && leptons[ilep2].id_M) || (leptons[ilep2].id_T && leptons[ilep1].id_M);
            dilep.id_MHWW = (leptons[ilep1].id_M && leptons[ilep2].id_HWW) || (leptons[ilep2].id_M && leptons[ilep1].id_HWW);
            dilep.id_TL = (leptons[ilep1].id_T && leptons[ilep2].id_L) || (leptons[ilep2].id_T && leptons[ilep1].id_L);
            dilep.id_TM = (leptons[ilep1].id_T && leptons[ilep2].id_M) || (leptons[ilep2].id_T && leptons[ilep1].id_M);
            dilep.id_TT = leptons[ilep1].id_T && leptons[ilep2].id_T;
            dilep.id_THWW = (leptons[ilep1].id_T && leptons[ilep2].id_HWW) || (leptons[ilep2].id_T && leptons[ilep1].id_HWW);
            dilep.id_HWWL = (leptons[ilep1].id_HWW && leptons[ilep2].id_L) || (leptons[ilep2].id_HWW && leptons[ilep1].id_L);
            dilep.id_HWWM = (leptons[ilep1].id_HWW && leptons[ilep2].id_M) || (leptons[ilep2].id_HWW && leptons[ilep1].id_M);
            dilep.id_HWWT = (leptons[ilep1].id_HWW && leptons[ilep2].id_T) || (leptons[ilep2].id_HWW && leptons[ilep1].id_T);
            dilep.id_HWWHWW = leptons[ilep1].id_HWW && leptons[ilep2].id_HWW;
            dilep.iso_LL = leptons[ilep1].iso_L && leptons[ilep2].iso_L;
            dilep.iso_LT = (leptons[ilep1].iso_L && leptons[ilep2].iso_T) || (leptons[ilep2].iso_L && leptons[ilep1].iso_T);
            dilep.iso_LHWW = (leptons[ilep1].iso_L && leptons[ilep2].iso_HWW) || (leptons[ilep2].iso_L && leptons[ilep1].iso_HWW);
            dilep.iso_TL = (leptons[ilep1].iso_T && leptons[ilep2].iso_L) || (leptons[ilep2].iso_T && leptons[ilep1].iso_L);
            dilep.iso_TT = leptons[ilep1].iso_T && leptons[ilep2].iso_T;
            dilep.iso_THWW = (leptons[ilep1].iso_T && leptons[ilep2].iso_HWW) || (leptons[ilep2].iso_T && leptons[ilep1].iso_HWW);
            dilep.iso_HWWL = (leptons[ilep1].iso_HWW && leptons[ilep2].iso_L) || (leptons[ilep2].iso_HWW && leptons[ilep1].iso_L);
            dilep.iso_HWWT = (leptons[ilep1].iso_HWW && leptons[ilep2].iso_T) || (leptons[ilep2].iso_HWW && leptons[ilep1].iso_T);
            dilep.iso_HWWHWW = leptons[ilep1].iso_HWW && leptons[ilep2].iso_HWW;
            dilep.DR_l_l = ROOT::Math::VectorUtil::DeltaR(leptons[ilep1].p4, leptons[ilep2].p4);
            dilep.DPhi_l_l = fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep1].p4, leptons[ilep2].p4));
            dilep.ht_l_l = leptons[ilep1].p4.Pt() + leptons[ilep2].p4.Pt();
            dilep.hlt_idxs = std::make_pair(leptons[ilep1].hlt_idx, leptons[ilep2].hlt_idx);
            dilep.gen_matched = leptons[ilep1].gen_matched && leptons[ilep2].gen_matched;
            dilep.gen_p4 = dilep.gen_matched ? leptons[ilep1].gen_p4 + leptons[ilep2].gen_p4 : null_p4;
            dilep.gen_DR = dilep.gen_matched ? ROOT::Math::VectorUtil::DeltaR(dilep.p4, dilep.gen_p4) : -1.;
            dilep.gen_DPtOverPt = dilep.gen_matched ? (dilep.p4.Pt() - dilep.gen_p4.Pt()) / dilep.p4.Pt() : -10.;
            if (event.isRealData()) {
               dilep.trigger_efficiency = 1.;
               dilep.trigger_efficiency_downVariated = 1.;
               dilep.trigger_efficiency_downVariated_Arun = 1.;
               dilep.trigger_efficiency_upVariated = 1.;
               dilep.trigger_efficiency_upVariated_Arun = 1.;
            }
            else {
               fillTriggerEfficiencies(leptons[ilep1], leptons[ilep2], dilep);
            }

            //if (!dilep.isOS)
            //    continue;

            // Counters
            tmp_count_has2leptons = event_weight;
            if (dilep.isElEl)
                tmp_count_has2leptons_elel = event_weight;
            if (dilep.isElMu)
                tmp_count_has2leptons_elmu = event_weight;
            if (dilep.isElEl)
                tmp_count_has2leptons_muel = event_weight;
            if (dilep.isElEl)
                tmp_count_has2leptons_mumu = event_weight;
            // Fill
            ll.push_back(dilep); 
        }
    }

    std::sort(ll.begin(), ll.end(), [](const HH::Dilepton& ll1, const HH::Dilepton& ll2) { return ll1.ht_l_l > ll2.ht_l_l; });     

    // ***** 
    // Adding MET(s)
    // ***** 
    met.clear();
    llmet.clear();
    const METProducer& pf_met = producers.get<METProducer>(m_met_producer);
    HH::Met mymet;
    mymet.p4 = pf_met.p4;
    mymet.isNoHF = false;
    mymet.gen_matched = false;
    mymet.gen_p4 = null_p4;
    mymet.gen_DR = -1.;
    mymet.gen_DPhi = -1.;
    mymet.gen_DPtOverPt = -10.;
    if (!event.isRealData())
    { // genMet is not constructed in the framework, so construct it manually out of the neutrinos hanging around the mc particles
        for (unsigned int ip = 0; ip < gp->pruned_p4.size(); ip++) {
            std::bitset<15> flags (gp->pruned_status_flags[ip]);
            if (!flags.test(13)) continue; // take the last copies
            if (abs(gp->pruned_pdg_id[ip]) == 12 || abs(gp->pruned_pdg_id[ip]) == 14 || abs(gp->pruned_pdg_id[ip]) == 16)
            {
                mymet.gen_matched = true;
                mymet.gen_p4 += gp->pruned_p4[ip];
            }
        }
        mymet.gen_DR = mymet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mymet.p4, mymet.gen_p4) : -1.;
        mymet.gen_DPhi = mymet.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(mymet.p4, mymet.gen_p4)) : -1.;
        mymet.gen_DPtOverPt = mymet.gen_matched ? (mymet.p4.Pt() - mymet.gen_p4.Pt()) / mymet.p4.Pt() : -10.;
    }
    met.push_back(mymet);
    //const METProducer& nohf_met = producers.get<METProducer>(m_nohf_met_producer);  // so that nohfmet is available in the tree
    //const METProducer& puppi_met = producers.get<METProducer>("puppimet");
    // TODO: adding puppi met will require changing the Met AND DileptonMet struct

    for (unsigned int imet = 0; imet < met.size(); imet++)
    {
        for (unsigned int ill = 0; ill < ll.size(); ill++)
        {
            HH::DileptonMet myllmet;
            // DileptonMet inherits from Dilepton struct, initalize everything properly
            myllmet.p4 = ll[ill].p4 + met[imet].p4;
            // blind copy of the ll content
            myllmet.idxs = std::make_pair(ll[ill].idxs.first, ll[ill].idxs.second);
            myllmet.ilep1 = ll[ill].ilep1;
            myllmet.ilep2 = ll[ill].ilep2;
            myllmet.isOS = ll[ill].isOS;
            myllmet.isMuMu = ll[ill].isMuMu;
            myllmet.isElEl = ll[ill].isElEl;
            myllmet.isElMu = ll[ill].isElMu;
            myllmet.isMuEl = ll[ill].isMuEl;
            myllmet.isSF = ll[ill].isSF;
            myllmet.id_LL = ll[ill].id_LL;
            myllmet.id_LM = ll[ill].id_LM;
            myllmet.id_LT = ll[ill].id_LT;
            myllmet.id_LHWW = ll[ill].id_LHWW;
            myllmet.id_ML = ll[ill].id_ML;
            myllmet.id_MM = ll[ill].id_MM;
            myllmet.id_MT = ll[ill].id_MT;
            myllmet.id_MHWW = ll[ill].id_MHWW;
            myllmet.id_TL = ll[ill].id_TL;
            myllmet.id_TM = ll[ill].id_TM;
            myllmet.id_TT = ll[ill].id_TT;
            myllmet.id_THWW = ll[ill].id_THWW;
            myllmet.id_HWWL = ll[ill].id_HWWL;
            myllmet.id_HWWM = ll[ill].id_HWWM;
            myllmet.id_HWWT = ll[ill].id_HWWT;
            myllmet.id_HWWHWW = ll[ill].id_HWWHWW;
            myllmet.iso_LL = ll[ill].iso_LL;
            myllmet.iso_LT = ll[ill].iso_LT;
            myllmet.iso_LHWW = ll[ill].iso_LHWW;
            myllmet.iso_TL = ll[ill].iso_TL;
            myllmet.iso_TT = ll[ill].iso_TT;
            myllmet.iso_THWW = ll[ill].iso_THWW;
            myllmet.iso_HWWL = ll[ill].iso_HWWL;
            myllmet.iso_HWWT = ll[ill].iso_HWWT;
            myllmet.iso_HWWHWW = ll[ill].iso_HWWHWW;
            myllmet.DR_l_l = ll[ill].DR_l_l;
            myllmet.DPhi_l_l = ll[ill].DPhi_l_l;
            // content specific to HH:DileptonMet
            myllmet.ill = ill;
            myllmet.imet = imet;
            myllmet.isNoHF = met[imet].isNoHF;
            float dphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(ll[ill].p4, met[imet].p4));
            myllmet.DPhi_ll_met = dphi;
            float mindphi = std::min(fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ll[ill].idxs.first].p4, met[imet].p4)), fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ll[ill].idxs.second].p4, met[imet].p4)));
            myllmet.minDPhi_l_met = mindphi; 
            float maxdphi = std::max(fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ll[ill].idxs.first].p4, met[imet].p4)), fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ll[ill].idxs.second].p4, met[imet].p4)));
            myllmet.maxDPhi_l_met = maxdphi;
            myllmet.MT = (ll[ill].p4 + met[imet].p4).M();
            myllmet.MT_formula = std::sqrt(2 * ll[ill].p4.Pt() * met[imet].p4.Pt() * (1-std::cos(dphi)));
            myllmet.projectedMet = mindphi >= M_PI ? met[imet].p4.Pt() : met[imet].p4.Pt() * std::sin(mindphi);
            myllmet.gen_matched = ll[ill].gen_matched && met[imet].gen_matched;
            myllmet.gen_p4 = myllmet.gen_matched ? ll[ill].gen_p4 + met[imet].gen_p4 : null_p4;
            myllmet.gen_DR = myllmet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myllmet.p4, myllmet.gen_p4) : -1.;
            myllmet.gen_DPhi = myllmet.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(myllmet.p4, myllmet.gen_p4)) : -1.;
            myllmet.gen_DPtOverPt = myllmet.gen_matched ? (myllmet.p4.Pt() - myllmet.gen_p4.Pt()) / myllmet.p4.Pt() : -10.;
            llmet.push_back(myllmet);
        }
    }


    // ***** 
    // Jets and dijets 
    // ***** 
    const JetsProducer& alljets = producers.get<JetsProducer>(m_jets_producer);
    jets.clear();
    nJetsL = 0;
    nJetsT = 0;
    nBJetsL = 0;
    nBJetsM = 0;
    nBJetsT = 0;
    for (unsigned int ijet = 0; ijet < alljets.p4.size(); ijet++)
    {
        float correctionFactor = m_applyBJetRegression ? alljets.regPt[ijet] / alljets.p4[ijet].Pt() : 1.;
        if (!((fabs(alljets.p4[ijet].Eta()) < m_jetEtaCut) && alljets.passLooseID[ijet]))
                continue;

        HH::Jet myjet;
        myjet.passPtCut =  (alljets.p4[ijet].Pt() * correctionFactor > m_jetPtCut);
        // Jet veto based on DR with all selected leptons
        float DRj_l1 = -1;
        float DRj_l2 = -1;
        if (ll.size() > 0) {
            DRj_l1 = ROOT::Math::VectorUtil::DeltaR(myjet.p4, leptons[ll[0].ilep1].p4);
            DRj_l2 = ROOT::Math::VectorUtil::DeltaR(myjet.p4, leptons[ll[0].ilep2].p4); 
        }
        myjet.minDR_jet_selLeptons = std::min(DRj_l1, DRj_l2);

        float DRjl_allLep = std::numeric_limits<float>::max();
        for (auto lepton : leptons) {
            float DRjl_allLep_temp = ROOT::Math::VectorUtil::DeltaR(myjet.p4, lepton.p4);
            if (DRjl_allLep_temp < DRjl_allLep)
                DRjl_allLep = DRjl_allLep_temp;
        }
        myjet.minDR_jet_allLeptons = DRjl_allLep;
        if (myjet.minDR_jet_allLeptons < m_minDR_l_j_Cut)
            continue;

        myjet.p4 = alljets.p4[ijet] * correctionFactor;
        myjet.idx = ijet;
        myjet.id_L = alljets.passLooseID[ijet];
        myjet.id_T = alljets.passTightID[ijet];
        myjet.id_TLV = alljets.passTightLeptonVetoID[ijet];
        myjet.CSV = alljets.getBTagDiscriminant(ijet, "pfCombinedInclusiveSecondaryVertexV2BJetTags");
        myjet.JP = alljets.getBTagDiscriminant(ijet, "pfJetProbabilityBJetTags");
        float mybtag = alljets.getBTagDiscriminant(ijet, m_jet_bDiscrName);
        myjet.btag_L = mybtag > m_jet_bDiscrCut_loose;
        myjet.btag_M = mybtag > m_jet_bDiscrCut_medium;
        myjet.btag_T = mybtag > m_jet_bDiscrCut_tight;
        // GEN MATCHING
        myjet.tt_matched = false;
        myjet.tt_gen_p4 = null_p4;
        myjet.tt_gen_DR = std::numeric_limits<float>::max();
        myjet.tt_parton_p4 = null_p4;
        myjet.tt_parton_DR = std::numeric_limits<float>::max();
        if (!event.isRealData() && gen_b_beforeFSR != 0 && gen_bbar_beforeFSR != 0 && gen_b !=0 && gen_bbar != 0) {
            float temp_DR_jet_t = std::numeric_limits<float>::max(); 
            float temp_DR_jet_antit = std::numeric_limits<float>::max(); 
            float DR_jetParton = std::numeric_limits<float>::max(); 
            uint16_t parton_id = 0;
            uint16_t parton_id_beforeFSR = 0;
            bool goodID = false;
            if (gp->pruned_pdg_id[gen_b] == 5) { // should be useless check
                temp_DR_jet_t = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[gen_b], myjet.p4);
                goodID = true;
            }
            if (gp->pruned_pdg_id[gen_bbar] == -5) { 
                temp_DR_jet_antit = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[gen_bbar], myjet.p4);
                goodID = true;
            }
            if (goodID) {
                parton_id = (temp_DR_jet_t < temp_DR_jet_antit) ? gen_b : gen_bbar;
                DR_jetParton = (temp_DR_jet_t < temp_DR_jet_antit) ? temp_DR_jet_t : temp_DR_jet_antit;
                parton_id_beforeFSR = (temp_DR_jet_t < temp_DR_jet_antit) ? gen_b_beforeFSR : gen_bbar_beforeFSR;
                myjet.tt_gen_DR = DR_jetParton;
                myjet.tt_gen_p4 = gp->pruned_p4[parton_id];
                myjet.tt_parton_p4 = gp->pruned_p4[parton_id_beforeFSR];
                myjet.tt_parton_DR = ROOT::Math::VectorUtil::DeltaR(gp->pruned_p4[parton_id_beforeFSR], myjet.p4);
                if (myjet.tt_gen_DR < tt_genJetMatching_maxDR)
                    myjet.tt_matched = true;
            }
        }
        myjet.gen_matched_bParton = (std::abs(alljets.partonFlavor[ijet]) == 5);
        myjet.gen_matched_bHadron = (std::abs(alljets.hadronFlavor[ijet]) == 5);
        myjet.gen_matched = alljets.matched[ijet];
        myjet.gen_p4 = myjet.gen_matched ? alljets.gen_p4[ijet] : null_p4;
        myjet.gen_DR = myjet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjet.p4, myjet.gen_p4) : -1.;
        myjet.gen_DPtOverPt = myjet.gen_matched ? (myjet.p4.Pt() - myjet.gen_p4.Pt()) / myjet.p4.Pt() : -10.;
        myjet.gen_b = (std::abs(alljets.hadronFlavor[ijet]) == 5); // redundant with gen_matched_bHadron defined above
        myjet.gen_c = (std::abs(alljets.hadronFlavor[ijet]) == 4);
        myjet.gen_l = (std::abs(alljets.hadronFlavor[ijet]) < 4);
        // counters
        if (myjet.id_L)
            nJetsL++;
        if (myjet.id_T)
            nJetsT++;
        if (myjet.btag_L)
            nBJetsL++;
        if (myjet.btag_M)
            nBJetsM++;
        if (myjet.btag_T)
            nBJetsT++;
        jets.push_back(myjet);
    }
    nJets = jets.size();

    jj.clear();

    for (unsigned int ijet1 = 0; ijet1 < jets.size(); ijet1++)
    {
        for (unsigned int ijet2 = ijet1 + 1; ijet2 < jets.size(); ijet2++)
        {
            HH::Dijet myjj;
            myjj.p4 = jets[ijet1].p4 + jets[ijet2].p4;
            myjj.idxs = std::make_pair(jets[ijet1].idx, jets[ijet2].idx);
            myjj.ijet1 = ijet1;
            myjj.ijet2 = ijet2;
            for (unsigned int ijet3 = 0; ijet3 < jets.size(); ijet3++) {
                if (ijet3 != ijet1 && ijet3 != ijet2) {
                    myjj.iExtraJets.push_back(ijet3);
                    myjj.min_DR_j_extraj.push_back(std::min(ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, jets[ijet3].p4), ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, jets[ijet3].p4)));
                    myjj.M_j_j_extraj.push_back((jets[ijet1].p4 + jets[ijet2].p4 + jets[ijet3].p4).M());
                }
            }
            myjj.btag_LL = jets[ijet1].btag_L && jets[ijet2].btag_L;
            myjj.btag_LM = (jets[ijet1].btag_L && jets[ijet2].btag_M) || (jets[ijet2].btag_L && jets[ijet1].btag_M);
            myjj.btag_LT = (jets[ijet1].btag_L && jets[ijet2].btag_T) || (jets[ijet2].btag_L && jets[ijet1].btag_T);
            myjj.btag_ML = (jets[ijet1].btag_M && jets[ijet2].btag_L) || (jets[ijet2].btag_M && jets[ijet1].btag_L);
            myjj.btag_MM = jets[ijet1].btag_M && jets[ijet2].btag_M;
            myjj.btag_MT = (jets[ijet1].btag_M && jets[ijet2].btag_T) || (jets[ijet2].btag_M && jets[ijet1].btag_T);
            myjj.btag_TL = (jets[ijet1].btag_T && jets[ijet2].btag_L) || (jets[ijet2].btag_T && jets[ijet1].btag_L);
            myjj.btag_TM = (jets[ijet1].btag_T && jets[ijet2].btag_M) || (jets[ijet2].btag_T && jets[ijet1].btag_M);
            myjj.btag_TT = jets[ijet1].btag_T && jets[ijet2].btag_T;
            myjj.sumCSV = jets[ijet1].CSV + jets[ijet2].CSV;
            myjj.sumJP = jets[ijet1].JP + jets[ijet2].JP;
            myjj.DR_j_j = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, jets[ijet2].p4);
            myjj.DPhi_j_j = fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[ijet1].p4, jets[ijet2].p4));
            myjj.ht_j_j = jets[ijet1].p4.Pt() + jets[ijet2].p4.Pt();
            myjj.gen_matched_bbPartons = jets[ijet1].gen_matched_bParton && jets[ijet2].gen_matched_bParton; 
            myjj.gen_matched_bbHadrons = jets[ijet1].gen_matched_bHadron && jets[ijet2].gen_matched_bHadron; 
            myjj.gen_matched = jets[ijet1].gen_matched && jets[ijet2].gen_matched;
            myjj.gen_p4 = myjj.gen_matched ? jets[ijet1].gen_p4 + jets[ijet2].gen_p4 : null_p4;
            myjj.gen_DR = myjj.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjj.p4, myjj.gen_p4) : -1.;
            myjj.gen_DPtOverPt = myjj.gen_matched ? (myjj.p4.Pt() - myjj.gen_p4.Pt()) / myjj.p4.Pt() : -10.;
            myjj.gen_bb = (jets[ijet1].gen_b && jets[ijet2].gen_b);
            myjj.gen_bc = (jets[ijet1].gen_b && jets[ijet2].gen_c) || (jets[ijet1].gen_c && jets[ijet2].gen_b);
            myjj.gen_bl = (jets[ijet1].gen_b && jets[ijet2].gen_l) || (jets[ijet1].gen_l && jets[ijet2].gen_b);
            myjj.gen_cc = (jets[ijet1].gen_c && jets[ijet2].gen_c);
            myjj.gen_cl = (jets[ijet1].gen_c && jets[ijet2].gen_l) || (jets[ijet1].gen_l && jets[ijet2].gen_c);
            myjj.gen_ll = (jets[ijet1].gen_l && jets[ijet2].gen_l);
            jj.push_back(myjj);
        }
    }
    // Sort Dijet such as the jet are pt ordered
    std::sort(jj.begin(), jj.end(), [](const HH::Dijet& jj1, const HH::Dijet& jj2) { return jj1.ht_j_j > jj2.ht_j_j; });     

    // ********** 
    // lljj, llbb, +pf_met
    // ********** 
    llmetjj.clear();
    for (unsigned int illmet = 0; illmet < llmet.size(); illmet++)
    {
        for (unsigned int ijj = 0; ijj < jj.size(); ijj++)
        {
            unsigned int imet = llmet[illmet].imet;
            unsigned int ill = llmet[illmet].ill;
            unsigned int ijet1 = jj[ijj].ijet1;
            unsigned int ijet2 = jj[ijj].ijet2;
            unsigned int ilep1 = ll[ill].ilep1;
            unsigned int ilep2 = ll[ill].ilep2;
            HH::DileptonMetDijet myllmetjj;
            myllmetjj.p4 = ll[ill].p4 + jj[ijj].p4 + met[imet].p4;
            myllmetjj.lep1_p4 = leptons[ilep1].p4;
            myllmetjj.lep2_p4 = leptons[ilep2].p4;
            myllmetjj.jet1_p4 = jets[ijet1].p4;
            myllmetjj.jet2_p4 = jets[ijet2].p4;
            myllmetjj.met_p4 = met[imet].p4;
            myllmetjj.ll_p4 = ll[ill].p4;
            myllmetjj.jj_p4 = jj[ijj].p4;
            myllmetjj.ht_j_j = jj[ijj].ht_j_j;
            myllmetjj.ht_l_l = ll[ill].ht_l_l;
            myllmetjj.lljj_p4 = ll[ill].p4 + jj[ijj].p4;
            // gen info
            myllmetjj.gen_matched = ll[ill].gen_matched && jj[ijj].gen_matched && met[imet].gen_matched;
            myllmetjj.gen_p4 = myllmetjj.gen_matched ? ll[ill].gen_p4 + jj[ijj].gen_p4 + met[imet].gen_p4 : null_p4;
            myllmetjj.gen_DR = myllmetjj.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myllmetjj.p4, myllmetjj.gen_p4) : -1.;
            myllmetjj.gen_DPhi = myllmetjj.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(myllmetjj.p4, myllmetjj.gen_p4)) : -1.;
            myllmetjj.gen_DPtOverPt = myllmetjj.gen_matched ? (myllmetjj.p4.Pt() - myllmetjj.gen_p4.Pt()) / myllmetjj.p4.Pt() : -10.;
            myllmetjj.gen_lep1_p4 = leptons[ilep1].gen_p4;
            myllmetjj.gen_lep2_p4 = leptons[ilep2].gen_p4;
            myllmetjj.gen_jet1_p4 = leptons[ijet1].gen_p4;
            myllmetjj.gen_jet2_p4 = leptons[ijet2].gen_p4;
            myllmetjj.gen_met_p4 = met[imet].gen_p4;
            myllmetjj.gen_ll_p4 = ll[ill].gen_p4;
            myllmetjj.gen_jj_p4 = jj[ijj].gen_p4;
            myllmetjj.gen_lljj_p4 = ll[ill].gen_p4 + jj[ijj].gen_p4;
            // blind copy of the jj content
            myllmetjj.ijet1 = jj[ijj].ijet1;
            myllmetjj.ijet2 = jj[ijj].ijet2;
            myllmetjj.btag_LL = jj[ijj].btag_LL;
            myllmetjj.btag_LM = jj[ijj].btag_LM;
            myllmetjj.btag_LT = jj[ijj].btag_LT;
            myllmetjj.btag_ML = jj[ijj].btag_ML;
            myllmetjj.btag_MM = jj[ijj].btag_MM;
            myllmetjj.btag_MT = jj[ijj].btag_MT;
            myllmetjj.btag_TL = jj[ijj].btag_TL;
            myllmetjj.btag_TM = jj[ijj].btag_TM;
            myllmetjj.btag_TT = jj[ijj].btag_TT;
            myllmetjj.sumCSV = jj[ijj].sumCSV;
            myllmetjj.sumJP = jj[ijj].sumJP;
            myllmetjj.DR_j_j = jj[ijj].DR_j_j;
            myllmetjj.DPhi_j_j = jj[ijj].DPhi_j_j;
            myllmetjj.gen_matched_bbPartons = jj[ijj].gen_matched_bbPartons;
            myllmetjj.gen_matched_bbHadrons = jj[ijj].gen_matched_bbHadrons;
            myllmetjj.gen_bb = jj[ijj].gen_bb;
            myllmetjj.gen_bc = jj[ijj].gen_bc;
            myllmetjj.gen_bl = jj[ijj].gen_bl;
            myllmetjj.gen_cc = jj[ijj].gen_cc;
            myllmetjj.gen_cl = jj[ijj].gen_cl;
            myllmetjj.gen_ll = jj[ijj].gen_ll;
            myllmetjj.iExtraJets = jj[ijj].iExtraJets;
            myllmetjj.min_DR_j_extraj = jj[ijj].min_DR_j_extraj;
            myllmetjj.M_j_j_extraj = jj[ijj].M_j_j_extraj;
            // blind copy of the llmet content
            myllmetjj.ilep1 = ll[ill].ilep1;
            myllmetjj.ilep2 = ll[ill].ilep2;
            myllmetjj.isOS = ll[ill].isOS;
            myllmetjj.isMuMu = ll[ill].isMuMu;
            myllmetjj.isElEl = ll[ill].isElEl;
            myllmetjj.isElMu = ll[ill].isElMu;
            myllmetjj.isMuEl = ll[ill].isMuEl;
            myllmetjj.isSF = ll[ill].isSF;
            myllmetjj.id_LL = ll[ill].id_LL;
            myllmetjj.id_LM = ll[ill].id_LM;
            myllmetjj.id_LT = ll[ill].id_LT;
            myllmetjj.id_LHWW = ll[ill].id_LHWW;
            myllmetjj.id_ML = ll[ill].id_ML;
            myllmetjj.id_MM = ll[ill].id_MM;
            myllmetjj.id_MT = ll[ill].id_MT;
            myllmetjj.id_MHWW = ll[ill].id_MHWW;
            myllmetjj.id_TL = ll[ill].id_TL;
            myllmetjj.id_TM = ll[ill].id_TM;
            myllmetjj.id_TT = ll[ill].id_TT;
            myllmetjj.id_THWW = ll[ill].id_THWW;
            myllmetjj.id_HWWL = ll[ill].id_HWWL;
            myllmetjj.id_HWWM = ll[ill].id_HWWM;
            myllmetjj.id_HWWT = ll[ill].id_HWWT;
            myllmetjj.id_HWWHWW = ll[ill].id_HWWHWW;
            myllmetjj.iso_LL = ll[ill].iso_LL;
            myllmetjj.iso_LT = ll[ill].iso_LT;
            myllmetjj.iso_LHWW = ll[ill].iso_LHWW;
            myllmetjj.iso_TL = ll[ill].iso_TL;
            myllmetjj.iso_TT = ll[ill].iso_TT;
            myllmetjj.iso_THWW = ll[ill].iso_THWW;
            myllmetjj.iso_HWWL = ll[ill].iso_HWWL;
            myllmetjj.iso_HWWT = ll[ill].iso_HWWT;
            myllmetjj.iso_HWWHWW = ll[ill].iso_HWWHWW;
            myllmetjj.DR_l_l = ll[ill].DR_l_l;
            myllmetjj.DPhi_l_l = ll[ill].DPhi_l_l;
            myllmetjj.trigger_efficiency = ll[ill].trigger_efficiency;
            myllmetjj.trigger_efficiency_downVariated = ll[ill].trigger_efficiency_downVariated;
            myllmetjj.trigger_efficiency_downVariated_Arun = ll[ill].trigger_efficiency_downVariated_Arun;
            myllmetjj.trigger_efficiency_upVariated = ll[ill].trigger_efficiency_upVariated;
            myllmetjj.trigger_efficiency_upVariated_Arun = ll[ill].trigger_efficiency_upVariated_Arun;
            myllmetjj.ill = ill;
            myllmetjj.imet = imet;
            myllmetjj.isNoHF = met[imet].isNoHF;
            myllmetjj.DPhi_ll_met = llmet[illmet].DPhi_ll_met;
            myllmetjj.minDPhi_l_met = llmet[illmet].minDPhi_l_met; 
            myllmetjj.maxDPhi_l_met = llmet[illmet].maxDPhi_l_met;
            myllmetjj.MT = llmet[illmet].MT;
            myllmetjj.MT_formula = llmet[illmet].MT_formula;
            myllmetjj.projectedMet = llmet[illmet].projectedMet;
            // content specific to HH::DijetMet
            // NB: computed for the first time here, no intermediate jjmet collection
            myllmetjj.DPhi_jj_met = fabs(ROOT::Math::VectorUtil::DeltaPhi(jj[ijj].p4, met[imet].p4));
            myllmetjj.minDPhi_j_met = std::min(fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[jj[ijj].ijet1].p4, met[imet].p4)), fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[jj[ijj].ijet2].p4, met[imet].p4)));
            myllmetjj.maxDPhi_j_met = std::max(fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[jj[ijj].ijet1].p4, met[imet].p4)), fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[jj[ijj].ijet2].p4, met[imet].p4)));
            // content specific to HH::DileptonMetDijet
            myllmetjj.illmet = illmet;
            myllmetjj.ijj = ijj;
            float DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2;
            DR_j1l1 = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, leptons[ilep1].p4);
            DR_j1l2 = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, leptons[ilep2].p4);
            DR_j2l1 = ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, leptons[ilep1].p4);
            DR_j2l2 = ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, leptons[ilep2].p4);
            myllmetjj.maxDR_l_j = std::max({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            myllmetjj.minDR_l_j = std::min({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            myllmetjj.DR_ll_jj = ROOT::Math::VectorUtil::DeltaR(ll[ill].p4, jj[ijj].p4);
            myllmetjj.DPhi_ll_jj = fabs(ROOT::Math::VectorUtil::DeltaPhi(ll[ill].p4, jj[ijj].p4));
            myllmetjj.DR_llmet_jj = ROOT::Math::VectorUtil::DeltaR(llmet[illmet].p4, jj[ijj].p4);
            myllmetjj.DPhi_llmet_jj = fabs(ROOT::Math::VectorUtil::DeltaPhi(llmet[illmet].p4, jj[ijj].p4));
            myllmetjj.cosThetaStar_CS = fabs(getCosThetaStar_CS(llmet[illmet].p4, jj[ijj].p4));
            myllmetjj.MT_fullsystem = myllmetjj.p4.Mt();
            // Some selection
            if (myllmetjj.minDR_l_j < m_minDR_l_j_Cut)
                continue;
            // Counters
            tmp_count_has2leptons_1llmetjj = event_weight;
            if (myllmetjj.isElEl)
                tmp_count_has2leptons_elel_1llmetjj = event_weight;
            if (myllmetjj.isElMu)
                tmp_count_has2leptons_elmu_1llmetjj = event_weight;
            if (myllmetjj.isElEl)
                tmp_count_has2leptons_muel_1llmetjj = event_weight;
            if (myllmetjj.isElEl)
                tmp_count_has2leptons_mumu_1llmetjj = event_weight;
            if (myllmetjj.btag_MM)
            {
                tmp_count_has2leptons_1llmetjj_2btagM = event_weight;
                if (myllmetjj.isElEl)
                    tmp_count_has2leptons_elel_1llmetjj_2btagM = event_weight;
                if (myllmetjj.isElMu)
                    tmp_count_has2leptons_elmu_1llmetjj_2btagM = event_weight;
                if (myllmetjj.isElEl)
                    tmp_count_has2leptons_muel_1llmetjj_2btagM = event_weight;
                if (myllmetjj.isElEl)
                    tmp_count_has2leptons_mumu_1llmetjj_2btagM = event_weight;
            }
            // Fill
            llmetjj.push_back(myllmetjj);
        }
    }
    // Different llmetjj candidate
    // PT ORDERED asking the first two jets to be b-tagged
    llmetjj_HWWleptons_nobtag_pt.clear();
    llmetjj_HWWleptons_btagL_pt.clear();
    llmetjj_HWWleptons_btagM_pt.clear();
    llmetjj_HWWleptons_btagMT_pt.clear();
    llmetjj_HWWleptons_btagML_pt.clear();
    llmetjj_HWWleptons_btagT_pt.clear();

    if (llmetjj.size()>0) {
        auto llmetjj_cand = llmetjj.at(0);
        llmetjj_HWWleptons_nobtag_pt.push_back(llmetjj_cand);
        if (llmetjj_cand.btag_LL)
            llmetjj_HWWleptons_btagL_pt.push_back(llmetjj_cand);
        if (llmetjj_cand.btag_MM)
            llmetjj_HWWleptons_btagM_pt.push_back(llmetjj_cand);
        if (llmetjj_cand.btag_ML || llmetjj_cand.btag_LM)
            llmetjj_HWWleptons_btagML_pt.push_back(llmetjj_cand);
        if (llmetjj_cand.btag_MT || llmetjj_cand.btag_TM)
            llmetjj_HWWleptons_btagMT_pt.push_back(llmetjj_cand);
        if (llmetjj_cand.btag_TT)
            llmetjj_HWWleptons_btagT_pt.push_back(llmetjj_cand);
        // PT ORDERED asking at least two jets to be b-tagged
        for (auto llmetjj_cand_incl : llmetjj){
            if (llmetjj_cand_incl.btag_MM && llmetjj_HWWleptons_btagM_pt_inclusive.size() == 0)
                llmetjj_HWWleptons_btagM_pt_inclusive.push_back(llmetjj_cand_incl);
        }

        // CSV ORDERED
        std::sort(llmetjj.begin(), llmetjj.end(), [](const HH::DileptonMetDijet& llmetjj1, const HH::DileptonMetDijet& llmetjj2) { return llmetjj1.sumCSV > llmetjj2.sumCSV; });     
        llmetjj_HWWleptons_nobtag_csv.clear();
        llmetjj_HWWleptons_btagL_csv.clear();
        llmetjj_HWWleptons_btagM_csv.clear();
        llmetjj_HWWleptons_btagMT_csv.clear();
        llmetjj_HWWleptons_btagML_csv.clear();
        llmetjj_HWWleptons_btagT_csv.clear();

        llmetjj_HWWleptons_nobtag_csv.push_back(llmetjj.at(0));
        for (auto llmetjj_cand : llmetjj){
            if (llmetjj_cand.btag_LL && llmetjj_HWWleptons_btagL_csv.size()==0)
                llmetjj_HWWleptons_btagL_csv.push_back(llmetjj_cand);
            if (llmetjj_cand.btag_MM && llmetjj_HWWleptons_btagM_csv.size()==0)
                llmetjj_HWWleptons_btagM_csv.push_back(llmetjj_cand);
            if ((llmetjj_cand.btag_ML || llmetjj_cand.btag_LM) && llmetjj_HWWleptons_btagML_pt.size()==0)
                llmetjj_HWWleptons_btagML_csv.push_back(llmetjj_cand);
            if ((llmetjj_cand.btag_MT || llmetjj_cand.btag_TM) && llmetjj_HWWleptons_btagMT_pt.size()==0)
                llmetjj_HWWleptons_btagMT_csv.push_back(llmetjj_cand);
            if (llmetjj_cand.btag_TT && llmetjj_HWWleptons_btagT_pt.size()==0)
                llmetjj_HWWleptons_btagT_csv.push_back(llmetjj_cand);
        }
    }

    count_has2leptons += tmp_count_has2leptons;
    count_has2leptons_elel += tmp_count_has2leptons_elel;
    count_has2leptons_elmu += tmp_count_has2leptons_elmu;
    count_has2leptons_muel += tmp_count_has2leptons_muel;
    count_has2leptons_mumu += tmp_count_has2leptons_mumu;
    count_has2leptons_1llmetjj += tmp_count_has2leptons_1llmetjj;
    count_has2leptons_elel_1llmetjj += tmp_count_has2leptons_elel_1llmetjj;
    count_has2leptons_elmu_1llmetjj += tmp_count_has2leptons_elmu_1llmetjj;
    count_has2leptons_muel_1llmetjj += tmp_count_has2leptons_muel_1llmetjj;
    count_has2leptons_mumu_1llmetjj += tmp_count_has2leptons_mumu_1llmetjj;
    count_has2leptons_1llmetjj_2btagM += tmp_count_has2leptons_1llmetjj_2btagM;
    count_has2leptons_elel_1llmetjj_2btagM += tmp_count_has2leptons_elel_1llmetjj_2btagM;
    count_has2leptons_elmu_1llmetjj_2btagM += tmp_count_has2leptons_elmu_1llmetjj_2btagM;
    count_has2leptons_muel_1llmetjj_2btagM += tmp_count_has2leptons_muel_1llmetjj_2btagM;
    count_has2leptons_mumu_1llmetjj_2btagM += tmp_count_has2leptons_mumu_1llmetjj_2btagM;
}

float HHAnalyzer::getCosThetaStar_CS(const LorentzVector & h1, const LorentzVector & h2, float ebeam /*= 6500*/)
{// cos theta star angle in the Collins Soper frame
    LorentzVector p1, p2;
    p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
    p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

    LorentzVector hh = h1 + h2;
    ROOT::Math::Boost boost(-hh.X() / hh.T(), -hh.Y() / hh.T(), -hh.Z() / hh.T());
    p1 = boost(p1);
    p2 = boost(p2);
    LorentzVector newh1 = boost(h1);
    ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> CSaxis(p1.Vect().Unit() - p2.Vect().Unit());

    return cos(ROOT::Math::VectorUtil::Angle(CSaxis.Unit(), newh1.Vect().Unit()));
}

void HHAnalyzer::fillTriggerEfficiencies(const Lepton & lep1, const Lepton & lep2, Dilepton & dilep) {

    float eff_lep1_leg1 = 1.;
    float eff_lep1_leg2 = 1.;
    float eff_lep1_tkleg2 = 0.;
    float eff_lep2_leg1 = 1.;
    float eff_lep2_leg2 = 1.;
    float eff_lep2_tkleg2 = 0.;

    if (lep1.isMu && lep2.isMu) {
        eff_lep1_leg1 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        eff_lep1_leg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8orIsoTkMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        //eff_lep1_tkleg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_TkMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        eff_lep2_leg1 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
        eff_lep2_leg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8orIsoTkMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
        //eff_lep2_tkleg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_TkMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
    }
    else if (lep1.isMu && lep2.isEl) {
        eff_lep1_leg1 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        eff_lep1_leg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        eff_lep2_leg1 = m_hlt_efficiencies["Ele17_12Leg1"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
        eff_lep2_leg2 = m_hlt_efficiencies["Ele17_12Leg2"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
    }
    else if (lep1.isEl && lep2.isMu) {
        eff_lep1_leg1 = m_hlt_efficiencies["Ele17_12Leg1"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        eff_lep1_leg2 = m_hlt_efficiencies["Ele17_12Leg2"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        eff_lep2_leg1 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
        eff_lep2_leg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
    }
    else if (lep1.isEl && lep2.isEl){
        eff_lep1_leg1 = m_hlt_efficiencies["Ele17_12Leg1"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        eff_lep1_leg2 = m_hlt_efficiencies["Ele17_12Leg2"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
        eff_lep2_leg1 = m_hlt_efficiencies["Ele17_12Leg1"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
        eff_lep2_leg2 = m_hlt_efficiencies["Ele17_12Leg2"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
    }
    else 
        std::cout << "We have something else then el or mu !!" << std::endl;

    float error_eff_lep1_leg1_up = 0.;
    float error_eff_lep1_leg2_up = 0.;
    float error_eff_lep1_tkleg2_up = 0.;
    float error_eff_lep2_leg1_up = 0.;
    float error_eff_lep2_leg2_up = 0.;
    float error_eff_lep2_tkleg2_up = 0.;

    if (lep1.isMu && lep2.isMu) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8orIsoTkMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        //error_eff_lep1_tkleg2_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_TkMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8orIsoTkMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
        //error_eff_lep2_tkleg2_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_TkMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
    }
    else if (lep1.isMu && lep2.isEl) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies["Ele17_12Leg1"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies["Ele17_12Leg2"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
    }
    else if (lep1.isEl && lep2.isMu) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies["Ele17_12Leg1"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies["Ele17_12Leg2"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
    }
    else if (lep1.isEl && lep2.isEl){
        error_eff_lep1_leg1_up = m_hlt_efficiencies["Ele17_12Leg1"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies["Ele17_12Leg2"].get({lep1.p4.Eta(), lep1.p4.Pt()})[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies["Ele17_12Leg1"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies["Ele17_12Leg2"].get({lep2.p4.Eta(), lep2.p4.Pt()})[2];
    }

    float error_eff_lep1_leg1_down = 0.;
    float error_eff_lep1_leg2_down = 0.;
    float error_eff_lep1_tkleg2_down = 0.;
    float error_eff_lep2_leg1_down = 0.;
    float error_eff_lep2_leg2_down = 0.;
    float error_eff_lep2_tkleg2_down = 0.;

    if (lep1.isMu && lep2.isMu) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8orIsoTkMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        //error_eff_lep1_tkleg2_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_TkMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8orIsoTkMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
        //error_eff_lep2_tkleg2_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_TkMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
    }
    else if (lep1.isMu && lep2.isEl) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies["Ele17_12Leg1"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies["Ele17_12Leg2"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
    }
    else if (lep1.isEl && lep2.isMu) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies["Ele17_12Leg1"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies["Ele17_12Leg2"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
    }
    else if (lep1.isEl && lep2.isEl){
        error_eff_lep1_leg1_down = m_hlt_efficiencies["Ele17_12Leg1"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies["Ele17_12Leg2"].get({lep1.p4.Eta(), lep1.p4.Pt()})[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies["Ele17_12Leg1"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies["Ele17_12Leg2"].get({lep2.p4.Eta(), lep2.p4.Pt()})[1];
    }


    float nominal = -(eff_lep1_leg1 * eff_lep2_leg1) +
        (1 - (1 - eff_lep1_leg2) * (1 - eff_lep1_tkleg2)) * eff_lep2_leg1 +
        eff_lep1_leg1 * (1 - (1 - eff_lep2_leg2) * (1 - eff_lep2_tkleg2));

    float error_squared_up =
        std::pow(1 - eff_lep2_leg1 - (1 - eff_lep2_leg2) * (1 - eff_lep2_tkleg2), 2) *
        std::pow(error_eff_lep1_leg1_up, 2) +
        std::pow(1 - eff_lep1_tkleg2, 2) * std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_leg2_up, 2) +
        std::pow(1 - eff_lep1_leg2, 2) * std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_tkleg2_up, 2) +
        std::pow(1 - eff_lep1_leg1 - (1 - eff_lep1_leg2) * (1 - eff_lep1_tkleg2), 2) *
        std::pow(error_eff_lep2_leg1_up, 2) +
        std::pow(eff_lep1_leg1, 2) * std::pow(1 - eff_lep2_tkleg2, 2) *
        std::pow(error_eff_lep2_leg2_up, 2) +
        std::pow(eff_lep1_leg1, 2) * std::pow(1 - eff_lep2_leg2, 2) *
        std::pow(error_eff_lep2_tkleg2_up, 2);

    float error_squared_down = 
        std::pow(1 - eff_lep2_leg1 - (1 - eff_lep2_leg2) * (1 - eff_lep2_tkleg2), 2) *
        std::pow(error_eff_lep1_leg1_down, 2) +
        std::pow(1 - eff_lep1_tkleg2, 2) * std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_leg2_down, 2) +
        std::pow(1 - eff_lep1_leg2, 2) * std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_tkleg2_down, 2) +
        std::pow(1 - eff_lep1_leg1 - (1 - eff_lep1_leg2) * (1 - eff_lep1_tkleg2), 2) *
        std::pow(error_eff_lep2_leg1_down, 2) +
        std::pow(eff_lep1_leg1, 2) * std::pow(1 - eff_lep2_tkleg2, 2) *
        std::pow(error_eff_lep2_leg2_down, 2) +
        std::pow(eff_lep1_leg1, 2) * std::pow(1 - eff_lep2_leg2, 2) *
        std::pow(error_eff_lep2_tkleg2_down, 2);

    dilep.trigger_efficiency = nominal;
    dilep.trigger_efficiency_upVariated = ((nominal + std::sqrt(error_squared_up)) > 1.)? 1. : nominal + std::sqrt(error_squared_up);
    dilep.trigger_efficiency_downVariated = ((nominal - std::sqrt(error_squared_down)) < 0.)? 0. : nominal - std::sqrt(error_squared_down);
    
    // Arun's method (not using the proper derivative formula)
    float X = eff_lep1_leg1 * eff_lep2_leg2 * std::sqrt((std::pow((error_eff_lep1_leg1_up/eff_lep1_leg1),2) + std::pow((error_eff_lep2_leg2_up/eff_lep2_leg2),2) ));
    float Y = eff_lep2_leg1 * eff_lep1_leg2 * std::sqrt((std::pow((error_eff_lep2_leg1_up/eff_lep2_leg1),2) + std::pow((error_eff_lep1_leg2_up/eff_lep1_leg2),2) ));
    float Z = eff_lep1_leg1 * eff_lep2_leg1 * std::sqrt((std::pow((error_eff_lep1_leg1_up/eff_lep1_leg1),2) + std::pow((error_eff_lep2_leg1_up/eff_lep2_leg1),2) ));
    float error_squared_up_Arun = X*X + Y*Y + Z*Z ;
    dilep.trigger_efficiency_upVariated_Arun = ((nominal + std::sqrt(error_squared_up_Arun)) > 1.)? 1. : nominal + std::sqrt(error_squared_up_Arun);

    X = eff_lep1_leg1 * eff_lep2_leg2 * std::sqrt((std::pow((error_eff_lep1_leg1_down/eff_lep1_leg1),2) + std::pow((error_eff_lep2_leg2_down/eff_lep2_leg2),2) ));
    Y = eff_lep2_leg1 * eff_lep1_leg2 * std::sqrt((std::pow((error_eff_lep2_leg1_down/eff_lep2_leg1),2) + std::pow((error_eff_lep1_leg2_down/eff_lep1_leg2),2) ));
    Z = eff_lep1_leg1 * eff_lep2_leg1 * std::sqrt((std::pow((error_eff_lep1_leg1_down/eff_lep1_leg1),2) + std::pow((error_eff_lep2_leg1_down/eff_lep2_leg1),2) ));
    float error_squared_down_Arun = X*X + Y*Y + Z*Z ;
    dilep.trigger_efficiency_downVariated_Arun = ((nominal - std::sqrt(error_squared_down_Arun)) < 0.)? 0. : nominal - std::sqrt(error_squared_down_Arun);

}

void HHAnalyzer::endJob(MetadataManager& metadata)
{
    metadata.add(this->m_name + "count_has2leptons", count_has2leptons);
    metadata.add(this->m_name + "count_has2leptons_elel", count_has2leptons_elel);
    metadata.add(this->m_name + "count_has2leptons_elmu", count_has2leptons_elmu);
    metadata.add(this->m_name + "count_has2leptons_muel", count_has2leptons_muel);
    metadata.add(this->m_name + "count_has2leptons_mumu", count_has2leptons_mumu);
    metadata.add(this->m_name + "count_has2leptons_1llmetjj", count_has2leptons_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_elel_1llmetjj", count_has2leptons_elel_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_elmu_1llmetjj", count_has2leptons_elmu_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_muel_1llmetjj", count_has2leptons_muel_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_mumu_1llmetjj", count_has2leptons_mumu_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_1llmetjj_2btagM", count_has2leptons_1llmetjj_2btagM);
    metadata.add(this->m_name + "count_has2leptons_elel_1llmetjj_2btagM", count_has2leptons_elel_1llmetjj_2btagM);
    metadata.add(this->m_name + "count_has2leptons_elmu_1llmetjj_2btagM", count_has2leptons_elmu_1llmetjj_2btagM);
    metadata.add(this->m_name + "count_has2leptons_muel_1llmetjj_2btagM", count_has2leptons_muel_1llmetjj_2btagM);
    metadata.add(this->m_name + "count_has2leptons_mumu_1llmetjj_2btagM", count_has2leptons_mumu_1llmetjj_2btagM);
}
