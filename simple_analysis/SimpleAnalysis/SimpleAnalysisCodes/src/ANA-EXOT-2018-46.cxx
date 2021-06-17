#include "SimpleAnalysisFramework/AnalysisClass.h"
#include <cmath>

DefineAnalysis(METHbb2018)

  void METHbb2018::Init() {
    // Signal region definitions
    // signal regions with 2 b-jets
    addRegions({"MET150200_2b", "MET200350_2b", "MET350500_2b"});
    addRegions({"MET500750_2b", "MET750_2b"});

    // signal regions with 3 or more b-jets
    addRegions({"MET150200_3b", "MET200350_3b", "MET350500_3b"});
    addRegions({"MET500_3b"});

    // Cutflow histograms
    // reco level cutflow merged (22 cuts)
    // - Initial
    // - PassGRL
    // - passLArTile
    // - Trigger
    // - HasVtx
    // - BadJet
    // - CosmicMuon
    // - BadMuon
    // - PFlow Electron veto
    // - IsMETTrigPassed
    // - 0 baseline electrons
    // - 0 baseline muons
    // - Tau Veto
    // - Extended Tau Veto
    // - MetTST>500
    // - >=1 fat-jets
    // - >= 2 b-tagged track jets
    // - (mJ > 40 || mJ_corr > 40)
    // - (mJ> 50 && mJ < 270)
    // - N_associated_trkJets>=2
    // - (TrackJet_1passOR = true && TrackJet_2passOR = true)
    // - |DeltaPhiMin3|>20deg
    addHistogram("Cutflow_merged", 22, -0.5, 21.5);
    // reco level cutflow resolved (26 cuts)
    // - Initial
    // - PassGRL
    // - passLArTile
    // - Trigger
    // - HasVtx
    // - BadJet
    // - CosmicMuon
    // - BadMuon
    // - PFlow Electron veto
    // - IsMETTrigPassed
    // - 0 baseline electrons
    // - 0 baseline muons
    // - Tau Veto
    // - Extended Tau Veto
    // - MetTST>150
    // - MetTST<=500
    // - >=2 jets
    // - >=2 b-tags
    // - (mjj > 40 || mjj_corr > 40)
    // - (mjj > 50 && mjj < 280)
    // - |DeltaPhiMin3|>20deg
    // - METSig>12
    // - (pt_jj > 100000 || pt_jj_corr > 100000)
    // - mT_METclosestBJet > 170000
    // - mT_METfurthestBJet > 200000
    // - <=4 jets
    addHistogram("Cutflow_resolved", 26, -0.5, 25.5);

    // Inclusive histograms
    addHistogram("MET", 40, 0., 2000.);
    addHistogram("mHiggs_resolved", 22, 50., 280.);
    addHistogram("mHiggs_merged", 11, 50., 270.);
    addHistogram("N_Bjets_resolved", 6, -0.5, 5.5);
    addHistogram("N_Bjets_merged_inside", 6, -0.5, 5.5);
    addHistogram("N_Bjets_merged_outside", 6, -0.5, 5.5);

    // Signal region histograms
    std::vector<float> binning_MET750 = {50., 90., 150., 270.};
    std::vector<float> binning_MET350500_3b = {50., 110., 150., 280.};
    addHistogram("MET150200_2b", 46, 50., 280.);
    addHistogram("MET200350_2b", 46, 50., 280.);
    addHistogram("MET350500_2b", 23, 50., 280.);
    addHistogram("MET500750_2b", 11, 50., 270.);
    addHistogram("MET750_2b", binning_MET750);
    addHistogram("MET150200_3b", 23, 50., 280.);
    addHistogram("MET200350_3b", 23, 50., 280.);
    addHistogram("MET350500_3b", binning_MET350500_3b);
    addHistogram("MET500_3b", binning_MET750);
}

void METHbb2018::ProcessEvent(AnalysisEvent *event) {
  // Object definitions
  // base leptons (used in veto definitions)
  auto baselineElectrons = event->getElectrons(7, 2.47, ELooseBLLH | EIsoFixedCutLoose);
  auto baselineMuons = event->getMuons(7, 2.5, MuLoose | MuIsoFixedCutLoose);
  auto baselineTaus = event->getTaus(20, 2.5, TauRNNVeryLoose);
  // central jets: pt > 20 GeV, |eta| < 2.5
  auto centralJets = event->getJets(20., 2.5, PFlowJet | JVTMedium);
  // forward jets: pt > 30 GeV, 2.5 < |eta| < 4.5 -> need to get creative with "filterCrack"
  auto forwardJets = event->getJets(30., 4.5);
  forwardJets = filterCrack(forwardJets, 2.5, 4.5);
 	sortObjectsByPt(centralJets);
	sortObjectsByPt(forwardJets);
  auto allJets = centralJets + forwardJets;
  auto fatJets = event->getFatJets(200., 2.0);
  auto trackJets = event->getTrackJets(10., 2.5);

  // missing transverse momentum and missing transverse momentum significance
  auto metVec = event->getMET();
  double metSignificance = event->getMETSignificance();
  double met = metVec.Et();

  // Overlap removal - including with object Pt-dependent radius calculation
  auto radiusCalcJet = [](const AnalysisObject &, const AnalysisObject &muon) {
    return std::min(0.4, 0.04 + 10 / muon.Pt());
  };
  auto radiusCalcMuon = [](const AnalysisObject &muon, const AnalysisObject &) {
    return std::min(0.4, 0.04 + 10 / muon.Pt());
  };
  baselineTaus = overlapRemoval(baselineTaus, baselineElectrons, 0.2);
  baselineElectrons = overlapRemoval(baselineElectrons, baselineMuons, 0.01);
  centralJets = overlapRemoval(centralJets, baselineElectrons, 0.2, NOT(BTag77DL1));
  baselineElectrons = overlapRemoval(baselineElectrons, centralJets, 0.4);
  centralJets = overlapRemoval(centralJets, baselineMuons, radiusCalcJet, LessThan3Tracks);
  baselineMuons = overlapRemoval(baselineMuons, centralJets, radiusCalcMuon);
  centralJets = overlapRemoval(centralJets, baselineTaus, 0.2);
  fatJets = overlapRemoval(fatJets, baselineElectrons, 1.0);

  // Advanced object definitions
  auto bjets = filterObjects(centralJets, 20., 2.5, BTag77DL1);
  // must resort to truth b-jet information for merged
  auto btrackjets = filterObjects(trackJets, 10., 2.5, BTag77DL1);
  auto baselineLeptons = baselineElectrons + baselineMuons;
  sortObjectsByPt(bjets);
  sortObjectsByPt(btrackjets);
	sortObjectsByPt(centralJets);
	sortObjectsByPt(fatJets);
  sortObjectsByPt(baselineLeptons);

  // Event-level observable reconstruction
  int nBjets = bjets.size();
  float dphiMin3 = minDphi(metVec, allJets, 3);
  float mt_min = 0;
  float mt_max = 0;
  float mindR = 999.;
  float maxdR = -1;
  
  for (const auto &jet : bjets) {
    if (metVec.DeltaR(jet) < mindR) {
      mindR = metVec.DeltaR(jet);
      mt_min = calcMT(jet, metVec);
    }
    if (metVec.DeltaR(jet) > maxdR) {
      maxdR = metVec.DeltaR(jet);
      mt_max = calcMT(jet, metVec);
    }
  }
  
  // Higgs candidate reconstruction
  float mHiggs = 0.;
  float ptHiggs = 0.;
  int nJetsMerged_inside = 0;
  int nBJetsMerged_inside = 0;
  int nBJetsMerged_outside = 0;
  if (met > 500) {
    if (fatJets.size() > 0) {
      mHiggs = fatJets[0].M();
      ptHiggs = fatJets[0].Pt();
      for (const auto &jet : trackJets) {
        if (jet.DeltaR(fatJets[0]) < 1.) {
          nJetsMerged_inside++;
          // only consider leading two associated track jets for b-tagging
          if (jet.pass(BTag77DL1) && nJetsMerged_inside < 3) {
            nBJetsMerged_inside++;
          }
        } else {
          if (jet.pass(BTag77DL1)) {
            nBJetsMerged_outside++;
          }
        }
      }
    }
  } else {
    if (nBjets >=2) {
      mHiggs = (bjets[0]+bjets[1]).M();
      ptHiggs = (bjets[0]+bjets[1]).Pt();
    }
  }

  // Event selection
  bool passMerged = false;
  bool passResolved = false;

  // common cuts in merged and resolved event selections
  bool evtsel_electronVeto = (baselineElectrons.size() == 0);
  bool evtsel_muonVeto = (baselineMuons.size() == 0);
  bool evtsel_tauVeto = (baselineTaus.size() == 0);
  bool evtsel_extendedTauVeto = true; // not implemented on truth level
  bool evtsel_minDPhi20 = (dphiMin3 > 20. * M_PI / 180.);

  // cuts in merged selection
  bool evtsel_met500 = (met > 500.);
  bool evtsel_nfatjets = (fatJets.size() >= 1);
  bool evtsel_nbtrackjets = (nBJetsMerged_inside >= 2);
  bool evtsel_mass40_merged = (mHiggs > 40.);
  bool evtsel_massRange_merged = (mHiggs > 50. && mHiggs < 270.);
  bool evtsel_ntrackjets_associated = (nJetsMerged_inside >= 2);

  // cuts in resolved selection
  bool evtsel_met150 = (met > 150.);
  bool evtsel_metleq500 = (met <= 500.);
  bool evtsel_njets = (centralJets.size() >=2);
  bool evtsel_nbjets = (nBjets >=2);
  bool evtsel_mass40_resolved = (mHiggs > 40.);
  bool evtsel_massRange_resolved = (mHiggs > 50. && mHiggs < 280.);
  bool evtsel_metSig = (metSignificance > 12.);

  bool evtsel_mt_mindR = (mt_min > 170.);
  bool evtsel_mt_maxdR = (mt_max > 200.);
  bool evtsel_njets_max = ((nBjets == 2 && centralJets.size() <= 4) || (nBjets >=3 && centralJets.size() <= 5));
  bool evtsel_ptHiggs = ((met <= 350. && ptHiggs > 100) || (met > 350. && ptHiggs > 300));

  // fill cutflows
  bool passCutFlow_Preselection = true;
  bool passCutFlow_Merged = true;
  bool passCutFlow_Resolved = true;



  // preselection
  if (true) {
    // reco level cuts which are not implemented at truth level
    fill("Cutflow_merged", 0.1); // Initial
    fill("Cutflow_merged", 1);   // PassGRL
    fill("Cutflow_merged", 2);   // passLArTile
    fill("Cutflow_merged", 3);   // Trigger
    fill("Cutflow_merged", 4);   // HasVtx
    fill("Cutflow_merged", 5);   // BadJet
    fill("Cutflow_merged", 6);   // CosmicMuon
    fill("Cutflow_merged", 7);   // BadMuon
    fill("Cutflow_merged", 8);   // PFlow Electron veto
    fill("Cutflow_merged", 9);   // IsMETTrigPassed

    fill("Cutflow_resolved", 0.1); // Initial
    fill("Cutflow_resolved", 1);   // PassGRL
    fill("Cutflow_resolved", 2);   // passLArTile
    fill("Cutflow_resolved", 3);   // Trigger
    fill("Cutflow_resolved", 4);   // HasVtx
    fill("Cutflow_resolved", 5);   // BadJet
    fill("Cutflow_resolved", 6);   // CosmicMuon
    fill("Cutflow_resolved", 7);   // BadMuon
    fill("Cutflow_resolved", 8);   // PFlow Electron veto
    fill("Cutflow_resolved", 9);   // IsMETTrigPassed
  } else{
    passCutFlow_Preselection = false;
  }
  if (passCutFlow_Preselection && evtsel_electronVeto) {
    fill("Cutflow_merged", 10);   // 0 baseline electrons
    fill("Cutflow_resolved", 10); // 0 baseline electrons
  } else{
    passCutFlow_Preselection = false;
  }
  if (passCutFlow_Preselection && evtsel_muonVeto) {
    fill("Cutflow_merged", 11);   // 0 baseline muons
    fill("Cutflow_resolved", 11); // 0 baseline muons
  } else{
    passCutFlow_Preselection = false;
  }
  if (passCutFlow_Preselection && evtsel_tauVeto) {
    fill("Cutflow_merged", 12);   // Tau Veto
    fill("Cutflow_resolved", 12); // Tau Veto
  } else{
    passCutFlow_Preselection = false;
  }
  if (passCutFlow_Preselection && evtsel_extendedTauVeto) {
    fill("Cutflow_merged", 13);   // Extended Tau Veto
    fill("Cutflow_resolved", 13); // Extended Tau Veto
  } else{
    passCutFlow_Preselection = false;
  }

  // merged cutflow
  if (passCutFlow_Preselection && passCutFlow_Merged && evtsel_met500) {
    fill("Cutflow_merged", 14); // MetTST>500
  } else{
    passCutFlow_Merged = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Merged && evtsel_nfatjets) {
    fill("Cutflow_merged", 15); // >=1 fat-jets
  } else{
    passCutFlow_Merged = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Merged && evtsel_nbtrackjets) {
    fill("Cutflow_merged", 16); // >= 2 b-tagged track jets
  } else{
    passCutFlow_Merged = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Merged && evtsel_mass40_merged) {
    fill("Cutflow_merged", 17); // (mJ > 40 || mJ_corr > 40)
  } else{
    passCutFlow_Merged = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Merged && evtsel_massRange_merged) {
    fill("Cutflow_merged", 18); // (mJ> 50 && mJ < 270)
  } else{
    passCutFlow_Merged = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Merged && evtsel_ntrackjets_associated) {
    fill("Cutflow_merged", 19); // N_associated_trkJets>=2
    fill("Cutflow_merged", 20); // (TrackJet_1passOR = true && TrackJet_2passOR = true)
  } else{
    passCutFlow_Merged = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Merged && evtsel_minDPhi20) {
    fill("Cutflow_merged", 21); // |DeltaPhiMin3|>20deg
  } else{
    passCutFlow_Merged = false;
  }

  // resolved cutflow
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_met150) {
    fill("Cutflow_resolved", 14); // MetTST>150
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_metleq500) {
    fill("Cutflow_resolved", 15); // MetTST<=500
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_njets) {
    fill("Cutflow_resolved", 16); // >=2 jets
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_nbjets) {
    fill("Cutflow_resolved", 17); // >=2 b-tags
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_mass40_resolved) {
    fill("Cutflow_resolved", 18); // (mjj > 40 || mjj_corr > 40)
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_massRange_resolved) {
    fill("Cutflow_resolved", 19); // (mjj > 50 && mjj < 280)
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_minDPhi20) {
    fill("Cutflow_resolved", 20); // |DeltaPhiMin3|>20deg
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_metSig) {
    fill("Cutflow_resolved", 21); // METSig>12
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_ptHiggs) {
    fill("Cutflow_resolved", 22); // (pt_jj > 100000 || pt_jj_corr > 100000)
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_mt_mindR) {
    fill("Cutflow_resolved", 23); // mT_METclosestBJet > 170000
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_mt_maxdR) {
    fill("Cutflow_resolved", 24); // mT_METfurthestBJet > 200000
  } else{
    passCutFlow_Resolved = false;
  }
  if (passCutFlow_Preselection && passCutFlow_Resolved && evtsel_njets_max) {
    fill("Cutflow_resolved", 25); // <=4 jets
  } else{
    passCutFlow_Resolved = false;
  }

 
  // Preselection
  if (!passCutFlow_Preselection) return;

  // Merged selection
  passMerged = passCutFlow_Preselection && passCutFlow_Merged;

  // Resolved selection
  passResolved = passCutFlow_Preselection && passCutFlow_Resolved;

  // Fill histograms
  if (passCutFlow_Preselection) {
    fill("MET", met);
  }

  if (passMerged) {
    fill("N_Bjets_merged_inside", nBJetsMerged_inside);
    fill("N_Bjets_merged_outside", nBJetsMerged_outside);
    if (nBJetsMerged_inside == 2 && nBJetsMerged_outside == 0)
      fill("mHiggs_merged", mHiggs);
  }

  if (passResolved) {
    fill("N_Bjets_resolved", nBjets);
    fill("mHiggs_resolved", mHiggs);
  }

  // Fill results
  if (passMerged) {
    if (nBJetsMerged_inside == 2 && nBJetsMerged_outside == 0 && met > 500 && met <= 750) {
      accept("MET500750_2b");
      fill("MET500750_2b", mHiggs);
    }
    if (nBJetsMerged_inside == 2 && nBJetsMerged_outside == 0 && met > 750) {
      accept("MET750_2b");
      fill("MET750_2b", mHiggs);
    }
    if (nBJetsMerged_inside == 2 && nBJetsMerged_outside > 0 && met > 500 ) {
      accept("MET500_3b");
      fill("MET500_3b", mHiggs);
    }
  }

  if (passResolved) {
    if (nBjets == 2 && met > 150 && met <= 200) {
      accept("MET150200_2b");
      fill("MET150200_2b", mHiggs);
    }
    if (nBjets == 2 && met > 200 && met <= 350) {
      accept("MET200350_2b");
      fill("MET200350_2b", mHiggs);
    }
    if (nBjets == 2 && met > 350 && met <= 500) {
      accept("MET350500_2b");
      fill("MET350500_2b", mHiggs);
    }

    if (nBjets >= 3 && met > 150 && met <= 200) {
      accept("MET150200_3b");
      fill("MET150200_3b", mHiggs);
    }
    if (nBjets >= 3 && met > 200 && met <= 350) {
      accept("MET200350_3b");
      fill("MET200350_3b", mHiggs);
    }
    if (nBjets >= 3 && met > 350 && met <= 500) {
      accept("MET350500_3b");
      fill("MET350500_3b", mHiggs);
    }
  }

  return;
}
