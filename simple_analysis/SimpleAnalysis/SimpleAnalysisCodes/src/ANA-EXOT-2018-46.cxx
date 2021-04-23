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
  int nBJetsMerged_inside = 0;
  int nBJetsMerged_outside = 0;
  if (met > 500) {
    if (fatJets.size() > 0) {
      mHiggs = fatJets[0].M();
      ptHiggs = fatJets[0].Pt();
      for (const auto &jet : btrackjets) {
        if (jet.DeltaR(fatJets[0]) < 1.)
          nBJetsMerged_inside++;
        else
          nBJetsMerged_outside++;        
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
  bool evtsel_met150 = (met > 150.);
  bool evtsel_leptonVeto = (baselineLeptons.size() == 0);
  bool evtsel_tauVeto = (baselineTaus.size() == 0);
  bool evtsel_extendedTauVeto = true; // not implemented on truth level
  bool evtsel_minDPhi20 = (dphiMin3 > 20. * M_PI / 180.);

  // cuts in merged selection
  bool evtsel_met500 = (met > 500.);
  bool evtsel_massRange_merged = (mHiggs > 50. && mHiggs < 270.);

  // cuts in resolved selection
  bool evtsel_metleq500 = (met <= 500.);
  bool evtsel_njets = (centralJets.size() >=2);
  bool evtsel_nbjets = (nBjets >=2);
  bool evtsel_ptHiggs = ((met <= 350. && ptHiggs > 100) || (met > 350. && ptHiggs > 300));
  bool evtsel_mt_mindR = (mt_min > 170.);
  bool evtsel_mt_maxdR = (mt_max > 200.);
  bool evtsel_metSig = (metSignificance > 12.);
  bool evtsel_njets_max = ((nBjets == 2 && centralJets.size() <= 4) || (nBjets >=3 && centralJets.size() <= 5));
  bool evtsel_massRange_resolved = (mHiggs > 50. && mHiggs < 280.);

  // Preselection
  if (!evtsel_met150) return;
  if (!evtsel_leptonVeto) return;
  if (!evtsel_tauVeto) return;
  if (!evtsel_extendedTauVeto) return;
  if (!evtsel_minDPhi20) return;

  // Merged selection
  passMerged = evtsel_met500 && evtsel_massRange_merged;

  // Resolved selection
  passResolved = evtsel_metleq500 && evtsel_njets && evtsel_nbjets && evtsel_ptHiggs && \
                 evtsel_mt_mindR && evtsel_mt_maxdR && evtsel_metSig && evtsel_njets_max && \
                 evtsel_massRange_resolved;

  // Fill histograms
  if (passMerged || passResolved) {
    fill("MET", met);
  }

  if (passMerged) {
    fill("N_Bjets_merged_inside", nBJetsMerged_inside);
    fill("N_Bjets_merged_outside", nBJetsMerged_outside);
    if (nBJetsMerged_inside == 2)
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
