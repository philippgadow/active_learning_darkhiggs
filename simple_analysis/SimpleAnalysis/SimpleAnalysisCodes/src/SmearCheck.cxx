#include "SimpleAnalysisFramework/AnalysisClass.h"

/* Simple routine for checking basic distribution - primarily for checking
 * smearing functions */

DefineAnalysis(SmearCheck) std::vector<std::string> names = {
    "MuonLoose",      "MuonTight",     "MuonHighPt", "ElectronLoose",
    "ElectronMedium", "ElectronTight", "TauLoose",   "TauMedium",
    "TauTight",       "Photon",        "Jet",        "BJet85",
    "BJet70"};

void SmearCheck::Init() {
  // Book 1/2D histograms
  addHistogram("MET", 100, 0, 2000);

  for (unsigned int ii = 0; ii < names.size(); ii++) {
    addHistogram(names[ii] + "Size", 10, 0, 10);
    addHistogram(names[ii] + "PT", 1000, 0, 1000);
    addHistogram(names[ii] + "Eta", 100, -5., 5.);
  }
}

void SmearCheck::ProcessEvent(AnalysisEvent *event) {
  auto electrons =
      event->getElectrons(2, 5., ELooseLH); // Filter on pT, eta and "ID"
  auto elecMed =
      event->getElectrons(2, 5., EMediumLH); // Filter on pT, eta and "ID"
  auto elecTight =
      event->getElectrons(2, 5., ETightLH); // Filter on pT, eta and "ID"
  auto muons = event->getMuons(2, 5., MuLoose);
  auto muonTight = event->getMuons(2, 5., MuTight);
  auto muonHighPt = event->getMuons(2, 5., MuHighPt);
  auto taus = event->getTaus(2, 5., TauLoose);
  auto tauMedium = event->getTaus(2, 5., TauMedium);
  auto tauTight = event->getTaus(2, 5., TauTight);
  auto photons = event->getPhotons(2, 5., PhotonLoose);
  auto jets = event->getJets(20., 5.);
  auto bjets85 = event->getJets(20., 5., BTag85MV2c10);
  auto bjets70 = event->getJets(20., 5., BTag70MV2c10);
  auto metVec = event->getMET();
  auto HSTruth = event->getHSTruth(); // for accessing BSM particles
  double met = metVec.Et();

  fill("MET", met);
  std::vector<AnalysisObjects *> objs = {
      &muons,     &muonTight, &muonHighPt, &electrons, &elecMed,
      &elecTight, &taus,      &tauMedium,  &tauTight,  &photons,
      &jets,      &bjets85,   &bjets70};
  for (unsigned int ii = 0; ii < names.size(); ii++) {
    fill(names[ii] + "Size", objs[ii]->size());
    for (const auto obj : *(objs[ii])) {
      fill(names[ii] + "PT", obj.Pt());
      fill(names[ii] + "Eta", obj.Eta());
    }
  }

  return;
}
