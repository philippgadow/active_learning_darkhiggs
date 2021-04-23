#include "SimpleAnalysisFramework/NtupleMaker.h"

// This is a special instance of the analysis class for writing out slimmed
// ntuples

void NtupleMaker::ProcessEvent(AnalysisEvent *event) {
  auto electrons =
      event->getElectrons(_minElecPt, 4.2, 0); // Filter on pT, eta and "ID"
  auto muons = event->getMuons(_minMuonPt, 4.2, 0);
  auto taus = event->getTaus(_minTauPt, 4.2, 0);
  auto photons = event->getPhotons(_minPhotonPt, 4.2, 0);
  auto jets = event->getJets(_minJetPt, 5.2, 0);
  auto fatjets = event->getFatJets(_minFatJetPt, 5.2, 0);
  auto bhadrons = event->getBHadrons(0, 99., 0);
  auto trackjets = event->getTrackJets(_minTrackJetPt, 2.5, 0);
  auto HSTruth = event->getHSTruth();
  auto met = event->getMET();
  float sumet = event->getSumET();

  ntupVar("el", electrons);
  ntupVar("mu", muons);
  ntupVar("tau", taus);
  ntupVar("ph", photons);
  ntupVar("jet", jets, true);
  ntupVar("fatjet", fatjets, true);
  ntupVar("bhadron", bhadrons, true);
  ntupVar("trackjet", trackjets, true);
  ntupVar("HSTruth", HSTruth, true, false, true);
  ntupVar("met", met);
  ntupVar("sumet", sumet);

  ntupVar("mcChannel", event->getMCNumber());
  ntupVar("mcVetoCode", int(event->getMCVeto()));
  ntupVar("susyChannel", event->getSUSYChannel());
  ntupVar("mcWeights", event->getMCWeights());
  ntupVar("genMET", event->getGenMET());
  ntupVar("genHT", event->getGenHT());
  ntupVar("pdf_id1", event->getPDF_id1());
  ntupVar("pdf_x1", event->getPDF_x1());
  ntupVar("pdf_pdf1", event->getPDF_pdf1());
  ntupVar("pdf_id2", event->getPDF_id2());
  ntupVar("pdf_x2", event->getPDF_x2());
  ntupVar("pdf_pdf2", event->getPDF_pdf2());
  ntupVar("pdf_scale", event->getPDF_scale());
  return;
}
