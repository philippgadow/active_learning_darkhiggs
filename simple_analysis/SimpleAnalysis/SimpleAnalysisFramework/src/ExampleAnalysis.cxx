#include "SimpleAnalysisFramework/AnalysisClass.h"

DefineAnalysis(ExampleAnalysis)

    void ExampleAnalysis::Init() {
  // Define signal/control regions
  addRegions({"SR1", "SR2"}); // Use same name as in paper, but avoid "-"
  addRegion("CR_B");          // and similar non-alphanumeric characters

  // Book 1/2D histograms
  addHistogram("MET", 100, 0, 2000);
  addHistogram("MTvsMET", 100, 0, 2000, 100, 0, 300);
}

void ExampleAnalysis::ProcessEvent(AnalysisEvent *event) {
  // Retrieve basic object lists
  // PLEASE NOTE UNITS ARE ALL IN GEV, AND NOT ATLAS STANDARD MEV!!!
  auto electrons =
      event->getElectrons(20, 2.47, ELooseLH); // Filter on pT, eta and "ID"
  auto muons = event->getMuons(20, 2.5, MuMedium);
  auto candJets = event->getJets(20., 2.8);
  auto metVec = event->getMET();
  double metSignificance =
      event->getMETSignificance();    // using object pT resolutions
  auto HSTruth = event->getHSTruth(); // for accessing BSM particles
  auto Charginos = event->getHSTruth(
      20, 5,
      Chargino1); // for accessing specific BSM particle with Pt>20 and |eta|<5
  double met = metVec.Et();

  // can get generator MET and HT (if available) to handle filters
  double generatorMET = event->getGenMET();
  double generatorHT = event->getGenHT();
  if (generatorMET < 100 || generatorHT > 10000)
    adjustEventWeight(
        10.4); // weight is multiplied to the default MC weight
               // it should be called before histogram filling and pass() calls

  // Fill in histogram
  fill("MET", met);

  // Overlap removal - including with object Pt-dependent radius calculation
  auto radiusCalcJet = [](const AnalysisObject &, const AnalysisObject &muon) {
    return std::min(0.4, 0.04 + 10 / muon.Pt());
  };
  auto radiusCalcMuon = [](const AnalysisObject &muon, const AnalysisObject &) {
    return std::min(0.4, 0.04 + 10 / muon.Pt());
  };
  electrons = overlapRemoval(electrons, muons, 0.01);
  candJets = overlapRemoval(candJets, electrons, 0.2, NOT(BTag85MV2c20));
  electrons = overlapRemoval(electrons, candJets, 0.4);
  candJets = overlapRemoval(candJets, muons, radiusCalcJet, LessThan3Tracks);
  muons = overlapRemoval(muons, candJets, radiusCalcMuon);

  // Jets can be reclusters into fat jets
  auto fatJets =
      reclusterJets(candJets, 1.0, 300, 0.2,
                    0.05); // input objects, radius, min pT, rclus and ptfrac
  fatJets = filterObjects(fatJets, 300, 2.0);

  // Basic filtering by pT, eta and "ID"
  auto signalJets = filterObjects(candJets, 30);
  auto signalElectrons = filterObjects(
      electrons, 20, 2.47, ETightLH | ED0Sigma5 | EZ05mm | EIsoBoosted);
  auto signalMuons = filterObjects(
      muons, 20, 2.5, MuD0Sigma3 | MuZ05mm | MuIsoBoosted | MuNotCosmic);
  auto bjets = filterObjects(signalJets, 30., 2.5, BTag85MV2c20);

  // Lists of objects can be merget by simple addition
  auto signalLeptons = signalElectrons + signalMuons;

  // Object counting
  int numSignalLeptons =
      signalJets
          .size(); // Object lists are essentially std::vectors so .size() works
  int nBjets = bjets.size();
  int numSignalJets50 = countObjects(
      candJets, 50,
      2.4); // Shorthand to avoid making new list when only number is important

  // Calculate a few kinematic variables - add your favorites to
  // AnalysisClass.cxx
  float meff4j =
      met + sumObjectsPt(signalJets,
                         4); // Meff from MET and (up to) four leading jets
  float meffIncl =
      met + sumObjectsPt(signalJets) +
      sumObjectsPt(signalLeptons); // Meff from MET, jets and leptons
  float mT = 0.;
  if (signalLeptons.size() > 0)
    mT = calcMT(signalLeptons[0], metVec); // MT from leading lepton
  float mTbmin =
      calcMTmin(bjets, metVec, 3); // Minimum MT of leading three b-jets
  float dphiMin4 =
      minDphi(metVec, signalJets,
              4); // Smallest Delta Phi between leading 4-jets and MET

  // Preselection
  if (numSignalJets50 < 2)
    return;

  // Fill another histogram
  fill("MTvsMET", met, mT);

  // Flag signal regions that pass selections
  if (met > 200 && mT > 100 && nBjets > 2)
    accept("SR1");
  if (meff4j > 1000 && mTbmin > 200 && dphiMin4 > 1.0)
    accept("SR2");
  if (numSignalLeptons == 2 && dphiMin4 < 0.4)
    accept("CR_B");

  // Access to BSM particles
  float neutralinoSum = 0;
  for (auto bsmParticle : HSTruth) {
    // Note that in most samples the BSM particles are stored multiple times, so
    // one should check either parent type or status code
    if (abs(bsmParticle.pdgId()) == 1000022 &&
        bsmParticle.status() == StablePart)
      neutralinoSum += bsmParticle.Pt();
    // For BSM particles the production and decay vertex is stored if available
    // otherwise set to (0,0,0)
    if (bsmParticle.prodVtx().Mag2() != 0) {
      std::cout << "Prod vertex: (" << bsmParticle.prodVtx().X() << ","
                << bsmParticle.prodVtx().Y() << "," << bsmParticle.prodVtx().Z()
                << ")" << std::endl;
    }
    if (bsmParticle.decayVtx().Mag2() != 0) {
      std::cout << "Decay vertex: (" << bsmParticle.decayVtx().X() << ","
                << bsmParticle.decayVtx().Y() << ","
                << bsmParticle.decayVtx().Z() << ")" << std::endl;
    }
  }

  // Access to parentID
  //  - for now only available for muons, electrons and photons and BSM
  //  particles and depends on input
  int leptonsFromZ = 0;
  for (auto lepton : signalLeptons) {
    if (lepton.motherID() == 23)
      leptonsFromZ++;
  }

  // Fill in optional ntuple variables
  ntupVar("met", met); // can be simple variables
  ntupVar("metSig", metSignificance);
  ntupVar("meff", meffIncl);
  ntupVar("nBjets", nBjets);
  ntupVar("signalJets",
          signalJets); // or even a list of objects (or single object)
  if (nBjets > 0)
    ntupVar("leadingBjet", bjets[0]);

  // The following calls might disappear again once their use
  // is better understood and can be incorporated in the framework
  ntupVar("susyProcess", event->getSUSYChannel());
  ntupVar("mcDSID", event->getMCNumber());
  ntupVar("mcWeights", event->getMCWeights());

  // example of accessing truth information - this only works for smeared events
  std::vector<float> jetpt, truthjetpt;
  for (auto jet : candJets) {
    jetpt.push_back(jet.Pt());
    auto truth = event->getTruthParticle(jet);
    if (truth.valid()) {
      truthjetpt.push_back(truth.Pt());
    } else
      truthjetpt.push_back(0);
  }
  ntupVar("jetpt", jetpt);
  ntupVar("truthjetpt", truthjetpt);
  auto truthMet = event->getTruthParticle(metVec);
  ntupVar("truthMET", truthMet);
  return;
}
