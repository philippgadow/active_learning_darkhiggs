#include <algorithm>
#include <cmath>
#include <iostream>

#include "SUSYTools/SUSYCrossSection.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "xAODCore/AuxContainerBase.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruthReader.h"
#include <TH2.h>
#include <TStyle.h>

#include "MCEventVetoHelper.h"
#include "SimpleAnalysisFramework/OutputHandler.h"
#include "SimpleAnalysisFramework/TruthEvent.h"

using std::vector;

#define WARN_ONCE(warning)                                                     \
  do {                                                                         \
    static bool first = true;                                                  \
    if (first)                                                                 \
      std::cout << warning << std::endl;                                       \
    first = false;                                                             \
  } while (0)

DefineReader(xAODTruthReader);

xAODTruthReader::xAODTruthReader() : Reader() {
  if (!xAOD::Init().isSuccess()) {
    throw std::runtime_error("Cannot initialise xAOD access !");
  }
}

void xAODTruthReader::AddOptions(po::options_description &desc) {
  desc.add_options()("useTrueTau,T",
                     "use true tau 4-vector instead of visible tau 4-vector")(
      "ignoreTruthBSM", "ignore BSM truth blocks and use directly "
                        "TruthParticles (needs TRUTH1 input)");
}

bool xAODTruthReader::isFileType(TFile *fh, po::variables_map &vm) {
  if (fh->FindKey("CollectionTree") && !vm.count("readReco")) {
    std::cout << "Reading xAOD truth input" << std::endl;
    _useVisTau = vm.count("useTrueTau") == 0;
    _useTruthBSM = vm.count("ignoreTruthBSM") == 0;
    return true;
  }
  return false;
}

void xAODTruthReader::Init() {
  _event = new xAOD::TEvent(xAOD::TEvent::kClassAccess);
  _susytools = new ST::SUSYObjDef_xAOD("mySUSYTools");
  _mctool = new MCTruthClassifier("myTruthFinder");
}

int xAODTruthReader::getTruthOrigin(const xAOD::TruthParticle *part) {
  if (part->isAvailable<unsigned int>("classifierParticleOrigin")) {
    return part->auxdata<unsigned int>("classifierParticleOrigin");
  }
  if (part->isAvailable<unsigned int>("particleOrigin")) {
    return part->auxdata<unsigned int>("particleOrigin");
  }

  const ElementLink<xAOD::TruthParticleContainer> origPart =
      part->auxdata<ElementLink<xAOD::TruthParticleContainer>>(
          "originalTruthParticle");
  if (origPart.isValid()) {
    const auto result = _mctool->particleTruthClassifier(*origPart);
    return result.second;
  }
  return 0;
}

int xAODTruthReader::getTruthType(const xAOD::TruthParticle *part) {
  if (part->isAvailable<unsigned int>("classifierParticleType")) {
    return part->auxdata<unsigned int>("classifierParticleType");
  }
  if (part->isAvailable<unsigned int>("particleType")) {
    return part->auxdata<unsigned int>("particleType");
  }

  if (part->isAvailable<ElementLink<xAOD::TruthParticleContainer>>(
          "originalTruthParticle")) {
    const ElementLink<xAOD::TruthParticleContainer> origPart =
        part->auxdata<ElementLink<xAOD::TruthParticleContainer>>(
            "originalTruthParticle");
    if (origPart.isValid()) {
      const auto result = _mctool->particleTruthClassifier(*origPart);
      return result.first;
    }
  }
  return 0;
}

static SG::AuxElement::Accessor<int> acc_motherID("motherID");

int xAODTruthReader::getMotherID(const xAOD::TruthParticle *part,
                                 bool traverseProdChain) {
  if (acc_motherID.isAvailable(*part))
    return acc_motherID(*part);

  if (part->nParents() > 0) {
    auto p = part->parent(0);
    if (p) {
      if (traverseProdChain && (p->pdgId() == part->pdgId())) {
        return getMotherID(p, true);
      }

      return p->pdgId();
    }
  }

  return 0;
}

xAOD::TruthParticleContainer *xAODTruthReader::findTruthParticles(
    xAOD::TStore *store, const xAOD::TruthParticleContainer *truthparticles,
    std::vector<int> pdgIds, int status) {
  xAOD::TruthParticleContainer *truth = new xAOD::TruthParticleContainer;
  xAOD::AuxContainerBase *truthAux = new xAOD::AuxContainerBase();
  truth->setStore(truthAux);
  SG::AuxElement::Decorator<ElementLink<xAOD::TruthParticleContainer>>
      linkDecorator("originalTruthParticle");
  int idx = 0;
  for (xAOD::TruthParticleContainer::const_iterator it =
           truthparticles->begin();
       it != truthparticles->end(); ++it) {

    if (std::find(pdgIds.begin(), pdgIds.end(), abs((*it)->pdgId())) !=
            pdgIds.end() &&
        (*it)->status() == status && (*it)->barcode() < 200000) {
      const auto result = _mctool->particleTruthClassifier(*it);
      MCTruthPartClassifier::ParticleOutCome outcome =
          _mctool->getParticleOutCome();
      if (outcome == MCTruthPartClassifier::DecaytoElectron ||
          outcome == MCTruthPartClassifier::DecaytoMuon)
        continue;
      xAOD::TruthParticle *part = new xAOD::TruthParticle();
      truth->push_back(part);
      *part = **it;
      part->auxdecor<unsigned int>("classifierParticleType") = result.first;
      ElementLink<xAOD::TruthParticleContainer> eltp(*truthparticles, idx);
      linkDecorator(*part) = eltp;
      if (status == 2) { // tau
        int numPart = 0;
        if (outcome == MCTruthPartClassifier::OneProng)
          numPart = 1;
        if (outcome == MCTruthPartClassifier::ThreeProng)
          numPart = 3;
        if (outcome == MCTruthPartClassifier::FiveProng)
          numPart = 5;
        part->auxdecor<char>("IsHadronicTau") = (char)1;
        part->auxdecor<size_t>("numCharged") = numPart;
      }
      if ((*it)->hasProdVtx()) {
        if ((*it)->prodVtx()->nIncomingParticles() > 0)
          part->auxdecor<int>("motherID") =
              (*it)->prodVtx()->incomingParticle(0)->pdgId();
      }
    }
    idx++;
  }

  std::string pid;
  for (auto ipdg : pdgIds)
    pid += std::to_string(ipdg);

  if (!store->record(truth, (std::string("Good") + pid).c_str()).isSuccess())
    throw std::runtime_error("Could not record truth particles");
  if (!store->record(truthAux, (std::string("Good") + pid + "Aux.").c_str())
           .isSuccess())
    throw std::runtime_error("Could not record truth particles Aux");
  return truth;
}

xAOD::TruthParticleContainer *xAODTruthReader::findTruthBSMParticles(
    xAOD::TStore *store, const xAOD::TruthParticleContainer *truthparticles) {
  xAOD::TruthParticleContainer *truth = new xAOD::TruthParticleContainer;
  xAOD::AuxContainerBase *truthAux = new xAOD::AuxContainerBase();
  truth->setStore(truthAux);
  SG::AuxElement::Decorator<ElementLink<xAOD::TruthParticleContainer>>
      linkDecorator("originalTruthParticle");
  int idx = 0;
  for (xAOD::TruthParticleContainer::const_iterator it =
           truthparticles->begin();
       it != truthparticles->end(); ++it) {
    int pdgId = abs((*it)->pdgId());
    if ((31 < pdgId && pdgId < 38) || pdgId == 39 || pdgId == 41 ||
        pdgId == 42 || pdgId == 8 || (1000000 < pdgId && pdgId < 1000040) ||
        (2000000 < pdgId && pdgId < 2000040) ||
        (3000000 < pdgId && pdgId < 3000040)) {
      const auto result = _mctool->particleTruthClassifier(*it);
      //      MCTruthPartClassifier::ParticleOutCome
      //      outcome=_mctool->getParticleOutCome();
      xAOD::TruthParticle *part = new xAOD::TruthParticle();
      truth->push_back(part);
      *part = **it;
      part->auxdecor<unsigned int>("classifierParticleType") = result.first;
      ElementLink<xAOD::TruthParticleContainer> eltp(*truthparticles, idx);
      linkDecorator(*part) = eltp;
      if ((*it)->hasProdVtx()) {
        if ((*it)->prodVtx()->nIncomingParticles() > 0)
          part->auxdecor<int>("motherID") =
              (*it)->prodVtx()->incomingParticle(0)->pdgId();
      }
    }
    idx++;
  }

  if (!store->record(truth, "GoodBSM").isSuccess())
    throw std::runtime_error("Could not record truth particles");
  if (!store->record(truthAux, "GoodBSMAux.").isSuccess())
    throw std::runtime_error("Could not record truth particles Aux");
  return truth;
}

static SG::AuxElement::Accessor<float> acc_filtHT("GenFiltHT");
static SG::AuxElement::Accessor<float> acc_filtMET("GenFiltMET");

bool xAODTruthReader::processEvent(xAOD::TEvent *xaodEvent,
                                   xAOD::TStore *store) {

  const xAOD::EventInfo *eventInfo = 0;
  if (!xaodEvent->retrieve(eventInfo, "EventInfo").isSuccess()) {
    throw std::runtime_error("Cannot read EventInfo");
  }

  int eventNumber = eventInfo->eventNumber();
  int mcChannel = eventInfo->mcChannelNumber();
  if (mcChannel == 0)
    mcChannel = eventInfo->runNumber();
  int susy_part_id1 = 0;
  int susy_part_id2 = 0;
  int susy_process = 0;

  const xAOD::TruthParticleContainer *truthparticles = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthParticles")) {
    if (!xaodEvent->retrieve(truthparticles, "TruthParticles").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthParticles");
    }
  } else {
    if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthBSM")) {
      if (!xaodEvent->retrieve(truthparticles, "TruthBSM").isSuccess()) {
        throw std::runtime_error(
            "Could not retrieve truth particles with key TruthBSM");
      }
    } else {
      xAOD::TruthParticleContainer *emptyContainer =
          new xAOD::TruthParticleContainer;
      xAOD::AuxContainerBase *emptyAux = new xAOD::AuxContainerBase();
      emptyContainer->setStore(emptyAux);
      truthparticles = emptyContainer;
      WARN_ONCE("Warning: No TruthParticles or TruthBSM container in input - "
                "can't identify process");
    }
  }
  _susytools->FindSusyHardProc(truthparticles, susy_part_id1, susy_part_id2);
  if (susy_part_id2 == 0 && truthparticles->size() > 1)
    susy_part_id2 = truthparticles->at(1)->pdgId();
  if (susy_part_id1 == 0 && truthparticles->size())
    susy_part_id1 = truthparticles->at(0)->pdgId();
  if ((abs(susy_part_id1) > 1000000) && (abs(susy_part_id1) > 1000000) &&
      (abs(susy_part_id1) < 3000000) &&
      (abs(susy_part_id1) < 3000000)) { // only consider SUSY  BSM particles
    int maxID =
        std::max(abs(susy_part_id1) % 1000000, abs(susy_part_id2) % 1000000);
    int minID =
        std::min(abs(susy_part_id1) % 1000000, abs(susy_part_id2) % 1000000);
    if (maxID <= 6 ||
        minID >
            6) // no final state defined for mixed squark+neutralino production
      susy_process = SUSY::finalState(susy_part_id1, susy_part_id2);
    else
      WARN_ONCE("Warning: Sample includes mixed squark-neutralino production!\n"
                "         No finalstate number is defined for this - harmless "
                "unless SUSY finalstate is used");
  }

  const xAOD::MissingETContainer *metCont = 0;
  if (!xaodEvent->retrieve(metCont, "MET_Truth").isSuccess()) {
    throw std::runtime_error("Could not retrieve truth met with key MET_Truth");
  }
  const xAOD::MissingET *met = (*metCont)["NonInt"];
  const xAOD::MissingET *sumet = (*metCont)["Int"];

  bool *mcaccept = new bool(true);
  unsigned int *veto = new unsigned int(0);
  float filtMET = 0., filtHT = 0.;

  filtMET = eventInfo->auxdata<float>("GenFiltMET");
  filtHT = eventInfo->auxdata<float>("GenFiltHT");

  *mcaccept = MCEventVetoHelper::mc16accept(*veto, mcChannel, filtMET, filtHT,
                                            true, truthparticles);
  xAOD::TStore *vetostore = xAOD::TActiveStore::store();
  RETURN_CHECK("MCEventVeto::processEvent",
               vetostore->record<bool>(mcaccept, "mcAccept"));
  RETURN_CHECK("MCEventVeto::processEvent",
               vetostore->record<unsigned int>(veto, "mcVetoCode"));

  TruthEvent *event = new TruthEvent(sumet->sumet() / 1000., met->mpx() / 1000.,
                                     met->mpy() / 1000.);
  event->setChannelInfo(mcChannel, *veto, susy_process);

  TLorentzVector tlv(0., 0., 0., 0.);

  int idx = 0;
  const xAOD::TruthParticleContainer *truthelectrons = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthElectrons")) {
    if (!xaodEvent->retrieve(truthelectrons, "TruthElectrons").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthElectrons");
    }
  } else {
    truthelectrons = findTruthParticles(store, truthparticles, {11});
  }

  if (!truthelectrons->getConstStore()) {
    WARN_ONCE("No Aux store found for TruthElectrons\n - cannot read "
              "electrons\n - check input if needed");
  } else {
    for (xAOD::TruthParticleContainer::const_iterator it =
             truthelectrons->begin();
         it != truthelectrons->end(); ++it) {
      const auto electron = *it;
      int iso = getTruthType(electron) == MCTruthPartClassifier::IsoElectron;
      tlv.SetPtEtaPhiM(electron->pt() / 1000., electron->eta(), electron->phi(),
                       electron->m() / 1000.);
      event->addElectron(tlv, electron->charge(), iso ? EIsoGood : 0,
                         getMotherID(electron), idx++);
    }
  }

  idx = 0;
  const xAOD::TruthParticleContainer *truthmuons = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthMuons")) {
    if (!xaodEvent->retrieve(truthmuons, "TruthMuons").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthMuons");
    }
  } else {
    truthmuons = findTruthParticles(store, truthparticles, {13});
  }
  for (xAOD::TruthParticleContainer::const_iterator it = truthmuons->begin();
       it != truthmuons->end(); ++it) {
    const auto muon = *it;
    int iso = getTruthType(muon) == MCTruthPartClassifier::IsoMuon;
    tlv.SetPtEtaPhiM(muon->pt() / 1000., muon->eta(), muon->phi(),
                     muon->m() / 1000.);
    event->addMuon(tlv, muon->charge(), iso ? MuIsoGood : 0, getMotherID(muon),
                   idx++);
  }

  idx = 0;
  const xAOD::TruthParticleContainer *truthtaus = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthTaus")) {
    if (!xaodEvent->retrieve(truthtaus, "TruthTaus").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthTaus");
    }
  } else {
    truthtaus = findTruthParticles(store, truthparticles, {15}, 2);
  }
  for (xAOD::TruthParticleContainer::const_iterator it = truthtaus->begin();
       it != truthtaus->end(); ++it) {
    const auto tau = *it;
    if (tau->auxdata<char>("IsHadronicTau")) {
      int iso = getTruthType(tau) == MCTruthPartClassifier::IsoTau;

      if (_useVisTau) {
        WARN_ONCE("Using visible tau 4-vector - to switch back to full truth "
                  "4-vector use "
                  "-T"
                  " option");
        tlv.SetPtEtaPhiM(tau->auxdata<double>("pt_vis") / 1000.,
                         tau->auxdata<double>("eta_vis"),
                         tau->auxdata<double>("phi_vis"),
                         tau->auxdata<double>("m_vis") / 1000.);
      } else {
        tlv.SetPtEtaPhiM(tau->pt() / 1000., tau->eta(), tau->phi(),
                         tau->m() / 1000.);
      }
      int tauId = iso ? TauIsoGood : 0;
      if (tau->auxdata<unsigned long>("numCharged") == 3)
        tauId |= TauThreeProng;
      else
        tauId |= TauOneProng;
      event->addTau(tlv, tau->charge(), tauId, getMotherID(tau), idx++);
    }
  }

  idx = 0;
  const xAOD::TruthParticleContainer *truthphotons = 0;
  std::string photonName = "TruthPhotons";
  if (!xaodEvent->contains<xAOD::TruthParticleContainer>(photonName))
    photonName = "Truth3Photons";

  if (xaodEvent->contains<xAOD::TruthParticleContainer>(photonName)) {
    if (!xaodEvent->retrieve(truthphotons, photonName).isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthPhotons");
    }
  } else {
    truthphotons = findTruthParticles(store, truthparticles, {22});
  }

  for (xAOD::TruthParticleContainer::const_iterator it = truthphotons->begin();
       it != truthphotons->end(); ++it) {
    const auto photon = *it;
    int iso = getTruthType(photon) == MCTruthPartClassifier::IsoPhoton;
    tlv.SetPtEtaPhiM(photon->pt() / 1000., photon->eta(), photon->phi(), 0);
    event->addPhoton(tlv, iso ? PhotonIsoGood : 0, getMotherID(photon), idx++);
  }

  // Generator Filter HT (e.g. for ttbar/singleTop samples)
  float gen_ht = 0.;
  if (acc_filtHT.isAvailable(*(eventInfo))) {
    gen_ht = eventInfo->auxdata<float>("GenFiltHT");
  } else {
    WARN_ONCE("Warning : No GenFiltHT decoration available. Setting HT to 0 "
              "for now...");
  }
  event->setGenHT(gen_ht / 1000.);

  // Generator Filter MET (e.g. for ttbar/singleTop samples)
  float gen_met = 0.;
  if (acc_filtMET.isAvailable(*(eventInfo))) {
    gen_met = eventInfo->auxdata<float>("GenFiltMET");
  } else { // recompute from particle containers!
    idx = 0;
    const xAOD::TruthParticleContainer *truthneutrinos = 0;
    std::string neutrinoName = "TruthNeutrinos";
    if (!xaodEvent->contains<xAOD::TruthParticleContainer>(neutrinoName))
      neutrinoName = "TruthNeutrinos";

    if (xaodEvent->contains<xAOD::TruthParticleContainer>(neutrinoName)) {
      if (!xaodEvent->retrieve(truthneutrinos, neutrinoName).isSuccess()) {
        throw std::runtime_error(
            "Could not retrieve truth particles with key TruthNeutrinos");
      }
    } else {
      truthneutrinos = findTruthParticles(store, truthparticles, {12, 14, 16});
    }

    tlv.SetPtEtaPhiM(0., 0., 0., 0.);
    for (xAOD::TruthParticleContainer::const_iterator it =
             truthneutrinos->begin();
         it != truthneutrinos->end(); ++it) {
      const auto nu = *it;
      int iPartOrig = getTruthOrigin(nu);

      switch (iPartOrig) {
      case MCTruthPartClassifier::PhotonConv:
      case MCTruthPartClassifier::DalitzDec:
      case MCTruthPartClassifier::ElMagProc:
      case MCTruthPartClassifier::Mu:
      case MCTruthPartClassifier::TauLep:
      case MCTruthPartClassifier::LightMeson:
      case MCTruthPartClassifier::StrangeMeson:
      case MCTruthPartClassifier::CharmedMeson:
      case MCTruthPartClassifier::BottomMeson:
      case MCTruthPartClassifier::CCbarMeson:
      case MCTruthPartClassifier::JPsi:
      case MCTruthPartClassifier::BBbarMeson:
      case MCTruthPartClassifier::LightBaryon:
      case MCTruthPartClassifier::StrangeBaryon:
      case MCTruthPartClassifier::CharmedBaryon:
      case MCTruthPartClassifier::BottomBaryon:
      case MCTruthPartClassifier::PionDecay:
      case MCTruthPartClassifier::KaonDecay:

      case MCTruthPartClassifier::NonDefined:
        continue;
      default:
        break;
      }
      tlv += nu->p4();
    }
    gen_met = tlv.Pt();
  }
  event->setGenMET(gen_met / 1000.);

  // std::cout << "MCEventVeto::processEvent GenFiltMET = " << filtMET << ",
  // GenFiltHT = " << filtHT << std::endl;

  // sample overlap removal

  idx = 0;
  const xAOD::JetContainer *truthjets = 0;
  std::string jetName = "AntiKt4TruthJets";
  if (!xaodEvent->contains<xAOD::JetContainer>(jetName))
    jetName = "AntiKt4TruthDressedWZJets";
  if (!xaodEvent->retrieve(truthjets, jetName).isSuccess()) {
    throw std::runtime_error("Could not retrieve truth particles with key " +
                             jetName);
  }
  for (xAOD::JetContainer::const_iterator it = truthjets->begin();
       it != truthjets->end(); ++it) {
    const auto jet = *it;
    tlv.SetPtEtaPhiM(jet->pt() / 1000., jet->eta(), jet->phi(),
                     jet->m() / 1000.);
    int flavor = 0;
    if (jet->isAvailable<int>("HadronConeExclTruthLabelID"))
      flavor = jet->auxdata<int>("HadronConeExclTruthLabelID");
    else if (jet->isAvailable<int>("ConeTruthLabelID"))
      flavor = jet->auxdata<int>("ConeTruthLabelID");
    else if (jet->isAvailable<int>("PartonTruthLabelID"))
      flavor = abs(jet->auxdata<int>("PartonTruthLabelID"));
    else if (jet->isAvailable<int>("GhostBHadronsFinalCount")) {
      if (jet->auxdata<int>("GhostBHadronsFinalCount")) {
        flavor = 5;
      } else if (jet->auxdata<int>("GhostCHadronsFinalCount")) {
        flavor = 4;
      } else
        flavor = 1;
    }
    int id = (flavor == 5) ? GoodBJet : GoodJet;
    if (flavor == 4)
      id |= TrueCJet;
    else if (flavor == 5)
      id |= TrueBJet;
    else if (flavor == 15)
      id |= TrueTau;
    else
      id |= TrueLightJet;
    event->addJet(tlv, id, idx++);
  }

  idx = 0;
  const xAOD::JetContainer *truthfatjets = 0;
  std::string fatjetName = "AntiKt10TruthTrimmedPtFrac5SmallR20Jets";
  if (!xaodEvent->contains<xAOD::JetContainer>(fatjetName))
    fatjetName = "TrimmedAntiKt10TruthJets";

  if (xaodEvent->contains<xAOD::JetContainer>(fatjetName)) {
    if (!xaodEvent->retrieve(truthfatjets, fatjetName).isSuccess()) {
      throw std::runtime_error("Could not retrieve truth particles with key " +
                               fatjetName);
    }
    for (xAOD::JetContainer::const_iterator it = truthfatjets->begin();
         it != truthfatjets->end(); ++it) {
      const auto jet = *it;
      tlv.SetPtEtaPhiM(jet->pt() / 1000., jet->eta(), jet->phi(),
                       jet->m() / 1000.);
      int flavor = 0; // FIXME check if there are more recent labels for fat
                      // jets
      if (jet->isAvailable<int>("PartonTruthLabelID"))
        flavor = abs(jet->auxdata<int>("PartonTruthLabelID"));
      else {
        if (jet->isAvailable<int>("GhostBHadronsFinalCount")) {
          if (jet->auxdata<int>("GhostBHadronsFinalCount")) {
            flavor = 5;
          } else if (jet->auxdata<int>("GhostCHadronsFinalCount")) {
            flavor = 4;
          } else
            flavor = 1;
        }
      }
      event->addFatJet(tlv, (flavor == 5) ? GoodBJet : GoodJet, idx++);
    }
  }

  idx = 0;
  const xAOD::TruthParticleContainer *truthBHadrons = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthParticles")) {
    if (!xaodEvent->retrieve(truthBHadrons, "TruthParticles").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth hadrons with key TruthParticles");
    }
  } else if (xaodEvent->contains<xAOD::TruthParticleContainer>(
                 "TruthHFWithDecayParticles")) {
    if (!xaodEvent->retrieve(truthBHadrons, "TruthHFWithDecayParticles")
             .isSuccess()) {
      throw std::runtime_error("Could not retrieve truth hadrons with key "
                               "TruthHFWithDecayParticles");
    }
  } else if (xaodEvent->contains<xAOD::TruthParticleContainer>(
                 "TruthBSMWithDecayParticles")) {
    if (!xaodEvent->retrieve(truthBHadrons, "TruthBSMWithDecayParticles")
             .isSuccess()) {
      throw std::runtime_error("Could not retrieve truth hadrons with key "
                               "TruthBSMWithDecayParticles");
    }
  } else {
    xAOD::TruthParticleContainer *emptyContainer =
        new xAOD::TruthParticleContainer;
    xAOD::AuxContainerBase *emptyAux = new xAOD::AuxContainerBase();
    emptyContainer->setStore(emptyAux);
    truthBHadrons = emptyContainer;
    WARN_ONCE("Warning: No TruthParticles or TruthBSM container in input - "
              "can't find b-hadrons");
  }

  for (xAOD::TruthParticleContainer::const_iterator it = truthBHadrons->begin();
       it != truthBHadrons->end(); ++it) {
    const auto part = *it;
    if (!(part->isBottomHadron()))
      continue;
    if (!part->hasDecayVtx())
      continue;

    bool weakDecay = true;
    auto decvtx = part->decayVtx();
    int nChildren = decvtx->nOutgoingParticles();
    for (int i = 0; i < nChildren; i++) {
      const auto *child = part->child(i);
      if (!child)
        continue;
      if (child->isBottomHadron())
        weakDecay = false;
    }
    if (weakDecay) {
      tlv.SetPtEtaPhiM(part->pt() / 1000., part->eta(), part->phi(),
                       part->m() / 1000.);
      event->addBHadron(tlv, 0, idx++);
    }
  }
  idx = 0;
  const xAOD::JetContainer *truthtrackjets = 0;
  std::string trackjetName = "AntiKtVR30Rmax4Rmin02TruthChargedJets";
  if (!xaodEvent->contains<xAOD::JetContainer>(trackjetName))
    trackjetName = "AntiKt2TruthChargedJets";

  if (xaodEvent->contains<xAOD::JetContainer>(trackjetName)) {
    if (!xaodEvent->retrieve(truthtrackjets, trackjetName).isSuccess()) {
      throw std::runtime_error("Could not retrieve truth particles with key " +
                               trackjetName);
    }
    for (xAOD::JetContainer::const_iterator it = truthtrackjets->begin();
         it != truthtrackjets->end(); ++it) {
      const auto jet = *it;
      tlv.SetPtEtaPhiM(jet->pt() / 1000., jet->eta(), jet->phi(),
                       jet->m() / 1000.);
      int flavor = 0;
      if (jet->isAvailable<int>("HadronConeExclTruthLabelID"))
        flavor = abs(jet->auxdata<int>("HadronConeExclTruthLabelID"));
      else {
        if (jet->isAvailable<int>("GhostBHadronsFinalCount")) {
          if (jet->auxdata<int>("GhostBHadronsFinalCount")) {
            flavor = 5;
          } else if (jet->auxdata<int>("GhostCHadronsFinalCount")) {
            flavor = 4;
          } else
            flavor = 1;
        }
      }
      int id = (flavor == 5) ? GoodBJet : GoodJet;
      if (flavor == 4)
        id |= TrueCJet;
      else if (flavor == 5)
        id |= TrueBJet;
      else
        id |= TrueLightJet;
      event->addTrackJet(tlv, id, idx++);
    }
  }

  idx = 0;
  const xAOD::TruthParticleContainer *truthBSM = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthBSM") &&
      _useTruthBSM) {
    if (!xaodEvent->retrieve(truthBSM, "TruthBSM").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthBSM");
    }
  } else {
    truthBSM = findTruthBSMParticles(store, truthparticles);
  }
  for (xAOD::TruthParticleContainer::const_iterator it = truthBSM->begin();
       it != truthBSM->end(); ++it) {
    const auto BSM = *it;
    tlv.SetPtEtaPhiM(BSM->pt() / 1000., BSM->eta(), BSM->phi(),
                     BSM->m() / 1000.);
    int pdgId = BSM->pdgId();
    int status = BSM->status();
    if (abs(pdgId) >= (1 << 24))
      WARN_ONCE("Event has particle with PDG ID above 16777216 - not supported "
                "at the moment");
    if (status < 0 || status > 127)
      WARN_ONCE("Events has particle with status<0 or >127 - not supported at "
                "the moment");
    int id = abs(pdgId) | (status << 24);
    if (pdgId < 0)
      id |= 1 << 31;

    auto obj =
        event->addHSTruth(tlv, BSM->charge(), id, getMotherID(BSM), idx++);
    if ((*it)->hasProdVtx()) {
      obj->setProdVtx((*it)->prodVtx()->v4().Vect());
    }
    if ((*it)->hasDecayVtx()) {
      obj->setDecayVtx((*it)->decayVtx()->v4().Vect());
    }
  }

  const xAOD::TruthParticleContainer *truthneutrinos = 0;
  std::string neutrinoName = "TruthNeutrinos";
  if (xaodEvent->contains<xAOD::TruthParticleContainer>(neutrinoName)) {
    if (!xaodEvent->retrieve(truthneutrinos, neutrinoName).isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthNeutrinos");
    }
  } else {
    truthneutrinos = findTruthParticles(store, truthparticles, {12, 14, 16});
  }

  for (xAOD::TruthParticleContainer::const_iterator it =
           truthneutrinos->begin();
       it != truthneutrinos->end(); ++it) {
    const auto neutrino = *it;
    tlv.SetPtEtaPhiM(neutrino->pt() / 1000., neutrino->eta(), neutrino->phi(),
                     0);
    int pdgId = neutrino->pdgId();
    int status = neutrino->status();
    if (abs(pdgId) >= (1 << 24))
      WARN_ONCE("Event has particle with PDG ID above 16777216 - not supported "
                "at the moment");
    if (status < 0 || status > 127)
      WARN_ONCE("Events has particle with status<0 or >127 - not supported at "
                "the moment");
    int id = abs(pdgId) | (status << 24);
    if (pdgId < 0)
      id |= 1 << 31;
    event->addHSTruth(tlv, 0, id, getMotherID(neutrino), idx++);
  }

  const xAOD::TruthParticleContainer *truthTop = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthTop")) {
    if (!xaodEvent->retrieve(truthTop, "TruthTop").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthTop");
    }
  } else {
    truthTop = findTruthParticles(store, truthparticles, {6});
  }

  for (xAOD::TruthParticleContainer::const_iterator it = truthTop->begin();
       it != truthTop->end(); ++it) {
    const auto Top = *it;
    bool top_isAtTheEndOfRadiativeCorrectionChain = true;
    for (unsigned int iChild = 0; iChild < Top->nChildren(); iChild++) {
      const xAOD::TruthParticle *child = Top->child(iChild);
      if (!child)
        continue;
      if (child->pdgId() == Top->pdgId()) {
        top_isAtTheEndOfRadiativeCorrectionChain = false;
        break;
      }
    }
    if (!top_isAtTheEndOfRadiativeCorrectionChain)
      continue;
    tlv.SetPtEtaPhiM(Top->pt() / 1000., Top->eta(), Top->phi(),
                     Top->m() / 1000.);
    int pdgId = Top->pdgId();
    int status = Top->status();
    if (abs(pdgId) >= (1 << 24))
      WARN_ONCE("Event has particle with PDG ID above 16777216 - not supported "
                "at the moment");
    if (status < 0 || status > 127)
      WARN_ONCE("Events has particle with status<0 or >127 - not supported at "
                "the moment");
    int id = abs(pdgId) | (status << 24);
    if (pdgId < 0)
      id |= 1 << 31;

    event->addHSTruth(tlv, Top->charge(), id, getMotherID(Top, true), idx++);
  }

  const xAOD::TruthParticleContainer *truthBottom = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthBottom")) {
    if (!xaodEvent->retrieve(truthBottom, "TruthBottom").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthBottom");
    }
  } else {
    truthBottom = findTruthParticles(store, truthparticles, {5});
  }

  for (xAOD::TruthParticleContainer::const_iterator it = truthBottom->begin();
       it != truthBottom->end(); ++it) {
    const auto Bottom = *it;
    bool bottom_isAtTheEndOfRadiativeCorrectionChain = true;
    for (unsigned int iChild = 0; iChild < Bottom->nChildren(); iChild++) {
      const xAOD::TruthParticle *child = Bottom->child(iChild);
      if (!child)
        continue;
      if (child->pdgId() == Bottom->pdgId()) {
        bottom_isAtTheEndOfRadiativeCorrectionChain = false;
        break;
      }
    }
    if (!bottom_isAtTheEndOfRadiativeCorrectionChain)
      continue;
    tlv.SetPtEtaPhiM(Bottom->pt() / 1000., Bottom->eta(), Bottom->phi(),
                     Bottom->m() / 1000.);
    int pdgId = Bottom->pdgId();
    int status = Bottom->status();
    if (abs(pdgId) >= (1 << 24))
      WARN_ONCE("Event has particle with PDG ID above 16777216 - not supported "
                "at the moment");
    if (status < 0 || status > 127)
      WARN_ONCE("Events has particle with status<0 or >127 - not supported at "
                "the moment");
    int id = abs(pdgId) | (status << 24);
    if (pdgId < 0)
      id |= 1 << 31;

    event->addHSTruth(tlv, Bottom->charge(), id, getMotherID(Bottom, true),
                      idx++);
  }

  const xAOD::TruthParticleContainer *truthBoson = 0;
  if (xaodEvent->contains<xAOD::TruthParticleContainer>("TruthBoson")) {
    if (!xaodEvent->retrieve(truthBoson, "TruthBoson").isSuccess()) {
      throw std::runtime_error(
          "Could not retrieve truth particles with key TruthBoson");
    }
  } else {
    truthBoson = findTruthParticles(store, truthparticles, {23, 24, 25});
  }

  for (xAOD::TruthParticleContainer::const_iterator it = truthBoson->begin();
       it != truthBoson->end(); ++it) {
    const auto Boson = *it;
    bool boson_isAtTheEndOfRadiativeCorrectionChain = true;
    for (unsigned int iChild = 0; iChild < Boson->nChildren(); iChild++) {
      const xAOD::TruthParticle *child = Boson->child(iChild);
      if (!child)
        continue;
      if (child->pdgId() == Boson->pdgId()) {
        boson_isAtTheEndOfRadiativeCorrectionChain = false;
        break;
      }
    }
    if (!boson_isAtTheEndOfRadiativeCorrectionChain)
      continue;
    tlv.SetPtEtaPhiM(Boson->pt() / 1000., Boson->eta(), Boson->phi(),
                     Boson->m() / 1000.);
    int pdgId = Boson->pdgId();
    int status = Boson->status();
    if (abs(pdgId) >= (1 << 24))
      WARN_ONCE("Event has particle with PDG ID above 16777216 - not supported "
                "at the moment");
    if (status < 0 || status > 127)
      WARN_ONCE("Events has particle with status<0 or >127 - not supported at "
                "the moment");
    int id = abs(pdgId) | (status << 24);
    if (pdgId < 0)
      id |= 1 << 31;

    event->addHSTruth(tlv, Boson->charge(), id, getMotherID(Boson), idx++);
  }

  // Special case for ttZ: in some cases (aMC@NLO, Sherpa) off-shell Z bosons
  // don't show up in the truth record. Instead we look for a pair of
  // fermion-anti-fermion attached to the ttbar vertex. Only do this if we
  // haven't picked up any Z in the loop above:
  auto Zcandidates =
      event->getHSTruth(0., 10., 23); // pT and eta cuts should cover all cases
  auto Tcandidates = event->getHSTruth(0., 10., 6);
  if (Zcandidates.size() < 1 &&
      Tcandidates.size() > 1) { // ttZ event should have 2 tops
    const xAOD::TruthParticleContainer *truthAll =
        0; // get our own, in case truthparticles is from BSM collection
    if (!truthAll && xaodEvent->contains<xAOD::TruthParticleContainer>(
                         "TruthParticles")) { // in TRUTH1
      if (!xaodEvent->retrieve(truthAll, "TruthParticles").isSuccess()) {
        throw std::runtime_error(
            "Could not retrieve truth particles with key TruthParticles");
      }
    } else {
      WARN_ONCE("No TruthParticles collection available! --> will check "
                "HardScatterParticles");
    }
    if (!truthAll && xaodEvent->contains<xAOD::TruthParticleContainer>(
                         "HardScatterParticles")) { // used to be in TRUTH3
      if (!xaodEvent->retrieve(truthAll, "HardScatterParticles").isSuccess()) {
        throw std::runtime_error(
            "Could not retrieve truth particles with key HardScatterParticles");
      }
    } else {
      WARN_ONCE("No HardScatterParticles collection available --> will build "
                "custom collection");
    }
    if (truthAll) {
      for (const auto &p : *truthAll) {
        if (abs(p->pdgId()) > 19)
          continue; // consider only elementary fermions
        if (p->pdgId() < 0)
          continue; // no anti-particle (to avoid double-counting)
        const auto &sibling = getFlavourSibling(
            p); // find opposite pdgId particle from same parent
        if (!sibling)
          continue;
        bool has_top_sibling{false}, has_antitop_sibling{false};
        const auto &parent =
            p->parent(0); // the parent of our candidate fermion
        if (!parent)
          continue;
        for (size_t i = 0; i < parent->nChildren(); ++i) {
          const auto *child = parent->child(i);
          if (!child)
            continue;
          if (child == p)
            continue; // don't look at our candidate!
          if (child->pdgId() == 6)
            has_top_sibling = true; // self-explanatory
          if (child->pdgId() == -6)
            has_antitop_sibling = true;
          if (has_top_sibling && has_antitop_sibling)
            break; // only possible if our candidate fermion is non-top!
        }
        if (!(has_top_sibling && has_antitop_sibling))
          continue;
        // now we can build a Z candidate from the selected fermion and its
        // flavour sibling
        TLorentzVector ZDecay1 = p->p4();
        TLorentzVector ZDecay2 = sibling->p4();
        TLorentzVector temp = ZDecay1 + ZDecay2;
        TLorentzVector Ztlv;
        Ztlv.SetPtEtaPhiE(temp.Pt() / 1000., temp.Eta(), temp.Phi(),
                          temp.E() / 1000.);
        // set particle ID and write to collection
        int pdgId = 23;
        int status = 22;
        int id = abs(pdgId) | (status << 24);
        if (pdgId < 0)
          id |= 1 << 31;
        event->addHSTruth(Ztlv, 0, id, 0, idx++);
        break;
      }
    }
  }

  // Get LHE3 weights
  const xAOD::TruthEventContainer *truthEvtCont;
  if (!xaodEvent->retrieve(truthEvtCont, "TruthEvents").isSuccess())
    throw std::runtime_error(
        "Could not retrieve truth event container with key TruthEvents");
  const xAOD::TruthEvent *truthevent = (*truthEvtCont)[0];
  std::vector<float> weights = eventInfo->mcEventWeights();
  if (truthevent->isAvailable<float>("Q")) {
    xAOD::TruthEvent::PdfInfo pdfInfo = truthevent->pdfInfo();
    if (pdfInfo.valid()) {
      event->setPDFInfo(pdfInfo.pdgId1, pdfInfo.x1, pdfInfo.xf1, pdfInfo.pdgId2,
                        pdfInfo.x2, pdfInfo.xf2, pdfInfo.Q);
    }
  } else {
    WARN_ONCE("No PDF information available");
    event->setPDFInfo(0, 0, 0, 0, 0, 0, 0);
  }
  event->setMCWeights(weights);

  _analysisRunner->processEvent(event, eventNumber);

  delete event;
  return true;
}

void xAODTruthReader::processFilesInternal(
    const std::vector<std::string> &inputFileNames, unsigned int nevents) {
  xAOD::TStore transientStorage;
  transientStorage.setActive();
  TFile *inFile = 0;
  unsigned int procEvents = 0;
  for (const auto &inName : inputFileNames) {
    delete inFile;
    std::cout << "Now reading: " << inName << std::endl;
    inFile = TFile::Open(inName.c_str());
    if (!_event->readFrom(inFile).isSuccess()) {
      throw std::runtime_error("Could not connect TEvent to file !");
    }
    Long64_t numEntries = _event->getEntries();
    for (Long64_t index = 0; index < numEntries; index++) {
      ++procEvents;
      if (procEvents > nevents)
        break;
      Long64_t entry = _event->getEntry(index);
      if (entry < 0)
        break;
      if (index % 10000 == 0)
        std::cout << "at: " << index << "/" << numEntries << std::endl;
      processEvent(_event, &transientStorage);
      transientStorage.clear();
    }
  }
}

const xAOD::TruthParticle *
xAODTruthReader::getFlavourSibling(const xAOD::TruthParticle *particle) {
  const auto &parent = particle->parent(0);
  if (!parent)
    return nullptr;

  for (size_t i = 0; i < parent->nChildren(); ++i) {
    const auto &sibling_candidate = parent->child(i);
    if (!sibling_candidate)
      continue;
    if (sibling_candidate->pdgId() == -particle->pdgId())
      return sibling_candidate;
  }

  return nullptr;
}

xAODTruthReader::~xAODTruthReader() {
  return;
  delete _mctool;
  delete _susytools;
  delete _event;
}
