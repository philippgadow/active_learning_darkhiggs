/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "MCEventVetoHelper.h"

#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"

bool MCEventVetoHelper::isHighPtDijet(const xAOD::JetContainer *jets) {

  float leadjetPt = -1.;
  float secondjetPt = -1.;
  for (xAOD::JetContainer::const_iterator it = jets->begin(); it != jets->end();
       ++it) {
    double pt = (*it)->pt();
    if (pt > leadjetPt) {
      secondjetPt = leadjetPt;
      leadjetPt = pt;
    } else if (pt > secondjetPt) {
      secondjetPt = pt;
    }
  }
  // std::cout << " leading thruth jet Pts " <<  leadjetPt << " " << secondjetPt
  // << std::endl;
  if (secondjetPt < 100000.)
    return false;

  return true;
}

bool MCEventVetoHelper::isHighPtJetMET(uint32_t mc_channel_number,
                                       const xAOD::JetContainer *jets,
                                       const xAOD::MissingETContainer *metc) {
  if (mc_channel_number < 107680 || mc_channel_number > 107720) {
    return false;
  }

  xAOD::MissingETContainer::const_iterator truthmet_it = metc->find("NonInt");
  if (truthmet_it == metc->end()) {
    throw(std::runtime_error("Could not find Truth MET with name NonInt"));
  }
  // std::cout << " MET Truth_NonInt " <<  (*truthmet_it)->met() << std::endl;
  if ((*truthmet_it)->met() < 100000.)
    return false;

  float leadjetPt = -1.;
  for (xAOD::JetContainer::const_iterator it = jets->begin(); it != jets->end();
       ++it) {
    double pt = (*it)->pt();
    if (pt > leadjetPt) {
      leadjetPt = pt;
    }
  }
  // std::cout << " leading thruth jet Pt " <<  leadjetPt << std::endl;
  if (leadjetPt < 100000.)
    return false;

  return true;
}

unsigned int
MCEventVetoHelper::vetoQEDFSR(uint32_t mc_channel_number,
                              const xAOD::TruthParticleContainer *mcparticles) {
  // Remove overlap between Sherpa Zll+gamma (Wlv+gamma) and inclusive Z (W)
  // samples Veto events from the Z or W inclusive samples if a QED FSR is
  // emitted with a large angle (>0.1) with respect to the charged lepton.
  unsigned int veto = 0;
  bool isSherpaWZ = false;
  isSherpaWZ = isSherpaWZ || (mc_channel_number >= 147774 &&
                              mc_channel_number <= 147776); // Inclusive W
  isSherpaWZ =
      isSherpaWZ || (mc_channel_number >= 144992 &&
                     mc_channel_number <= 144994); // Sherpa W + >= 3jets
  isSherpaWZ =
      isSherpaWZ || (mc_channel_number >= 157534 &&
                     mc_channel_number <= 157536); // Sherpa W + PtW > 200 GeV
  isSherpaWZ = isSherpaWZ || (mc_channel_number >= 147770 &&
                              mc_channel_number <= 147772); // Inclusive Zll
  if (!isSherpaWZ)
    return veto;

  for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
       it != mcparticles->end(); ++it) {
    // std::cout << " MCparticle " << (*it)->pdgId() << " " << (*it)->barcode()
    // << " " <<  (*it)->status() << " " << (*it)->pt() << std::endl;
    // FIXME is the order now such that we can break the loop in this case?
    if ((*it)->barcode() >= 200000)
      break; // GEANT particles, skip brem photons
    if ((*it)->pdgId() != 22)
      continue;
    if ((*it)->pt() < 10000.)
      continue;
    size_t nParents = 0;
    if ((*it)->hasProdVtx())
      nParents = (*it)->prodVtx()->nIncomingParticles();
    // std::cout << " MCparticle " << (*it)->pdgId() << " " << (*it)->barcode()
    // << (*it)->status() << " " << (*it)->pt() << " " << nParents << std::endl;
    if (nParents < 2)
      continue;

    unsigned int nlep = 0;
    float dRmin = 999.;
    const xAOD::TruthVertex *originVertex = (*it)->prodVtx();
    for (size_t i = 0; i < nParents; ++i) {
      const xAOD::TruthParticle *parent = originVertex->incomingParticle(i);
      if (parent->absPdgId() < 11 || parent->absPdgId() > 16)
        continue;
      nlep++;
      float dR = (*it)->p4().DeltaR(parent->p4());
      if (dR < dRmin)
        dRmin = dR;
    }
    if ((nlep < 2) || (dRmin < 0.1))
      continue;
    veto += 1000;
    break;
  }
  return veto;
}

bool MCEventVetoHelper::mc12AlpgenZjets_accept(
    uint32_t channel, const xAOD::TruthParticleContainer *mcparticles) {
  if (channel > 156803 && channel < 156829) // Alpgen Z+jets
  {
    float ptmin = 0.;
    float ptmax = -1.;
    switch (channel) {
    case 156804:
    case 156809:
    case 156814:
    case 156819:
    case 156824:
      ptmin = 0.f;
      ptmax = 80000.f;
      break;
    case 156808:
    case 156813:
    case 156818:
    case 156823:
    case 156828:
      ptmin = 80000.f;
      ptmax = 150000.f;
      break;
    case 156805:
    case 156810:
    case 156815:
    case 156820:
    case 156825:
      ptmin = 150000.f;
      ptmax = 290000.f;
      break;
    case 156806:
    case 156811:
    case 156816:
    case 156821:
    case 156826:
      ptmin = 290000.f;
      ptmax = 510000.f;
      break;
    case 156807:
    case 156812:
    case 156817:
    case 156822:
    case 156827:
      ptmin = 510000.f;
      ptmax = 1.e30;
      break;
    }
    if (ptmax < 0.)
      throw std::logic_error(
          "No Z pt cut defined for Algen Zjets MC12 sample !");

    float zpt = 0.f;
    size_t count = 0;
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      count++;
      if (count >= 20)
        break;

      if ((*it)->status() == 124 && (*it)->pdgId() == 23) {
        if (zpt > 0.)
          throw std::logic_error(
              "More that one Z with status=124 in Algen Z+jet sample !");
        zpt = (*it)->pt();
      }
    }
    if (zpt < ptmin || zpt >= ptmax)
      return false;
  }
  return true;
}

bool MCEventVetoHelper::mc12AlpgenYjets_accept(
    uint32_t channel, const xAOD::TruthParticleContainer *mcparticles) {
  if (channel >= 156839 && channel <= 156863) // Alpgen gamma+jets
  {
    float ptmin = 0.;
    float ptmax = -1.;
    switch (channel) {
    case 156841:
    case 156846:
    case 156851:
    case 156856:
    case 156861:
      ptmin = 45000.f;
      ptmax = 80000.f;
      break;
    case 156843:
    case 156848:
    case 156853:
    case 156858:
    case 156863:
      ptmin = 80000.f;
      ptmax = 150000.f;
      break;
    case 156839:
    case 156844:
    case 156849:
    case 156854:
    case 156859:
      ptmin = 150000.f;
      ptmax = 290000.f;
      break;
    case 156840:
    case 156845:
    case 156850:
    case 156855:
    case 156860:
      ptmin = 290000.f;
      ptmax = 510000.f;
      break;
    case 156842:
    case 156847:
    case 156852:
    case 156857:
    case 156862:
      ptmin = 510000.f;
      ptmax = 1.e30;
      break;
    }
    if (ptmax < 0.)
      throw std::logic_error(
          "No gamma pt cut defined for Algen gamm+jets MC12 sample !");

    float gampt = 0.f;
    size_t count = 0;
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      count++;
      if (count >= 20)
        break;

      if ((*it)->status() == 124 && (*it)->pdgId() == 22) {
        if (gampt > 0.)
          throw std::logic_error("More that one gamma with status=124 in Algen "
                                 "gamma+jet sample !");
        gampt = (*it)->pt();
      }
    }
    if (gampt < ptmin || gampt >= ptmax)
      return false;
  }
  return true;
}

void MCEventVetoHelper::mc12AlpgenJimmyW_accept(
    unsigned int &veto, uint32_t channel,
    const xAOD::TruthParticleContainer *mcparticles,
    const xAOD::MissingETContainer *metc) {
  if ((channel >= 107680 && channel <= 107685) ||
      (channel >= 107690 && channel <= 107695) ||
      (channel >= 172001 && channel <= 172006) ||
      (channel >= 172011 && channel <= 172016)) { // Alpgen+Jimmy We/munu+jets
    bool isbaseline = (channel >= 107680 && channel <= 107685) ||
                      (channel >= 107690 && channel <= 107695);
    bool ishighptnu = false;
    bool ishighptq = false;
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      if ((*it)->status() != 123 && (*it)->status() != 124)
        continue;
      if (((*it)->absPdgId() == 12 || (*it)->absPdgId() == 14) &&
          (*it)->pt() > 105000.)
        ishighptnu = true;
      if (((*it)->absPdgId() <= 6 || (*it)->absPdgId() == 21) &&
          (*it)->pt() > 85000.)
        ishighptq = true;
      if (ishighptnu && ishighptq)
        break;
    }
    if (isbaseline && ishighptnu && ishighptq)
      veto += 4;
    if (!isbaseline && !(ishighptnu && ishighptq))
      veto += 4;
  }

  else if ((channel >= 107700 && channel <= 107705) ||
           (channel >= 172021 &&
            channel <= 172026)) { // Alpgen+Jimmy Wtaunu+jets
    bool isbaseline = (channel >= 107700 && channel <= 107705);
    xAOD::MissingETContainer::const_iterator truthmet_it = metc->find("NonInt");
    if (truthmet_it == metc->end()) {
      throw(std::runtime_error("Could not find Truth MET with name NonInt"));
    }
    float MET_Truth_NonInt_et = (*truthmet_it)->met();

    bool ishighptnu = (MET_Truth_NonInt_et > 100000.);
    bool ishighptq = false;
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      if ((*it)->status() != 123 && (*it)->status() != 124)
        continue;
      if (((*it)->absPdgId() <= 6 || (*it)->absPdgId() == 21) &&
          (*it)->pt() > 85000.) {
        ishighptq = true;
        break;
      }
    }
    if (isbaseline && ishighptnu && ishighptq)
      veto += 4;
    if (!isbaseline && !(ishighptnu && ishighptq))
      veto += 4;
  }
}

void MCEventVetoHelper::mc12SherpaZnunu_accept(
    unsigned int &veto, uint32_t channel,
    const xAOD::TruthParticleContainer *mcparticles) {
  if (channel >= 157537 && channel <= 157540) { // Sherpa Zvv+jets
    float ptmin = 0.;
    float ptmax = -1.;
    switch (channel) {
    case 157537:
      ptmin = 0.f;
      ptmax = 150000.f;
      break;
    case 157538:
      ptmin = 150000.f;
      ptmax = 290000.f;
      break;
    case 157539:
      ptmin = 290000.f;
      ptmax = 510000.f;
      break;
    case 157540:
      ptmin = 510000.f;
      ptmax = 1.e30;
      break;
    }
    if (ptmax < 0.)
      throw std::logic_error(
          "No Z pt cut defined for Sherpa MC12 Znunu sample !");
    TLorentzVector tlvZ(0., 0., 0., 0.);
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      if ((*it)->status() == 3 &&
          ((*it)->absPdgId() == 12 || (*it)->absPdgId() == 14 ||
           (*it)->absPdgId() == 16)) {
        tlvZ += (*it)->p4();
      }
    }
    if (tlvZ.Pt() < ptmin || tlvZ.Pt() >= ptmax)
      veto += 1;
  }
}

void MCEventVetoHelper::mc12SherpaYjets_accept(
    unsigned int &veto, uint32_t channel,
    const xAOD::TruthParticleContainer *mcparticles) {
  if ((channel >= 113714 && channel <= 113717) || (channel == 126371) ||
      (channel >= 126955 && channel <= 126956)) { // Sherpa gamma+jets
    float ptmin = 0.;
    float ptmax = -1.;
    switch (channel) {
    case 113714:
      ptmin = 0.f;
      ptmax = 80000.f;
      break;
    case 113715:
      ptmin = 80000.f;
      ptmax = 150000.f;
      break;
    case 113716:
      ptmin = 150000.f;
      ptmax = 290000.f;
      break;
    case 113717:
      ptmin = 290000.f;
      ptmax = 510000.f;
      break;
    case 126371:
      ptmin = 510000.f;
      ptmax = 810000.f;
      break;
    case 126955:
      ptmin = 810000.f;
      ptmax = 1010000.f;
      break;
    case 126956:
      ptmin = 1010000.f;
      ptmax = 1.e30;
      break;
    }
    if (ptmax < 0.)
      throw std::logic_error(
          "No gamma pt cut defined for Sherpa MC12 sample !");

    float gampt = 0.f;
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      if ((*it)->status() == 3 && (*it)->pdgId() == 22) {
        gampt = (*it)->pt();
        break;
      }
    }
    if (gampt < ptmin || gampt >= ptmax)
      veto += 1;
  }
}

bool MCEventVetoHelper::trueBosonFromWorZplusJetsMCSample(
    TLorentzVector &trueBoson, uint32_t /*mc_channel_number*/,
    const xAOD::TruthParticleContainer *mcparticles) {
  std::vector<TLorentzVector> true_leptons;

  for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
       it != mcparticles->end(); ++it) {
    int id = (*it)->absPdgId();

    //============================================================
    // In Sherpa the true boson 4momentum is computed
    // from the 2 leptons with status==2
    //============================================================
    // FIXME comment says status==2, code status==3
    if (id >= 10 && id <= 19 && (*it)->status() == 3) {
      TLorentzVector lep = TLorentzVector();
      true_leptons.push_back((*it)->p4());
    }
    if (true_leptons.size() >= 2) {
      trueBoson = true_leptons[0] + true_leptons[1];
      return true;
    }

    //============================================================
    // In alpgen the true W or Z is the first boson with status=155
    //============================================================
    if ((id == 23 || id == 24) && (*it)->status() == 155) {
      trueBoson = (*it)->p4();
      return true;
    }
  }
  return true;
}

bool MCEventVetoHelper::mc12SherpaWZjets_accept(
    unsigned int &veto, uint32_t channel,
    const xAOD::TruthParticleContainer *mcparticles) {
  if (channel >= 167150 && channel <= 167155) { // Sherpa W+jets MassiveB
    TLorentzVector vboson;
    bool isbosonfound = MCEventVetoHelper::trueBosonFromWorZplusJetsMCSample(
        vboson, channel, mcparticles);
    if (!isbosonfound)
      return true;
    if (vboson.Pt() > 140000.)
      veto += 200;

  }

  else if ((channel >= 167740 && channel <= 167748) ||
           (channel >= 167749 &&
            channel <= 167760)) { // Sherpa W/Z+jets MassiveCB
    TLorentzVector vboson;
    bool isbosonfound = MCEventVetoHelper::trueBosonFromWorZplusJetsMCSample(
        vboson, channel, mcparticles);
    if (!isbosonfound)
      return true;
    if (vboson.Pt() > 70000.)
      veto += 300;

  }

  else if ((channel >= 144992 && channel <= 144994) ||
           (channel >= 147774 && channel <= 147776)) { // Sherpa W+jets
    TLorentzVector vboson;
    bool isbosonfound = MCEventVetoHelper::trueBosonFromWorZplusJetsMCSample(
        vboson, channel, mcparticles);
    if (!isbosonfound)
      return true;
    if (vboson.Pt() > 210000.)
      veto += 100;
    if (channel >= 144992 && channel <= 144994)
      return true;

    const xAOD::TruthVertex *originVertex = 0;
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      if ((*it)->absPdgId() >= 11 && (*it)->absPdgId() <= 16 &&
          (*it)->status() == 3) {
        (*it)->prodVtx();
        break;
      }
    }
    if (originVertex) {
      if (originVertex->nOutgoingParticles() - 2 >= 3)
        veto += originVertex->nOutgoingParticles() - 2;
    }
  }
  return true;
}

void MCEventVetoHelper::mc12HerwigVVjets_accept(
    unsigned int &veto, uint32_t channel,
    const xAOD::TruthParticleContainer *mcparticles) {
  if (channel == 105985 ||
      channel ==
          105987) { // Herwig WW or WZ with 1 el or 1 mu filter => skip events
                    // with a W->tau (taken into account by 157951 or 157952)
    bool foundtau = false;
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      if ((*it)->absPdgId() == 24 && (*it)->status() == 195) {
        const xAOD::TruthVertex *decayV = 0;
        if ((*it)->hasDecayVtx())
          decayV = (*it)->decayVtx();
        size_t nDaughters = decayV->nOutgoingParticles();
        for (size_t ichild = 0; ichild < nDaughters; ichild++) {
          if (decayV->outgoingParticle(ichild)->absPdgId() != 15)
            continue;
          foundtau = true;
          veto += 15;
          break;
        }
      }
      if (foundtau)
        break;
    }
  }

  else if (channel >= 161995 &&
           channel <= 161997) { // Herwig WW/WZ/ZZ without filter : skip events
                                // with tau and 1 lepton
    bool foundlep = false;
    for (xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
         it != mcparticles->end(); ++it) {
      if ((*it)->status() == 1 && (*it)->barcode() < 200000 &&
          std::abs((*it)->eta()) < 2.8 && (*it)->pt() > 10000 &&
          ((*it)->absPdgId() == 11 || (*it)->absPdgId() == 13)) {
        foundlep = true;
        break;
      }
    }
    bool foundtau = false;
    if (channel == 161995 || channel == 161996) {
      for (xAOD::TruthParticleContainer::const_iterator it =
               mcparticles->begin();
           it != mcparticles->end(); ++it) {
        if ((*it)->absPdgId() == 24 && (*it)->status() == 195) {
          const xAOD::TruthVertex *decayV = 0;
          if ((*it)->hasDecayVtx())
            decayV = (*it)->decayVtx();
          size_t nDaughters = decayV->nOutgoingParticles();
          for (size_t ichild = 0; ichild < nDaughters; ichild++) {
            if (decayV->outgoingParticle(ichild)->absPdgId() != 15)
              continue;
            foundtau = true;
            veto += 15;
            break;
          }
        }
        if (foundtau)
          break;
      }
    }
    if (foundtau || foundlep)
      veto += 1000;
  }
}

bool MCEventVetoHelper::mc12accept(
    unsigned int &veto, uint32_t mc_channel_number,
    const xAOD::TruthParticleContainer *mcparticles,
    const xAOD::MissingETContainer *metc) {
  veto += vetoQEDFSR(mc_channel_number, mcparticles);

  if (!MCEventVetoHelper::mc12AlpgenZjets_accept(mc_channel_number,
                                                 mcparticles))
    return false;

  if (!MCEventVetoHelper::mc12AlpgenYjets_accept(mc_channel_number,
                                                 mcparticles))
    return false;

  MCEventVetoHelper::mc12AlpgenJimmyW_accept(veto, mc_channel_number,
                                             mcparticles, metc);

  MCEventVetoHelper::mc12SherpaZnunu_accept(veto, mc_channel_number,
                                            mcparticles);
  MCEventVetoHelper::mc12SherpaYjets_accept(veto, mc_channel_number,
                                            mcparticles);

  if (!MCEventVetoHelper::mc12SherpaWZjets_accept(veto, mc_channel_number,
                                                  mcparticles))
    return false;

  MCEventVetoHelper::mc12HerwigVVjets_accept(veto, mc_channel_number,
                                             mcparticles);

  return true;
}

bool MCEventVetoHelper::mc14SherpaWZjets_accept(
    unsigned int &veto, uint32_t channel,
    const xAOD::TruthParticleContainer *mcparticles) {
  if ((channel >= 167740 && channel <= 167748) ||
      (channel >= 167749 && channel <= 167760)) { // Sherpa W/Z+jets MassiveCB
    TLorentzVector vboson;
    bool isbosonfound = MCEventVetoHelper::trueBosonFromWorZplusJetsMCSample(
        vboson, channel, mcparticles);
    if (!isbosonfound)
      return true;
    // W+jets has a 40-70 bin in mc14  but not Z+jet
    if (channel >= 167740 && channel <= 167748) {
      if (vboson.Pt() > 40000.)
        veto += 300;
    } else {
      if (vboson.Pt() > 70000.)
        veto += 300;
    }
  }
  return true;
}

bool MCEventVetoHelper::mc14accept(
    unsigned int &veto, uint32_t mc_channel_number,
    const xAOD::TruthParticleContainer *mcparticles,
    const xAOD::MissingETContainer * /*metc*/) {
  if (!MCEventVetoHelper::mc14SherpaWZjets_accept(veto, mc_channel_number,
                                                  mcparticles))
    return false;
  return true;
}

bool MCEventVetoHelper::mc15accept(
    unsigned int &veto, uint32_t mc_channel_number, const float filtMET,
    const float filtHT, const bool isTruth,
    const xAOD::TruthParticleContainer *truthPC) {

  // veto for filtered samples
  unsigned int filtveto = 0;

  // ttbar nominal (PowhegPythia6)
  vetoHTMETfilter(filtveto, mc_channel_number, 410000, 407012,
                  (bool)(filtHT < 600000.), (bool)(filtMET < 200000.));

  // ttbar Herwigpp
  vetoHTMETfilter(filtveto, mc_channel_number, 410004, 407040,
                  (bool)(filtHT < 600000.), (bool)(filtMET < 200000.));

  // ttbar PowhegPythia8
  if (isTruth) {
    vetoHTMETfilter(filtveto, mc_channel_number, 410500, 407044,
                    (bool)(filtHT < 600000.), (bool)(filtMET < 200000.));
  }

  // ttbar aMC@NLOHerwigpp
  if (isTruth) {
    vetoHTMETfilter(filtveto, mc_channel_number, 410003, 407048,
                    (bool)(filtHT < 600000.), (bool)(filtMET < 200000.));
  }

  // ttbar radHi
  vetoHTMETfilter(filtveto, mc_channel_number, 410001, 407032,
                  (bool)(filtHT < 600000.), (bool)(filtMET < 200000.));

  // ttbar radLo
  vetoHTMETfilter(filtveto, mc_channel_number, 410002, 407036,
                  (bool)(filtHT < 600000.), (bool)(filtMET < 200000.));

  // Wt inclusive top
  vetoHTMETfilter(filtveto, mc_channel_number, 410013, 407019,
                  (bool)(filtHT < 500000.), (bool)(filtMET < 200000.));

  // Wt inclusive antitop
  vetoHTMETfilter(filtveto, mc_channel_number, 410014, 407021,
                  (bool)(filtHT < 500000.), (bool)(filtMET < 200000.));

  // bug for  mc15_13TeV.363359.Sherpa_221_NNPDF30NNLO_WpqqWmlv.evgen.EVNT.e5583
  vetoWplvWmqqOverlap(filtveto, mc_channel_number, truthPC);

  veto += filtveto;
  if (false)
    std::cout << "MCEventVetoHelper::mc15accept : filtveto=" << filtveto
              << " filtHT=" << filtHT << " filtMET=" << filtMET << std::endl;

  return true;
}

bool MCEventVetoHelper::mc16accept(
    unsigned int &veto, uint32_t mc_channel_number, const float filtMET,
    const float filtHT, const bool /* isTruth */,
    const xAOD::TruthParticleContainer * /*truthPC*/) {

  // veto for filtered samples
  unsigned int filtveto = 0;

  // ttbar nominal (PowhegPythia8)
  if (mc_channel_number == 410470 && !(filtHT < 600000. && filtMET < 200000.))
    filtveto = 410470;

  if (mc_channel_number == 407345 && !(filtHT < 600000.))
    filtveto = 407345;

  if (mc_channel_number == 407346 && !(filtHT < 600000.))
    filtveto = 407346;

  if (mc_channel_number == 407347 && !(filtHT < 600000.))
    filtveto = 407347;

  // ttbar Mc@NLO + Pythia 8
  if (mc_channel_number == 410464 && !(filtHT < 600000. && filtMET < 200000.))
    filtveto = 410464;

  if (mc_channel_number == 410465 && !(filtHT < 600000. && filtMET < 200000.))
    filtveto = 410465;

  if (mc_channel_number == 407351 && !(filtHT < 600000.))
    filtveto = 407351;

  if (mc_channel_number == 407352 && !(filtHT < 600000.))
    filtveto = 407352;

  if (mc_channel_number == 407353 && !(filtHT < 600000.))
    filtveto = 407353;

  // ttbar Powheg + Herwig 7
  if (mc_channel_number == 410557 && !(filtHT < 600000. && filtMET < 200000.))
    filtveto = 410557;

  if (mc_channel_number == 410558 && !(filtHT < 600000. && filtMET < 200000.))
    filtveto = 410558;

  if (mc_channel_number == 407357 && !(filtHT < 600000.))
    filtveto = 407357;

  if (mc_channel_number == 407358 && !(filtHT < 600000.))
    filtveto = 407358;

  if (mc_channel_number == 407359 && !(filtHT < 600000.))
    filtveto = 407359;
  // Wt inclusive top  # Kenta comment: MET/HT sliced PowhegPythia8 samples
  // (mc16a, mc16d) are not ready. Therefore, I comment out these two lines.
  // vetoHTMETfilter(filtveto, mc_channel_number, 410646, 407019,
  // (bool)(filtHT<500000.) , (bool)(filtMET<200000.) );
  // Wt inclusive antitop
  // vetoHTMETfilter(filtveto, mc_channel_number, 410647, 407021,
  // (bool)(filtHT<500000.) , (bool)(filtMET<200000.) );

  veto += filtveto;
  if (false)
    std::cout << "MCEventVetoHelper::mc16accept : filtveto=" << filtveto
              << " filtHT=" << filtHT << " filtMET=" << filtMET << std::endl;

  return true;
}

void MCEventVetoHelper::vetoHTMETfilter(unsigned int &filtveto,
                                        uint32_t mc_channel_number,
                                        unsigned int mc_channel_nominal,
                                        unsigned int mc_channel_METfilter,
                                        bool passLowHTcut, bool passLowMETcut) {
  // veto nominal sample if HT>HTcut or MET>METcut
  if (mc_channel_number == mc_channel_nominal &&
      !(passLowMETcut && passLowHTcut))
    filtveto = mc_channel_nominal;
  // veto MET-filtered sample if HT>HTcut
  if (mc_channel_number == mc_channel_METfilter && !passLowHTcut)
    filtveto = mc_channel_METfilter;
}

void MCEventVetoHelper::vetoWplvWmqqOverlap(
    unsigned int &filtveto, uint32_t mc_channel_number,
    const xAOD::TruthParticleContainer *truthParticles) {

  // Only for the WplvWmqq
  if (mc_channel_number != 363360)
    return;

  int nWs = 0;
  for (auto truthParticle : *truthParticles) {
    if (truthParticle->isW()) {
      nWs++;
      if (truthParticle->nChildren() > 0) {
        if (truthParticle->charge() > 0 &&
            !truthParticle->child(0)->isLepton()) {
          filtveto = mc_channel_number;
        } else if (truthParticle->charge() < 0 &&
                   !truthParticle->child(0)->isQuark()) {
          filtveto = mc_channel_number;
        }
      } else {
        std::cout << "ERROR::Can't find a decay for the W-boson in the "
                     "WplvWmqq sample"
                  << std::endl;
        abort();
      }
    }
  }
  if (nWs != 2) {
    std::cout << "ERROR::Found " << nWs
              << "W Bosons instead of 2! in WplvWmqq sample" << std::endl;
    abort();
  }
  return;
}

// reimplement Run1 functions:
// SUSYTUtils::mc12accept(unsigned int& veto)
// SUSYObjDef::Sherpa_WW_veto
// SUSYTUtils::gammaSampleFilter()  NB: 7tev only
