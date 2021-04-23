#include "SimpleAnalysisFramework/AnalysisClass.h"
#include "MT2.h"
#include "TMatrixDSym.h"
#include "TMctLib.h"
#include "TVectorD.h"
#include "TopnessTool.h"

#include <algorithm>

#define FASTJET
#ifdef FASTJET
// Jets
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#endif

#include "PathResolver/PathResolver.h"

const float Z_Mass = 91.2;

std::vector<AnalysisClass *> *getAnalysisList() {
  static std::vector<AnalysisClass *> *list = new std::vector<AnalysisClass *>;
  return list;
}

MVA *AnalysisClass::addTMVA(const std::string &name,
                            const std::vector<std::string> &variableDefs,
                            const std::string fname1,
                            const std::string fname2) {
  if (m_MVAs.find(name) != m_MVAs.end())
    throw std::runtime_error("Duplicate MVA name");
  m_MVAs[name] = new TMVAReader(name, variableDefs, fname1, fname2);
  return m_MVAs[name];
}

MVA *AnalysisClass::addMVAUtilsBDT(const std::string &name,
                                   const std::string fname1,
                                   const std::string fname2) {
  if (m_MVAs.find(name) != m_MVAs.end())
    throw std::runtime_error("Duplicate MVA name");
  m_MVAs[name] = new MVAUtilsReader(name, fname1, fname2);
  return m_MVAs[name];
}

MVA *AnalysisClass::addONNX(const std::string &name, const std::string fname1,
                            const std::string fname2) {
  if (m_MVAs.find(name) != m_MVAs.end())
    throw std::runtime_error("Duplicate MVA name");
  m_MVAs[name] = new ONNXReader(name, fname1, fname2);
  return m_MVAs[name];
}

void AnalysisClass::ntupVar(const std::string &label, AnalysisObject &object,
                            bool saveMass, bool saveType, bool saveVtx,
                            bool saveEnergy) {
  _output->ntupVar(label + "_pt", object.Pt());
  _output->ntupVar(label + "_eta", object.Eta());
  _output->ntupVar(label + "_phi", object.Phi());
  if (saveEnergy)
    _output->ntupVar(label + "_e", object.E());
  _output->ntupVar(label + "_charge", object.charge());
  if (saveMass)
    _output->ntupVar(label + "_m", object.M());
  if (saveType)
    _output->ntupVar(label + "_type", object.type());
  _output->ntupVar(label + "_id", object.id());
  _output->ntupVar(label + "_motherID", object.motherID());
  if (saveVtx) {
    _output->ntupVar(label + "_prodVx", object.prodVtx().X());
    _output->ntupVar(label + "_prodVy", object.prodVtx().Y());
    _output->ntupVar(label + "_prodVz", object.prodVtx().Z());
    _output->ntupVar(label + "_decayVx", object.decayVtx().X());
    _output->ntupVar(label + "_decayVy", object.decayVtx().Y());
    _output->ntupVar(label + "_decayVz", object.decayVtx().Z());
  }
}

void AnalysisClass::ntupVar(const std::string &label, AnalysisObjects &objects,
                            bool saveMass, bool saveType, bool saveVtx,
                            bool saveEnergy) {
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<float> e;
  std::vector<int> charge;
  std::vector<int> type;
  std::vector<int> id;
  std::vector<int> motherid;
  std::vector<float> mass;
  std::vector<float> prodVx;
  std::vector<float> prodVy;
  std::vector<float> prodVz;
  std::vector<float> decayVx;
  std::vector<float> decayVy;
  std::vector<float> decayVz;
  for (const auto &object : objects) {
    pt.push_back(object.Pt());
    eta.push_back(object.Eta());
    phi.push_back(object.Phi());
    e.push_back(object.E());
    charge.push_back(object.charge());
    mass.push_back(object.M());
    type.push_back(object.type());
    id.push_back(object.id());
    motherid.push_back(object.motherID());
    prodVx.push_back(object.prodVtx().X());
    prodVy.push_back(object.prodVtx().Y());
    prodVz.push_back(object.prodVtx().Z());
    decayVx.push_back(object.decayVtx().X());
    decayVy.push_back(object.decayVtx().Y());
    decayVz.push_back(object.decayVtx().Z());
  }
  _output->ntupVar(label + "_pt", pt);
  _output->ntupVar(label + "_eta", eta);
  _output->ntupVar(label + "_phi", phi);
  if (saveEnergy)
    _output->ntupVar(label + "_e", e);
  _output->ntupVar(label + "_charge", charge);
  if (saveMass)
    _output->ntupVar(label + "_m", mass);
  if (saveType)
    _output->ntupVar(label + "_type", type);
  if (saveVtx) {
    _output->ntupVar(label + "_prodVx", prodVx);
    _output->ntupVar(label + "_prodVy", prodVy);
    _output->ntupVar(label + "_prodVz", prodVz);
    _output->ntupVar(label + "_decayVx", decayVx);
    _output->ntupVar(label + "_decayVy", decayVy);
    _output->ntupVar(label + "_decayVz", decayVz);
  }
  _output->ntupVar(label + "_id", id);
  _output->ntupVar(label + "_motherID", motherid);
}

std::vector<float> AnalysisClass::fakeJER(const AnalysisObjects &jets) {
  std::vector<float> jer;

  // Copied from
  // AnalysisBase/2.4.29/JetResolution/Root/JERTool.cxx::getRelResolutionMC
  const int m_nEtaBins = 7;
  double noise[m_nEtaBins] = {};
  double stochastic[m_nEtaBins] = {};
  double constant[m_nEtaBins] = {};

  noise[0] = 3.34;
  stochastic[0] = 0.627;
  constant[0] = 0.0234;
  noise[1] = 3.05;
  stochastic[1] = 0.693;
  constant[1] = 0.0224;
  noise[2] = 3.29;
  stochastic[2] = 0.658;
  constant[2] = 0.0300;
  noise[3] = 2.56;
  stochastic[3] = 0.607;
  constant[3] = 0.0250;
  noise[4] = 0.988;
  stochastic[4] = 0.753;
  constant[4] = 0.0228;
  noise[5] = 2.74;
  stochastic[5] = 0.783;
  constant[5] = 0.0465;
  noise[6] = 2.80;
  stochastic[6] = 0.623;
  constant[6] = 0.0000;
  double etaBins[m_nEtaBins + 1] = {0, 0.8, 1.2, 2.1, 2.8, 3.2, 3.6, 4.5};
  TAxis m_etaAxis = TAxis(m_nEtaBins, etaBins);

  for (const auto &jet : jets) {
    double pt = std::min(std::max(jet.Pt(), 10.), 1500.);
    int etaBin = m_etaAxis.FindBin(fabs(jet.Eta()));
    etaBin = std::min((int)m_nEtaBins, etaBin) - 1;

    float jerMC =
        pow(pow(noise[etaBin] / (pt), 2) + pow(stochastic[etaBin], 2) / (pt) +
                pow(constant[etaBin], 2),
            0.5);
    jer.push_back(jerMC);
  }

  return jer;
}

AnalysisObjects AnalysisClass::overlapRemoval(
    const AnalysisObjects &cands, const AnalysisObjects &others,
    std::function<float(const AnalysisObject &, const AnalysisObject &)>
        radiusFunc,
    int passId) {
  AnalysisObjects reducedList;

  for (const auto &cand : cands) {
    bool overlap = false;
    for (const auto &other : others) {
      if (cand.DeltaR(other) < radiusFunc(cand, other) && cand != other &&
          cand.pass(passId)) {
        overlap = true;
        break;
      }
    }
    if (!overlap)
      reducedList.push_back(cand);
  }
  return reducedList;
}
AnalysisObjects AnalysisClass::overlapRemoval(const AnalysisObjects &cands,
                                              const AnalysisObjects &others,
                                              float deltaR, int passId) {
  return overlapRemoval(
      cands, others,
      [deltaR](const AnalysisObject &, const AnalysisObject &) {
        return deltaR;
      },
      passId);
}

AnalysisObjects AnalysisClass::lowMassRemoval(
    const AnalysisObjects &cand,
    std::function<bool(const AnalysisObject &, const AnalysisObject &)>
        RemoveSwitch,
    float MinMass, float MaxMass, int type) {
  AnalysisObjects reducedList;
  unsigned int NObj = cand.size();
  std::vector<bool> Good(NObj, true);

  AnalysisObjects::const_iterator begin_cand = cand.begin(),
                                  end_cand = cand.end();
  std::vector<bool>::iterator good_begin = Good.begin();

  std::vector<bool>::iterator good = good_begin;
  for (AnalysisObjects::const_iterator obj = begin_cand; obj != end_cand;
       ++obj) {

    // Check if the object is already rejected
    // if (!Good.at(obj)) continue;
    const AnalysisObject &Object = (*obj);
    std::vector<bool>::iterator good1 = good_begin;
    for (AnalysisObjects::const_iterator obj1 = begin_cand; obj1 != obj;
         ++obj1) {
      // if (!Good.at(obj1)) continue;
      // Remove the pair
      const AnalysisObject &Object1 = (*obj1);
      float InvMass = (Object + Object1).M();
      if (MinMass < InvMass && InvMass < MaxMass &&
          RemoveSwitch(Object, Object1)) {
        (*good) = false;
        (*good1) = false;
      }
      ++good1;
    }
    ++good;
  }
  good = good_begin;

  for (AnalysisObjects::const_iterator obj = begin_cand; obj != end_cand;
       ++obj) {
    bool Continue = !(*good);
    ++good;
    if (Continue)
      continue;
    if (type == -1 || obj->type() == type)
      reducedList.push_back(*obj);
  }
  return reducedList;
}
AnalysisObjects AnalysisClass::filterObjects(const AnalysisObjects &cands,
                                             float ptCut, float etaCut, int id,
                                             unsigned int maxNum) {
  AnalysisObjects reducedList;
  unsigned int cnt = 0;
  for (const auto &cand : cands) {
    if ((cand.Pt() >= ptCut) && (fabs(cand.Eta()) < etaCut) && (cand.pass(id)))
      reducedList.push_back(cand);
    cnt++;
    if (cnt == maxNum)
      break;
  }
  return reducedList;
}

AnalysisObjects AnalysisClass::filterObjectsRange(
    const AnalysisObjects &cands, float ptMinCut, float ptMaxCut,
    float etaMinCut, float etaMaxCut, int id, unsigned int maxNum) {
  AnalysisObjects reducedList;
  unsigned int cnt = 0;
  for (const auto &cand : cands) {
    if ((cand.Pt() >= ptMinCut) && (cand.Pt() < ptMaxCut) &&
        (fabs(cand.Eta()) >= etaMinCut) && (fabs(cand.Eta()) < etaMaxCut) &&
        (cand.pass(id)))
      reducedList.push_back(cand);
    cnt++;
    if (cnt == maxNum)
      break;
  }
  return reducedList;
}

AnalysisObjects AnalysisClass::filterCrack(const AnalysisObjects &cands,
                                           float minEta, float maxEta) {
  AnalysisObjects reducedList;
  for (const auto &cand : cands) {
    if (fabs(cand.Eta()) < minEta || fabs(cand.Eta()) > maxEta)
      reducedList.push_back(cand);
  }
  return reducedList;
}

int AnalysisClass::countObjects(const AnalysisObjects &cands, float ptCut,
                                float etaCut, int id) {
  int count = 0;
  for (const auto &cand : cands) {
    if (cand.Pt() >= ptCut && fabs(cand.Eta()) < etaCut && cand.pass(id))
      count++;
  }
  return count;
}

float AnalysisClass::sumObjectsPt(const AnalysisObjects &cands,
                                  unsigned int maxNum, float ptCut) {
  float sum = 0;
  for (int ii = 0; ii < std::min((int)maxNum, (int)cands.size()); ii++)
    if (cands[ii].Pt() > ptCut)
      sum += cands[ii].Pt();
  return sum;
}

float AnalysisClass::sumObjectsM(const AnalysisObjects &cands,
                                 unsigned int maxNum, float mCut) {
  float sum = 0;
  for (int ii = 0; ii < std::min((int)maxNum, (int)cands.size()); ii++)
    if (cands[ii].M() > mCut)
      sum += cands[ii].M();
  return sum;
}

struct pt_sort {
  bool operator()(const TLorentzVector &v1, const TLorentzVector &v2) const {
    return v1.Pt() > v2.Pt();
  }
};
void AnalysisClass::sortObjectsByPt(AnalysisObjects &cands) {
  std::sort(cands.begin(), cands.end(), pt_sort());
}

float AnalysisClass::calcMCT(const AnalysisObject &o1,
                             const AnalysisObject &o2) {
  float mCT = pow(o1.Et() + o2.Et(), 2) - pow(o1.Px() - o2.Px(), 2) -
              pow(o1.Py() - o2.Py(), 2);
  mCT = (mCT >= 0.) ? sqrt(mCT) : sqrt(-mCT);
  return mCT;
}

float AnalysisClass::calcMCT(const AnalysisObject &o1, const AnalysisObject &o2,
                             const AnalysisObject &met, double ecm) {
  TVector2 m;
  m.TVector2::SetMagPhi(met.Et(), met.Phi());
  TLorentzVector *vds = new TLorentzVector();
  vds->SetPtEtaPhiM(0, 0, 0, 0);

  TMctLib mctcalc;
  float mCT = mctcalc.mctcorr(o1, o2, (*vds), m, ecm, 0);

  return mCT;
}

float AnalysisClass::calcMT(const AnalysisObject &lepton,
                            const AnalysisObject &met) {
  float mT =
      2 * lepton.Pt() * met.Et() * (1 - TMath::Cos(lepton.Phi() - met.Phi()));
  mT = (mT >= 0.) ? sqrt(mT) : sqrt(-mT);
  return mT;
}

float AnalysisClass::calcMTmin(const AnalysisObjects &cands,
                               const AnalysisObject &met, int maxNum) {
  float mtmin = 1e9;

  for (int ii = 0; ii < std::min((int)maxNum, (int)cands.size()); ii++)
    mtmin = std::min(mtmin, calcMT(cands[ii], met));
  return mtmin;
}

float AnalysisClass::calcMT2(const AnalysisObject &o1, const AnalysisObject &o2,
                             const AnalysisObject &met) {
  return calcAMT2(o1, o2, met, 0, 0);
}

float AnalysisClass::calcAMT2(const AnalysisObject &o1,
                              const AnalysisObject &o2,
                              const AnalysisObject &met, float m1, float m2) {
  asymm_mt2_lester_bisect::disableCopyrightMessage();
  return asymm_mt2_lester_bisect::get_mT2(o1.M(), o1.Px(), o1.Py(), o2.M(),
                                          o2.Px(), o2.Py(), met.Px(), met.Py(),
                                          m1, m2);
}

float AnalysisClass::calcMTauTau(const AnalysisObject &o1,
                                 const AnalysisObject &o2,
                                 const AnalysisObject &met) {
  float determinant = o1.Px() * o2.Py() - o1.Py() * o2.Px();
  float xi_1 = (met.Px() * o2.Py() - o2.Px() * met.Py()) / determinant;
  float xi_2 = (met.Py() * o1.Px() - o1.Py() * met.Px()) / determinant;

  float MSqTauTau = (1. + xi_1) * (1. + xi_2) * 2 * o1.Dot(o2);

  float MTauTau = 0.;
  if (MSqTauTau >= 0)
    MTauTau = sqrt(MSqTauTau);
  if (MSqTauTau < 0)
    MTauTau = -sqrt(fabs(MSqTauTau));

  return MTauTau;
}

float AnalysisClass::calcTopness(const AnalysisObject &lepton,
                                 const AnalysisObject &met,
                                 const AnalysisObject &jet1,
                                 const AnalysisObject &jet2) {
  TopnessTool mycalc(lepton, TVector2(met.Px(), met.Py()), jet1, jet2);
  return mycalc.minimize();
}

float AnalysisClass::calcTopness(const AnalysisObject &lepton,
                                 const AnalysisObject &met,
                                 const AnalysisObject &jet1,
                                 const AnalysisObject &jet2,
                                 const AnalysisObject &jet3) {
  TopnessTool mycalc(lepton, TVector2(met.Px(), met.Py()), jet1, jet2, jet3);
  return mycalc.minimize();
}

float AnalysisClass::minDphi(const AnalysisObject &met,
                             const AnalysisObjects &cands, unsigned int maxNum,
                             float ptCut) {
  float dphi_min = 999.;
  for (int ii = 0; ii < std::min((int)maxNum, (int)cands.size()); ii++) {
    float dphi = fabs(met.DeltaPhi(cands[ii]));
    if (dphi < dphi_min && cands[ii].Pt() > ptCut)
      dphi_min = dphi;
  }
  return dphi_min;
}

float AnalysisClass::minDphi(const AnalysisObject &met,
                             const AnalysisObject &cand) {
  return fabs(met.DeltaPhi(cand));
}

float AnalysisClass::minDR(const AnalysisObject &cand,
                           const AnalysisObjects &cands, unsigned int maxNum,
                           float ptCut) {
  float dr_min = 999.;
  for (int ii = 0; ii < std::min((int)maxNum, (int)cands.size()); ii++) {
    float dr = fabs(cands[ii].DeltaR(cand));
    if (dr < dr_min && cands[ii].Pt() > ptCut)
      dr_min = dr;
  }
  return dr_min;
}

float AnalysisClass::minDR(const AnalysisObjects &cands, unsigned int maxNum,
                           float ptCut) {
  float dr_min = 999.;
  for (int ii = 0; ii < std::min((int)maxNum, (int)cands.size()); ii++) {
    for (int ij = ii + 1; ij < std::min((int)maxNum, (int)cands.size()); ij++) {
      float dr = fabs(cands[ii].DeltaR(cands[ij]));
      if (dr < dr_min && cands[ii].Pt() > ptCut)
        dr_min = dr;
    }
  }
  return dr_min;
}

AnalysisObjects AnalysisClass::reclusterJets(const AnalysisObjects &jets,
                                             float radius, float ptmin,
                                             float rclus, float ptfrac) {
  std::vector<fastjet::PseudoJet> JetVector;
  for (const auto &jet : jets) {
    fastjet::PseudoJet Pjet(jet.Px(), jet.Py(), jet.Pz(), jet.E());
    JetVector.push_back(Pjet);
  }
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, radius,
                                 recomb_scheme, strategy);

  std::vector<fastjet::PseudoJet> fat_jets;
  fastjet::ClusterSequence clust_seq(JetVector, jet_def);
  if (!JetVector.empty()) {
    fat_jets = clust_seq.inclusive_jets(ptmin);
  }
  AnalysisObjects fatJets;
  fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, rclus),
                          fastjet::SelectorPtFractionMin(ptfrac));
  for (const auto &fat_jet : fat_jets) {
    if (ptfrac >= 0) {
      fastjet::PseudoJet trimmed_jet = trimmer(fat_jet);
      if (trimmed_jet.pt() > ptmin)
        fatJets.push_back(AnalysisObject(trimmed_jet.px(), trimmed_jet.py(),
                                         trimmed_jet.pz(), trimmed_jet.E(), 0,
                                         0, COMBINED, 0, 0));
    } else {
      fatJets.push_back(AnalysisObject(fat_jet.px(), fat_jet.py(), fat_jet.pz(),
                                       fat_jet.E(), 0, 0, COMBINED, 0, 0));
    }
  }
  sortObjectsByPt(fatJets);
  return fatJets;
}

static TVectorD calcEigenValues(const AnalysisObjects &jets) {
  TMatrixDSym momTensor(3);
  TVectorD eigenValues(3);
  momTensor.Zero();
  double norm = 0;
  for (const auto &jet : jets) {
    momTensor(0, 0) += jet.Px() * jet.Px();
    momTensor(0, 1) += jet.Px() * jet.Py();
    momTensor(0, 2) += jet.Px() * jet.Pz();
    momTensor(1, 0) += jet.Py() * jet.Px();
    momTensor(1, 1) += jet.Py() * jet.Py();
    momTensor(1, 2) += jet.Py() * jet.Pz();
    momTensor(2, 0) += jet.Pz() * jet.Px();
    momTensor(2, 1) += jet.Pz() * jet.Py();
    momTensor(2, 2) += jet.Pz() * jet.Pz();
    norm += jet.Vect().Mag2();
  }
  momTensor *= 1. / norm;
  momTensor.EigenVectors(eigenValues);
  return eigenValues;
}

float AnalysisClass::aplanarity(const AnalysisObjects &jets) {
  if (jets.size() < 2)
    return 0;
  TVectorD eigenValues = calcEigenValues(jets);
  return 1.5 * eigenValues(2);
}

float AnalysisClass::sphericity(const AnalysisObjects &jets) {
  if (jets.size() < 2)
    return 0;
  TVectorD eigenValues = calcEigenValues(jets);
  return 1.5 * (eigenValues(1) + eigenValues(2));
}

bool AnalysisClass::IsSFOS(const AnalysisObject &L, const AnalysisObject &L1) {
  return L1.type() == L.type() && L1.charge() * L.charge() < 0.;
}

std::pair<float, float>
AnalysisClass::DiZSelection(const AnalysisObjects &electrons,
                            const AnalysisObjects &muons) {
  std::pair<float, float> ZMasses(-1, -1);
  float BestDZ = 1.e25;
  auto Leptons = electrons + muons;
  AnalysisObjects::const_iterator Lep_begin = Leptons.begin();
  AnalysisObjects::const_iterator Lep_end = Leptons.end();
  unsigned int NumLep = Leptons.size();
  if (NumLep < 2)
    return ZMasses;
  for (AnalysisObjects::const_iterator L = Lep_begin + 1; L != Lep_end; ++L) {
    for (AnalysisObjects::const_iterator L1 = Lep_begin; L1 != L; ++L1) {
      // Loop over the first two leptons
      const AnalysisObject &FirstLep = (*L);
      const AnalysisObject &SecondLep = (*L1);
      // No SFOS pair
      if (!IsSFOS(FirstLep, SecondLep))
        continue;
      float M1 = (FirstLep + SecondLep).M();
      float dM1 = fabs(M1 - Z_Mass);
      // Find the second pair in the row
      for (AnalysisObjects::const_iterator L2 = Lep_begin + 1; L2 != Lep_end;
           ++L2) {
        // Do not use the same object twice
        if (L2 == L1 || L2 == L)
          continue;
        for (AnalysisObjects::const_iterator L3 = Lep_begin; L3 != L2; ++L3) {
          if (L3 == L1 || L3 == L)
            continue;
          const AnalysisObject &ThirdLep = (*L2);
          const AnalysisObject &FourthLep = (*L3);
          if (!IsSFOS(ThirdLep, FourthLep))
            continue;
          float M2 = (ThirdLep + FourthLep).M();
          float dM2 = fabs(M2 - Z_Mass);
          float dZTest = dM1 + dM2;
          // The sum of the distances to M_Z is smaller than anything known
          if (dZTest < BestDZ) {
            // Order them by distance closest -> leading
            if (dM1 < dM2)
              ZMasses = std::pair<float, float>(M1, M2);
            else
              ZMasses = std::pair<float, float>(M2, M1);
            BestDZ = dZTest;
          }
        } // End of the second pair finder

        // No way to get two Zs of the event. Find the best one nervertheless
        if (NumLep < 4) {
          if (dM1 < BestDZ) {
            BestDZ = dM1;
            ZMasses.first = M1;
          }
        }
      }
    }
  }
  return ZMasses;
}

bool AnalysisClass::PassZVeto(const AnalysisObjects &electrons,
                              const AnalysisObjects &muons, float Window) {
  auto Leptons = electrons + muons;
  AnalysisObjects::const_iterator Lep_begin = Leptons.begin();
  AnalysisObjects::const_iterator Lep_end = Leptons.end();
  if (Leptons.size() < 2)
    return true;
  for (AnalysisObjects::const_iterator L = Lep_begin + 1; L != Lep_end; ++L) {
    for (AnalysisObjects::const_iterator L1 = Lep_begin; L1 != L; ++L1) {
      // Loop over the first two leptons
      const AnalysisObject &FirstLep = (*L);
      const AnalysisObject &SecondLep = (*L1);
      AnalysisObject DiLep = (FirstLep + SecondLep);
      if (IsSFOS(FirstLep, SecondLep) && fabs(DiLep.M() - Z_Mass) < Window)
        return false;
      for (AnalysisObjects::const_iterator L2 = Lep_begin; L2 != L1; ++L2) {
        const AnalysisObject &ThirdLep = (*L2);
        AnalysisObject TriLep = DiLep + ThirdLep;
        // Radiated photons might be converted into leptons
        bool TriLepSFOS = IsSFOS(FirstLep, SecondLep) ||
                          IsSFOS(FirstLep, ThirdLep) ||
                          IsSFOS(SecondLep, ThirdLep);
        if (TriLepSFOS && fabs(ThirdLep.M() - Z_Mass) < Window)
          return false;
        // Check forZ ->llll
        for (AnalysisObjects::const_iterator L3 = Lep_begin; L3 != L2; ++L3) {
          const AnalysisObject &FourthLep = (*L3);
          bool FourLepSFOS = (IsSFOS(FirstLep, SecondLep) &&
                              IsSFOS(ThirdLep, FourthLep)) // eemumu
                             || (IsSFOS(FirstLep, ThirdLep) &&
                                 IsSFOS(SecondLep, FourthLep)) // e mu e mu
                             || (IsSFOS(FirstLep, FourthLep) &&
                                 IsSFOS(ThirdLep, SecondLep)); // e mu mu e
          if (FourLepSFOS && fabs((TriLep + FourthLep).M() - Z_Mass) < Window)
            return false;
        }
      }
    }
  }
  return true;
}

struct ClusteringHistory : public fastjet::PseudoJet::UserInfoBase {
  enum Status {
    GOOD,
    JET_TOO_SMALL,
    JET_TOO_LARGE,
    TOO_MANY_ITERATIONS,
    NONE,
  };

  struct Step {
    double pt;
    double r;
    size_t constit;
    Status status;
  };

  size_t id; // a per-event unique jet id that is needed for the event dump
  std::vector<Step> steps;

  static ClusteringHistory *AddStep(ClusteringHistory &history,
                                    const Step &step) {
    auto newHistory = new ClusteringHistory(history);
    newHistory->steps.push_back(step);
    return newHistory;
  }
};

// Return the history of a PseudoJet object, handling all the ugly casting.
ClusteringHistory &GetHistory(const fastjet::PseudoJet &jet) {
  auto shared_ptr = jet.user_info_shared_ptr();
  return *dynamic_cast<ClusteringHistory *>(shared_ptr.get());
}

inline double optimalRadius(const double pT, const double m) {
  return 2 * m / pT;
}
inline double minRadius(const double pT, const double m) {
  return optimalRadius(pT, m) - 0.3;
}
inline double maxRadius(const double pT, const double m) {
  return optimalRadius(pT, m) + 0.5;
}

static std::vector<fastjet::PseudoJet>
SortedByNConstit(std::vector<fastjet::PseudoJet> jets) {
  std::sort(jets.begin(), jets.end(),
            [](const fastjet::PseudoJet &a, const fastjet::PseudoJet &b) {
              if (a.constituents().size() != b.constituents().size())
                return a.constituents().size() > b.constituents().size();

              return a.pt() > b.pt();
            });

  return jets;
}

std::pair<bool, fastjet::PseudoJet>
RecursiveRecluster(const fastjet::PseudoJet &candidate, double candRadius,
                   const double mass, size_t step) {
  if (minRadius(candidate.pt(), mass) > candRadius) {
    GetHistory(candidate).steps.back().status =
        ClusteringHistory::JET_TOO_SMALL;
    return std::make_pair(false, candidate);
  } else if (maxRadius(candidate.pt(), mass) < candRadius) {
    const double newR =
        std::max(maxRadius(candidate.pt(), mass), candRadius / 2.);
    GetHistory(candidate).steps.back().status =
        ClusteringHistory::JET_TOO_LARGE;

    if (step > 10) {
      GetHistory(candidate).steps.back().status =
          ClusteringHistory::TOO_MANY_ITERATIONS;
      return std::make_pair(false, candidate);
    }

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, newR);
    auto cs = new fastjet::ClusterSequence(candidate.constituents(), jetDef);

    std::vector<fastjet::PseudoJet> reclusteredJets;
    reclusteredJets = SortedByNConstit(cs->inclusive_jets());

    if (reclusteredJets.size() == 0) {
      delete cs;
      return std::make_pair(false, fastjet::PseudoJet());
    }

    cs->delete_self_when_unused();
    auto newCandidate = reclusteredJets[0];

    auto newHistory = ClusteringHistory::AddStep(
        GetHistory(candidate),
        {newCandidate.pt(), newR, newCandidate.constituents().size(),
         ClusteringHistory::NONE});
    newCandidate.set_user_info(newHistory);

    return RecursiveRecluster(newCandidate, newR, mass, step + 1);
  } else {
    GetHistory(candidate).steps.back().status = ClusteringHistory::GOOD;
    return std::make_pair(true, candidate);
  }
}

AnalysisObject AnalysisClass::reclusteredParticle(const AnalysisObjects &jets,
                                                  const AnalysisObjects &bjets,
                                                  const double mass,
                                                  const bool useBJets) {
  AnalysisObject p =
      AnalysisObject(0., 0., 0., 0., 0, 0, AnalysisObjectType::JET, 0, 0);
  double r0 = 3.0;

  auto usejets = jets;
  if (useBJets && bjets.size())
    usejets = jets + bjets;

  std::vector<fastjet::PseudoJet> initialJets;
  for (const auto &jet : usejets) {
    fastjet::PseudoJet Pjet(jet.Px(), jet.Py(), jet.Pz(), jet.E());
    initialJets.push_back(Pjet);
  }

  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, r0);
  fastjet::ClusterSequence cs(initialJets, jetDef);

  auto candidates = fastjet::sorted_by_pt(cs.inclusive_jets());

  std::vector<fastjet::PseudoJet> selectedJets;
  selectedJets.reserve(candidates.size());
  std::vector<fastjet::PseudoJet> badJets;
  badJets.reserve(candidates.size());

  size_t i = 0;
  for (auto &cand : candidates) {
    auto history = new ClusteringHistory();
    history->id = i;
    history->steps.push_back(
        {cand.pt(), r0, cand.constituents().size(), ClusteringHistory::NONE});
    cand.set_user_info(history);
    ++i;
  }

  for (const auto &cand : candidates) {
    bool selected = false;
    fastjet::PseudoJet jet;

    std::tie(selected, jet) = RecursiveRecluster(cand, r0, mass, 0);

    if (selected)
      selectedJets.push_back(jet);
    else
      badJets.push_back(jet);
  }

  if (selectedJets.size() < 1) {
    return p;
  }

  AnalysisObjects aoSelectedJets;
  for (const auto &jet : selectedJets)
    aoSelectedJets.push_back(
        AnalysisObject(jet.px(), jet.py(), jet.pz(), jet.E(), 0, 0,
                       AnalysisObjectType::COMBINED, 0, 0));

  AnalysisClass::sortObjectsByPt(aoSelectedJets);
  p = aoSelectedJets[0];

  return p;
}

std::string FindFile(const std::string &name) {
  return PathResolverFindCalibFile("SimpleAnalysisCodes/" + name);
}

void printObject(const AnalysisObject &obj) {
  int id = obj.id();
  if (obj.type() == TRUTH)
    id = obj.pdgId();
  std::cout << " Pt: " << obj.Pt() << " Phi: " << obj.Phi()
            << " Eta:" << obj.Eta() << " charge: " << obj.charge()
            << " id: " << id << " type:" << obj.type() << std::endl;
}

void printObjects(const AnalysisObjects &objs) {
  int ii = 0;
  for (const auto &obj : objs) {
    std::cout << " " << ii;
    printObject(obj);
    ii++;
  }
}
