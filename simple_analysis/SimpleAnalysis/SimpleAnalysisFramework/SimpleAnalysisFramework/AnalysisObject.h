#ifndef ANALYSISOBJECT_H
#define ANALYSISOBJECT_H

#include "TLorentzVector.h"
#include <vector>

enum AnalysisObjectType { // Please don't change order
  ELECTRON,
  MUON,
  TAU,
  PHOTON,
  JET,
  FATJET,
  MET,
  COMBINED,
  NONE,
  TRUTH,
  BHADRON,
  TRACKJET
};

enum AnalysisCommonID { NotBit = 1 << 31 };
#define NOT(x) (NotBit | x)

enum AnalysisElectronID {
  EVeryLooseLH = 1 << 0,
  ELooseLH = 1 << 1,
  EMediumLH = 1 << 2,
  ETightLH = 1 << 3,
  ELooseBLLH = 1 << 4,
  EIsoGradientLoose = 1 << 8,
  EIsoBoosted = 1 << 9, // from 3-bjet paper
  EIsoFixedCutTight = 1 << 10,
  EIsoLooseTrack = 1 << 11,
  EIsoLoose = 1 << 12,
  EIsoGradient = 1 << 13,
  EIsoFixedCutLoose = 1 << 14,
  EIsoFixedCutTightTrackOnly = 1 << 15,
  ED0Sigma5 = 1 << 16,
  EZ05mm = 1 << 17,
  EIsoFCTight = 1 << 18,
  EIsoFCTightTrackOnly = 1 << 19,
  EIsoTightTrackOnly = 1 << 19,
  EIsoFCHighPtCaloOnly = 1 << 20,
  EIsoTightTrackOnly_FixedRad = 1 << 21,
  EIsoFCLoose = 1 << 22,
  EIsoPLVTight = 1 << 23,
  EGood = EVeryLooseLH | ELooseLH | EMediumLH | ETightLH | ELooseBLLH |
          ED0Sigma5 | EZ05mm,
  EIsoGood = EGood | EIsoGradientLoose | EIsoBoosted | EIsoFixedCutTight |
             EIsoFCTight | EIsoLooseTrack | EIsoLoose | EIsoGradient |
             EIsoFixedCutLoose | EIsoFixedCutTightTrackOnly |
             EIsoFCTightTrackOnly | EIsoFCHighPtCaloOnly |
             EIsoTightTrackOnly_FixedRad | EIsoFCLoose | EIsoPLVTight
};

enum AnalysisMuonID {
  MuLoose = 1 << 0,
  MuMedium = 1 << 1,
  MuTight = 1 << 2,
  MuVeryLoose = 1 << 3,
  MuHighPt = 1 << 4,
  MuIsoGradientLoose = 1 << 8,
  MuIsoBoosted = 1 << 9, // from 3-bjet paper
  MuIsoFixedCutTightTrackOnly = 1 << 10,
  MuIsoLooseTrack = 1 << 11,
  MuIsoLoose = 1 << 12,
  MuIsoGradient = 1 << 13,
  MuIsoFixedCutLoose = 1 << 14,
  MuD0Sigma3 = 1 << 16,
  MuZ05mm = 1 << 17,
  MuNotCosmic = 1 << 18,
  MuQoPSignificance = 1 << 19,
  MuCaloTaggedOnly = 1 << 20,
  MuIsoFCTightTrackOnly = 1 << 21,
  MuIsoFCTightFR = 1 << 22,
  MuIsoPflowTight_VarRad = 1 << 23,
  MuIsoPflowTight_FixedRad = 1 << 24,
  MuIsoPflowLoose_VarRad = 1 << 25,
  MuIsoPflowLoose_FixedRad = 1 << 26,
  MuIsoHighPtTrackOnly = 1 << 27,
  MuIsoTightTrackOnly_VarRad = 1 << 28,
  MuIsoTightTrackOnly_FixedRad = 1 << 29,
  MuIsoFCLoose = 1 << 30,
  MuIsoPLVTight = 1 << 15,

  MuGood = MuLoose | MuMedium | MuTight | MuHighPt | MuVeryLoose | MuD0Sigma3 |
           MuZ05mm | MuNotCosmic | MuQoPSignificance,
  MuIsoGood = MuGood | MuIsoGradientLoose | MuIsoBoosted |
              MuIsoFixedCutTightTrackOnly | MuIsoFCTightTrackOnly |
              MuIsoLooseTrack | MuIsoLoose | MuIsoGradient |
              MuIsoFixedCutLoose | MuIsoFCTightFR | MuIsoPflowTight_VarRad |
              MuIsoPflowTight_FixedRad | MuIsoPflowLoose_VarRad |
              MuIsoPflowLoose_FixedRad | MuIsoHighPtTrackOnly |
              MuIsoTightTrackOnly_VarRad | MuIsoTightTrackOnly_FixedRad |
              MuIsoFCLoose | MuIsoPLVTight
};

enum AnalysisTauID {
  TauLoose = 1 << 0,
  TauBDTLoose = 1 << 0,
  TauMedium = 1 << 1,
  TauBDTMedium = 1 << 1,
  TauTight = 1 << 2,
  TauBDTTight = 1 << 2,
  TauRNNVeryLoose = 1 << 3,
  TauRNNLoose = 1 << 4,
  TauRNNMedium = 1 << 5,
  TauRNNTight = 1 << 6,
  TauOneProng = 1 << 10,
  TauThreeProng = 1 << 11,
  TauGood = TauLoose | TauMedium | TauTight | TauRNNVeryLoose | TauRNNLoose |
            TauRNNMedium | TauRNNTight,
  TauIsoGood = TauGood
};

enum AnalysisPhotonID {
  PhotonLoose = 1 << 0,
  PhotonTight = 1 << 1,
  PhotonIsoFixedCutLoose = 1 << 8,
  PhotonIsoFixedCutTight = 1 << 9,
  PhotonIsoFixedCutTightCaloOnly = 1 << 10,
  PhotonGood = PhotonLoose | PhotonTight,
  PhotonIsoGood = PhotonGood | PhotonIsoFixedCutLoose | PhotonIsoFixedCutTight |
                  PhotonIsoFixedCutTightCaloOnly
};

enum AnalysisJetID {
  LooseBadJet = 1 << 8,
  TightBadJet = 1 << 9,
  JVT50Jet = 1 << 10,
  LessThan3Tracks = 1 << 11, // most signal jets should fail this
  JVT59Jet = 1 << 12,
  LessOrEq4Tracks = 1 << 13, // For ttbarMET0L2019 tau veto
  JVT120Jet = 1 << 14,
  PFlowJet = 1 << 15,
  JVTLoose = 1 << 16,
  JVTMedium = 1 << 17,
  JVTTight = 1 << 18,
  GoodJet = LooseBadJet | TightBadJet | JVT50Jet | JVT59Jet | JVT120Jet |
            PFlowJet | JVTLoose | JVTMedium | JVTTight,
  BTag85MV2c20 = 1 << 0,
  BTag80MV2c20 = 1 << 1,
  BTag77MV2c20 = 1 << 2,
  BTag70MV2c20 = 1 << 3,
  BTag85DL1 = 1 << 0,
  BTag77DL1 = 1 << 1,
  BTag70DL1 = 1 << 2,
  BTag60DL1 = 1 << 3,
  BTag85MV2c10 = 1 << 4,
  BTag77MV2c10 = 1 << 5,
  BTag70MV2c10 = 1 << 6,
  BTag60MV2c10 = 1 << 7,
  BTag85DL1r = 1 << 20,
  BTag77DL1r = 1 << 21,
  BTag70DL1r = 1 << 22,
  BTag60DL1r = 1 << 23,
  BTag40OnlineMV2 = 1 << 24,
  BTag60OnlineMV2 = 1 << 25,
  GoodBJet =
      BTag85MV2c20 | BTag80MV2c20 | BTag77MV2c20 |
      BTag70MV2c20 | // note DL1 reuse the same bits as cannot be in same data
      BTag85MV2c10 | BTag77MV2c10 | BTag70MV2c10 | BTag60MV2c10 | BTag85DL1r |
      BTag77DL1r | BTag70DL1r | BTag60DL1r | GoodJet,
  TrueLightJet = 1 << 27, // These should not be used for actual analysis
                          // selection, only for understanding
  TrueCJet = 1 << 28,
  TrueBJet = 1 << 29,
  TrueTau = 1 << 30,

};

enum AnalysisFatJetID { LooseFatJet = 1 << 8, GoodFatJet = LooseFatJet };

enum AnalysisTrackJetID {
  LooseTrackJet = 1 << 8,
  GoodTrackJet = LooseTrackJet
};

enum ParticleIDs {
  // PDGIDs
  Neutralino1 = 1000022,
  Neutralino2 = 1000023,
  Neutralino3 = 1000025,
  Neutralino4 = 1000035,
  Chargino1 = 1000024,
  Chargino2 = 1000037,
  // Status codes
  StablePart = 1 << 24,
  DecayedPart = 2 << 24
};

class AnalysisObject : public TLorentzVector {
public:
  AnalysisObject(double Px, double Py, double Pz, double E, int charge, int id,
                 AnalysisObjectType type, int motherID, int orgIndex)
      : TLorentzVector(Px, Py, Pz, E), _charge(charge), _id(id), _type(type),
        _motherID(motherID), _orgIndex(orgIndex){};
  AnalysisObject(TLorentzVector tlv, int charge, int id,
                 AnalysisObjectType type, int motherID, int orgIndex)
      : TLorentzVector(tlv), _charge(charge), _id(id), _type(type),
        _motherID(motherID), _orgIndex(orgIndex){};
  virtual bool pass(int id) const {
    if (_type == TRUTH) {
      int mask = (id >> 24) ? 0x7FFFFFFF : 0xFFFFFF;
      return ((_id & mask) == id) || id == 0;
    }
    if (id & NotBit)
      return (id & _id) == 0;
    else
      return (id & _id) == id;
  };
  virtual int charge() const { return _charge; };
  virtual AnalysisObjectType type() const { return _type; };
  virtual bool valid() const { return _type != NONE; };
  virtual int id() const { return _id; };
  virtual int pdgId() const {
    if (_type != TRUTH)
      throw std::runtime_error("Can only call pdgId() for truth objects");
    int pdg = _id & 0xFFFFFF;
    if (_id & (1 << 31))
      pdg = -1 * pdg;
    return pdg;
  };
  virtual int status() const {
    if (_type != TRUTH)
      throw std::runtime_error("Can only call status() for truth objects");
    return (_id >> 24) & 0x7F;
  };
  virtual int motherID() const {
    return _motherID;
  }; // not supposed to be used directly except to store in ntuples
  virtual int index() const { return _orgIndex; }; // for internal use
  virtual AnalysisObject transFourVect() const {
    TLorentzVector tlv;
    tlv.SetPtEtaPhiM(Pt(), 0.0, Phi(), M());
    return AnalysisObject(tlv, charge(), id(), type(), motherID(), _orgIndex);
  };
  virtual TVector3 prodVtx() const {
    TVector3 zero;
    if (_vertices.size() > 0)
      return _vertices[0];
    return zero;
  };
  virtual TVector3 decayVtx() const {
    TVector3 zero;
    if (_vertices.size() > 1)
      return _vertices[1];
    return zero;
  };
  virtual void setProdVtx(const TVector3 &vtx) {
    if (_vertices.size())
      _vertices[0] = vtx;
    else
      _vertices.push_back(vtx);
  };
  virtual void setId(const int id) { // not meant for analysis use
    _id = id;
  }
  virtual void setDecayVtx(const TVector3 &vtx) {
    if (_vertices.size() > 1)
      _vertices[1] = vtx;
    else if (_vertices.size() > 0)
      _vertices.push_back(vtx);
    else {
      TVector3 zero;
      _vertices.push_back(zero);
      _vertices.push_back(vtx);
    }
  };

private:
  int _charge; // not used for jets or photons
  int _id;
  AnalysisObjectType _type;
  int _motherID;
  int _orgIndex;
  std::vector<TVector3> _vertices;
};

AnalysisObject operator+(const AnalysisObject &lhs, const AnalysisObject &rhs);

// typedef std::vector<AnalysisObject> AnalysisObjects;

class AnalysisObjects : public std::vector<AnalysisObject> {
  using std::vector<AnalysisObject>::vector; // reuse vector initializers
  // force range check when accessing objects directly
public:
  AnalysisObject &operator[](std::vector<AnalysisObject>::size_type ii) {
    return this->at(ii);
  };
  const AnalysisObject &
  operator[](std::vector<AnalysisObject>::size_type ii) const {
    return this->at(ii);
  };
};

AnalysisObjects operator+(const AnalysisObjects &lhs,
                          const AnalysisObjects &rhs);

#endif
