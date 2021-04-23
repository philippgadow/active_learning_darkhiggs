#ifndef TRUTHEVENT_H
#define TRUTHEVENT_H

#include "SimpleAnalysisFramework/AnalysisClass.h"
#include "SimpleAnalysisFramework/METSignificanceCalculator.h"
#include <cmath>

class TruthDecayer;

class TruthEvent : public AnalysisEvent {
  friend TruthDecayer;

public:
  TruthEvent(float sumet, double Ex, double Ey)
      : _sumet(sumet),
        _met(Ex, Ey, 0, sqrt(Ex * Ex + Ey * Ey), 0, 0, MET, 0, 0),
        _applySmearing(false), _truth(0){};
  virtual AnalysisObjects getElectrons(float ptCut, float etaCut,
                                       int isolation) {
    return AnalysisClass::filterObjects(_baseElectrons, ptCut, etaCut,
                                        isolation);
  };
  virtual AnalysisObjects getMuons(float ptCut, float etaCut, int isolation) {
    return AnalysisClass::filterObjects(_baseMuons, ptCut, etaCut, isolation);
  };
  virtual AnalysisObjects getTaus(float ptCut, float etaCut, int isolation) {
    return AnalysisClass::filterObjects(_baseTaus, ptCut, etaCut, isolation);
  };
  virtual AnalysisObjects getPhotons(float ptCut, float etaCut, int isolation) {
    return AnalysisClass::filterObjects(_basePhotons, ptCut, etaCut, isolation);
  };
  virtual AnalysisObjects getJets(float ptCut, float etaCut, int btag) {
    return AnalysisClass::filterObjects(_baseJets, ptCut, etaCut, btag);
  };
  virtual AnalysisObjects getFatJets(float ptCut, float etaCut, int btag) {
    return AnalysisClass::filterObjects(_baseFatJets, ptCut, etaCut, btag);
  };
  virtual AnalysisObjects getTrackJets(float ptCut, float etaCut, int btag) {
    return AnalysisClass::filterObjects(_baseTrackJets, ptCut, etaCut, btag);
  };
  virtual AnalysisObjects getBHadrons(float ptCut, float etaCut, int pdgId) {
    return AnalysisClass::filterObjects(_baseBHadrons, ptCut, etaCut, pdgId);
  };
  virtual AnalysisObjects getHSTruth(float ptCut, float etaCut, int pid) {
    return AnalysisClass::filterObjects(_HSTruth, ptCut, etaCut, pid);
  };
  virtual AnalysisObject getMET() { return _met; };
  virtual unsigned int getMCVeto() { return _veto; };
  virtual float getMETSignificance(bool applyOverlapRemoval = false) {
    return calcMETSignificance(this, applyOverlapRemoval);
  };
  virtual float getMETSignificance(AnalysisObjects &electrons,
                                   AnalysisObjects &photons,
                                   AnalysisObjects &muons,
                                   AnalysisObjects &jets, AnalysisObjects &taus,
                                   AnalysisObject &metVec) {
    return calcMETSignificance(electrons, photons, muons, jets, taus, metVec);
  };

  virtual float getSumET() { return _sumet; };

  virtual float getGenMET() { return _genmet; };
  virtual float getGenHT() { return _genht; };

  virtual int getMCNumber() {
    return _mcChannel;
  }; // Temporary until better solution found
  virtual int getSUSYChannel() {
    return _susyChannel;
  }; // Temporary until better solution found
  virtual std::vector<float> getMCWeights() {
    return _mcWeights;
  }; // Temporary until better solution found

  virtual bool applyLocalSmearing() { return _applySmearing; };

  virtual AnalysisObject getTruthParticle(const AnalysisObject &obj) {
    if (_truth && (obj.index() >= 0))
      switch (obj.type()) {
      case ELECTRON:
        return _truth->_baseElectrons[obj.index()];
      case MUON:
        return _truth->_baseMuons[obj.index()];
      case TAU:
        return _truth->_baseTaus[obj.index()];
      case PHOTON:
        return _truth->_basePhotons[obj.index()];
      case JET:
        return _truth->_baseJets[obj.index()];
      case FATJET:
        return _truth->_baseFatJets[obj.index()];
      case BHADRON:
        return _truth->_baseBHadrons[obj.index()];
      case MET:
        return _truth->_met;
      case TRACKJET:
        return _truth->_baseTrackJets[obj.index()];
      case COMBINED: // No truth for these
      case NONE:
      case TRUTH:
        break;
      }
    return AnalysisObject(0, 0, 0, 0, 0, 0, NONE, 0, -1);
  };

  virtual void setChannelInfo(int mcChannel, int susyChannel) {
    _mcChannel = mcChannel;
    _susyChannel = susyChannel;
  };
  virtual void setChannelInfo(int mcChannel, unsigned int veto,
                              int susyChannel) {
    _mcChannel = mcChannel;
    _veto = veto;
    _susyChannel = susyChannel;
  };
  virtual int getPDF_id1() { return _id1; };
  virtual float getPDF_x1() { return _x1; };
  virtual float getPDF_pdf1() { return _pdf1; };
  virtual int getPDF_id2() { return _id2; };
  virtual float getPDF_x2() { return _x2; };
  virtual float getPDF_pdf2() { return _pdf2; };
  virtual float getPDF_scale() { return _scale; };
  virtual void setPDFInfo(int id1, float x1, float pdf1, int id2, float x2,
                          float pdf2, float scale) {
    _id1 = id1;
    _x1 = x1;
    _pdf1 = pdf1;
    _id2 = id2;
    _x2 = x2;
    _pdf2 = pdf2;
    _scale = scale;
  };

  void addElectron(double Px, double Py, double Pz, double E, int charge,
                   int iso, int motherID, int idx) {
    _baseElectrons.push_back(
        AnalysisObject(Px, Py, Pz, E, charge, iso, ELECTRON, motherID, idx));
  };
  void addElectron(TLorentzVector tlv, int charge, int iso, int motherID,
                   int idx) {
    _baseElectrons.push_back(
        AnalysisObject(tlv, charge, iso, ELECTRON, motherID, idx));
  };
  void addMuon(double Px, double Py, double Pz, double E, int charge, int iso,
               int motherID, int idx) {
    _baseMuons.push_back(
        AnalysisObject(Px, Py, Pz, E, charge, iso, MUON, motherID, idx));
  };
  void addMuon(TLorentzVector tlv, int charge, int iso, int motherID, int idx) {
    _baseMuons.push_back(AnalysisObject(tlv, charge, iso, MUON, motherID, idx));
  };
  void addTau(double Px, double Py, double Pz, double E, int charge, int iso,
              int motherID, int idx) {
    _baseTaus.push_back(
        AnalysisObject(Px, Py, Pz, E, charge, iso, TAU, motherID, idx));
  };
  void addTau(TLorentzVector tlv, int charge, int iso, int motherID, int idx) {
    _baseTaus.push_back(AnalysisObject(tlv, charge, iso, TAU, motherID, idx));
  };
  void addPhoton(double Px, double Py, double Pz, double E, int iso,
                 int motherID, int idx) {
    _basePhotons.push_back(
        AnalysisObject(Px, Py, Pz, E, 0, iso, PHOTON, motherID, idx));
  };
  void addPhoton(TLorentzVector tlv, int iso, int motherID, int idx) {
    _basePhotons.push_back(AnalysisObject(tlv, 0, iso, PHOTON, motherID, idx));
  };
  void addJet(double Px, double Py, double Pz, double E, int iso, int idx) {
    _baseJets.push_back(AnalysisObject(Px, Py, Pz, E, 0, iso, JET, 0, idx));
  };
  void addJet(TLorentzVector tlv, int iso, int idx) {
    _baseJets.push_back(AnalysisObject(tlv, 0, iso, JET, 0, idx));
  };
  void addFatJet(double Px, double Py, double Pz, double E, int iso, int idx) {
    _baseFatJets.push_back(
        AnalysisObject(Px, Py, Pz, E, 0, iso, FATJET, 0, idx));
  };
  void addFatJet(TLorentzVector tlv, int iso, int idx) {
    _baseFatJets.push_back(AnalysisObject(tlv, 0, iso, FATJET, 0, idx));
  };
  void addBHadron(TLorentzVector tlv, int pdgId, int idx) {
    _baseBHadrons.push_back(AnalysisObject(tlv, 0, pdgId, BHADRON, 0, idx));
  };
  void addTrackJet(double Px, double Py, double Pz, double E, int iso,
                   int idx) {
    _baseTrackJets.push_back(
        AnalysisObject(Px, Py, Pz, E, 0, iso, TRACKJET, 0, idx));
  };
  void addTrackJet(TLorentzVector tlv, int iso, int idx) {
    _baseTrackJets.push_back(AnalysisObject(tlv, 0, iso, TRACKJET, 0, idx));
  };

  AnalysisObject *addHSTruth(double Px, double Py, double Pz, double E,
                             int charge, int id, int motherID, int idx) {
    _HSTruth.push_back(
        AnalysisObject(Px, Py, Pz, E, charge, id, TRUTH, motherID, idx));
    return &_HSTruth[_HSTruth.size() - 1];
  };
  AnalysisObject *addHSTruth(TLorentzVector tlv, int charge, int id,
                             int motherID, int idx) {
    _HSTruth.push_back(AnalysisObject(tlv, charge, id, TRUTH, motherID, idx));
    return &_HSTruth[_HSTruth.size() - 1];
  };
  void addHSTruthList(AnalysisObjects truthList) {
    for (const auto &truth : truthList)
      _HSTruth.push_back(truth);
  };
  virtual void setMET(float Ex = 0., float Ey = 0.) {
    _met.SetPxPyPzE(Ex, Ey, 0, sqrt(Ex * Ex + Ey * Ey));
  };

  virtual void setGenMET(float genMET = 0.) { _genmet = genMET; };

  virtual void setGenHT(float genHT = 0.) { _genht = genHT; };

  virtual void setSmearing(bool smear = true) { _applySmearing = smear; };

  virtual void setMCWeights(std::vector<float> ws) { _mcWeights = ws; }

  virtual void setTruth(AnalysisEvent *event) {
    _truth = dynamic_cast<TruthEvent *>(event);
  };

private:
  AnalysisObjects _baseElectrons;
  AnalysisObjects _baseMuons;
  AnalysisObjects _baseTaus;
  AnalysisObjects _basePhotons;
  AnalysisObjects _baseJets;
  AnalysisObjects _baseFatJets;
  AnalysisObjects _baseBHadrons;
  AnalysisObjects _baseTrackJets;
  AnalysisObjects _HSTruth;
  float _sumet;
  AnalysisObject _met;
  unsigned int _veto;
  float _genmet;
  float _genht;

  int _mcChannel;
  int _susyChannel;
  int _id1;
  float _x1;
  float _pdf1;
  int _id2;
  float _x2;
  float _pdf2;
  float _scale;
  bool _applySmearing;

  std::vector<float> _mcWeights;
  TruthEvent *_truth;

public:
  void sortObjects() {
    AnalysisClass::sortObjectsByPt(_baseElectrons);
    AnalysisClass::sortObjectsByPt(_baseMuons);
    AnalysisClass::sortObjectsByPt(_baseTaus);
    AnalysisClass::sortObjectsByPt(_basePhotons);
    AnalysisClass::sortObjectsByPt(_baseJets);
    AnalysisClass::sortObjectsByPt(_baseFatJets);
    AnalysisClass::sortObjectsByPt(_baseBHadrons);
    AnalysisClass::sortObjectsByPt(_baseTrackJets);
  }
};

#endif
