#ifndef ANALYSISEVENT_H
#define ANALYSISEVENT_H
#include "SimpleAnalysisFramework/AnalysisObject.h"

class AnalysisEvent {
public:
  virtual AnalysisObjects getElectrons(float ptCut, float etaCut,
                                       int isolation = 1) = 0;
  virtual AnalysisObjects getMuons(float ptCut, float etaCut,
                                   int isolation = 1) = 0;
  virtual AnalysisObjects getTaus(float ptCut, float etaCut,
                                  int isolation = 1) = 0;
  virtual AnalysisObjects getPhotons(float ptCut, float etaCut,
                                     int isolation = 1) = 0;
  virtual AnalysisObjects getJets(float ptCut, float etaCut, int btag = 0) = 0;
  virtual AnalysisObjects getFatJets(float ptCut, float etaCut,
                                     int btag = 0) = 0;
  virtual AnalysisObjects getBHadrons(float ptCut = 0, float etaCut = 100,
                                      int pdgId = 0) = 0;
  virtual AnalysisObjects getTrackJets(float ptCut, float etaCut,
                                       int btag = 0) = 0;
  virtual AnalysisObjects getHSTruth(float ptCut = 0, float etaCut = 100,
                                     int pid = 0) = 0;
  virtual AnalysisObject getMET() = 0;
  virtual unsigned int getMCVeto() = 0;
  virtual float getMETSignificance(bool applyOverlapRemoval = false) = 0;
  virtual float getMETSignificance(AnalysisObjects &electrons,
                                   AnalysisObjects &photons,
                                   AnalysisObjects &muons,
                                   AnalysisObjects &jets, AnalysisObjects &taus,
                                   AnalysisObject &metVec) = 0;
  virtual float getSumET() = 0;
  virtual float getGenMET() = 0;
  virtual float getGenHT() = 0;
  virtual int getMCNumber() = 0;    // Temporary until better solution found
  virtual int getSUSYChannel() = 0; // Temporary until better solution found
  virtual bool applyLocalSmearing() = 0;
  virtual int getPDF_id1() = 0;
  virtual float getPDF_x1() = 0;
  virtual float getPDF_pdf1() = 0;
  virtual int getPDF_id2() = 0;
  virtual float getPDF_x2() = 0;
  virtual float getPDF_pdf2() = 0;
  virtual float getPDF_scale() = 0;
  virtual std::vector<float>
  getMCWeights() = 0; // Temporary until better solution found
  virtual AnalysisObject getTruthParticle(const AnalysisObject &obj) = 0;
  virtual ~AnalysisEvent(){};
};

#endif
