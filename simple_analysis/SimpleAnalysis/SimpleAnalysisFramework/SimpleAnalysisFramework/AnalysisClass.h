#ifndef ANALYSISCLASS_H
#define ANALYSISCLASS_H

#include "MVAUtils/BDT.h"
#include "SimpleAnalysisFramework/AnalysisEvent.h"
#include "SimpleAnalysisFramework/AnalysisObject.h"
#include "SimpleAnalysisFramework/OutputHandler.h"
#include <algorithm>
#include <functional>
#include <initializer_list>
#include <vector>

#include "SimpleAnalysisFramework/BDT.h"
#include "SimpleAnalysisFramework/MVA.h"

#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/NNLayerConfig.hh"
#include "lwtnn/lightweight_network_config.hh"
#include "lwtnn/parse_json.hh"

#include "SimpleAnalysisFramework/RestFramesHelper.h"

#include "TRandomGen.h"

#ifdef PACKAGE_BTaggingTruthTagging
#include "xAODBTaggingEfficiency/BTaggingTruthTaggingTool.h"
#include "xAODBTaggingEfficiency/TruthTagResults.h"
#endif

namespace lwt {
class LightweightNeuralNetwork;
}

#ifdef PACKAGE_BTaggingTruthTagging
class BTaggingTruthTaggingTool;
#endif

class AnalysisClass;

std::vector<AnalysisClass *> *
getAnalysisList(); // for automatically tracking which analyses have been
                   // defined

class AnalysisClass {
public:
  AnalysisClass(const std::string &name) : _name(name), _output(0) {
    getAnalysisList()->push_back(this);
  };
  AnalysisClass(){};
  virtual void Init(){};
  virtual void ProcessEvent(AnalysisEvent *event) = 0;
  virtual void Final(){};
  virtual ~AnalysisClass(){};
  virtual const std::string &name() { return _name; };
  virtual void setOutput(OutputHandler *output) { _output = output; };
  virtual void setEventNumber(int eventNumber) {
    _rndm.SetSeed(
        eventNumber * 1234 +
        2135); // Reproducible random numbers for use in some analysis codes
    for (auto mva : m_MVAs) {
      mva.second->setEventNumber(eventNumber);
    }
  };
  virtual TRandom *getRandom() { return &_rndm; };
  OutputHandler *getOutput() { return _output; };
  virtual void addRegion(const std::string &label) {
    _output->addEntry(label);
  };
  virtual void addRegions(const std::vector<std::string> &labels) {
    _output->addEntries(labels);
  };
  virtual void addHistogram(const std::string &label, int bins, float min,
                            float max) {
    _output->addHistogram(label, bins, min, max);
  };
  virtual void addHistogram(const std::string &label, int bins, float *edges) {
    _output->addHistogram(label, bins, edges);
  }
  virtual void addHistogram(const std::string &label,
                            std::vector<float> &edges) {
    _output->addHistogram(label, edges);
  }
  virtual void addHistogram(const std::string &label,
                            std::vector<std::string> charBins) {
    _output->addHistogram(label, charBins);
  }
  virtual void addHistogram(const std::string &label, int binsX, float minX,
                            float maxX, int binsY, float minY, float maxY) {
    _output->addHistogram(label, binsX, minX, maxX, binsY, minY, maxY);
  };

  virtual MVA *addTMVA(const std::string &name,
                       const std::vector<std::string> &variableDefs,
                       const std::string fname1, const std::string fname2 = "");

  virtual MVA *addMVAUtilsBDT(const std::string &name, const std::string fname1,
                              const std::string fname2 = "");

  virtual MVA *addONNX(const std::string &name, const std::string fname1,
                       const std::string fname2 = "");

  virtual MVA *getMVA(const std::string &name) {
    if (m_MVAs.find(name) == m_MVAs.end())
      throw std::runtime_error("Unknown MVA name");
    return m_MVAs[name];
  }

  virtual void accept(const std::string &name, double weight = 1) {
    _output->pass(name, weight);
  };
  virtual void fill(const std::string &name, double x) {
    _output->fillHistogram(name, x);
  };
  virtual void fill(const std::string &name, const char *bin) {
    _output->fillHistogram(name, bin);
  };
  virtual void fill(const std::string &name, double x, double y) {
    _output->fillHistogram(name, x, y);
  };
  virtual void fillWeighted(const std::string &name, double x, double y,
                            double val) {
    _output->fillHistogramWeighted(name, x, y, val);
  };
  void adjustEventWeight(double weight) { _output->adjustEventWeight(weight); };
  void ntupVar(const std::string &label, int value) {
    _output->ntupVar(label, value);
  };
  void ntupVar(const std::string &label, float value) {
    _output->ntupVar(label, value);
  };
  void ntupVar(const std::string &label, double value) {
    _output->ntupVar(label, (float)value);
  };
  void ntupVar(const std::string &label, std::vector<int> values) {
    _output->ntupVar(label, values);
  };
  void ntupVar(const std::string &label, std::vector<float> values) {
    _output->ntupVar(label, values);
  };
  void ntupVar(const std::string &label, AnalysisObject &object,
               bool saveMass = false, bool saveType = false,
               bool saveVtx = false, bool saveEnergy = false);
  void ntupVar(const std::string &label, AnalysisObjects &objects,
               bool saveMass = false, bool saveType = false,
               bool saveVtx = false, bool saveEnergy = false);

  // Helper functions
  std::vector<float> fakeJER(const AnalysisObjects &jets);
  static AnalysisObjects overlapRemoval(const AnalysisObjects &cands,
                                        const AnalysisObjects &others,
                                        float deltaR, int passId = 0);

  static AnalysisObjects overlapRemoval(
      const AnalysisObjects &cands, const AnalysisObjects &others,
      std::function<float(const AnalysisObject &, const AnalysisObject &)>
          radiusFunc,
      int passId = 0);

  static AnalysisObjects lowMassRemoval(
      const AnalysisObjects &cand,
      std::function<bool(const AnalysisObject &, const AnalysisObject &)>,
      float MinMass = 0, float MaxMass = FLT_MAX, int type = -1);

  static AnalysisObjects filterObjects(const AnalysisObjects &cands,
                                       float ptCut, float etaCut = 100.,
                                       int id = 0, unsigned int maxNum = 10000);

  static AnalysisObjects filterCrack(const AnalysisObjects &cands,
                                     float minEta = 1.37, float maxEta = 1.52);

  static AnalysisObjects filterObjectsRange(const AnalysisObjects &cands,
                                            float ptMinCut, float ptMaxCut,
                                            float etaMinCut, float etaMaxCut,
                                            int id = 0,
                                            unsigned int maxNum = 10000);

  static int countObjects(const AnalysisObjects &cands, float ptCut,
                          float etaCut = 100., int id = 0);

  static float sumObjectsPt(const AnalysisObjects &cands,
                            unsigned int maxNum = 10000, float ptCut = 0);

  static float sumObjectsM(const AnalysisObjects &cands,
                           unsigned int maxNum = 10000, float mCut = 0);

  static void sortObjectsByPt(AnalysisObjects &cands);
  static float minDphi(const AnalysisObject &met, const AnalysisObjects &cands,
                       unsigned int maxNum = 10000, float ptCut = 0);
  static float minDphi(const AnalysisObject &met, const AnalysisObject &cand);
  static float minDR(const AnalysisObject &cand, const AnalysisObjects &cands,
                     unsigned int maxNum = 10000, float ptCut = 0);
  static float minDR(const AnalysisObjects &cands, unsigned int maxNum = 10000,
                     float ptCut = 0);

  static float calcMCT(const AnalysisObject &o1, const AnalysisObject &o2);
  static float calcMCT(const AnalysisObject &o1, const AnalysisObject &o2,
                       const AnalysisObject &met, double ecm = 13000.0);
  static float calcMT(const AnalysisObject &lepton, const AnalysisObject &met);
  static float calcMTmin(const AnalysisObjects &cands,
                         const AnalysisObject &met, int maxNum = 10000);
  static float calcMT2(const AnalysisObject &o1, const AnalysisObject &o2,
                       const AnalysisObject &met);
  static float calcAMT2(const AnalysisObject &o1, const AnalysisObject &o2,
                        const AnalysisObject &met, float m1, float m2);
  static float calcMTauTau(const AnalysisObject &o1, const AnalysisObject &o2,
                           const AnalysisObject &met);
  static float calcTopness(const AnalysisObject &p_lep,
                           const AnalysisObject &p_met_xy,
                           const AnalysisObject &p_jet1,
                           const AnalysisObject &p_jet2,
                           const AnalysisObject &p_jet3);
  static float calcTopness(const AnalysisObject &p_lep,
                           const AnalysisObject &p_met_xy,
                           const AnalysisObject &p_jet1,
                           const AnalysisObject &p_jet2);

  static float aplanarity(const AnalysisObjects &jets);
  static float sphericity(const AnalysisObjects &jets);

  static bool IsSFOS(const AnalysisObject &L, const AnalysisObject &L1);
  static std::pair<float, float> DiZSelection(const AnalysisObjects &electrons,
                                              const AnalysisObjects &muons);
  static bool PassZVeto(const AnalysisObjects &electrons,
                        const AnalysisObjects &Muons, float Window = 10);

  static AnalysisObjects reclusterJets(const AnalysisObjects &jets,
                                       float radius, float ptmin,
                                       float rclus = -1, float ptfrac = -1);

  static AnalysisObject reclusteredParticle(const AnalysisObjects &jets,
                                            const AnalysisObjects &bjets,
                                            const double mass,
                                            const bool useBJets);

  std::vector<BDT::BDTReader *> m_BDTReaders;
  MVAUtils::BDT *m_MVAUtilsBDT;

  std::unique_ptr<lwt::LightweightGraph> RNN_graph;

  std::map<std::string, MVA *> m_MVAs;

  RestFramesHelper m_RF_helper;

#ifdef PACKAGE_BTaggingTruthTagging
  BTaggingTruthTaggingTool *m_btt;
#endif

protected:
  std::string _name;
  OutputHandler *_output;
  TRandomMT64 _rndm;
};

#define DefineAnalysis(ANALYSISNAME)                                           \
  class ANALYSISNAME : public AnalysisClass {                                  \
  public:                                                                      \
    ANALYSISNAME() : AnalysisClass(#ANALYSISNAME){};                           \
    void Init();                                                               \
    void ProcessEvent(AnalysisEvent *event);                                   \
  };                                                                           \
  static const AnalysisClass *ANALYSISNAME_instance __attribute__((used)) =    \
      new ANALYSISNAME();

std::string FindFile(const std::string &name);

// debugging functions

void printObject(const AnalysisObject &obj);

void printObjects(const AnalysisObjects &objs);

#define DebugObject(obj, num)                                                  \
  do {                                                                         \
    static int cnt = 0;                                                        \
    if (num || cnt < num)                                                      \
      printObject(obj);                                                        \
    cnt++;                                                                     \
  } while (0)

#define DebugObjects(objs, num)                                                \
  do {                                                                         \
    static int cnt = 0;                                                        \
    if (num || cnt < num)                                                      \
      printObjects(objs);                                                      \
    cnt++;                                                                     \
  } while (0)

#endif
