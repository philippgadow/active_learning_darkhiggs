#ifndef ANALYSISRUNNER_H
#define ANALYSISRUNNER_H

#include <string>
#include <vector>

#include "SimpleAnalysisFramework/AnalysisClass.h"
#include "SimpleAnalysisFramework/OutputHandler.h"
#include "SimpleAnalysisFramework/Reweight.h"
#include "SimpleAnalysisFramework/Smearing.h"
#include "SimpleAnalysisFramework/TruthEvent.h"

class AnalysisRunner {
public:
  AnalysisRunner(std::vector<AnalysisClass *> &analysisList)
      : _analysisList(analysisList), _smear(0), _reweighter(0), _runs(1){};
  void init() {
    for (const auto &analysis : _analysisList) {
      analysis->getOutput()->init();
      analysis->Init();
    }
  };
  void final() {
    for (const auto &analysis : _analysisList)
      analysis->Final();
  };
  void SetSmearing(Smearer *smear) { _smear = smear; };
  void SetMultiRuns(int runs) { _runs = runs; };
  void AddReweighting(Reweighter *reweighter) {
    _reweighter.push_back(reweighter);
  };
  void SetMCWeightIndex(int mcwidx) { _mcwindex = mcwidx; };

  int getMCWeightIndex() { return _mcwindex; };

  void processEvent(TruthEvent *event, int eventNumber);

private:
  std::vector<AnalysisClass *> &_analysisList;
  Smearer *_smear;
  std::vector<Reweighter *> _reweighter;
  int _mcwindex;
  int _runs;
};

#endif
