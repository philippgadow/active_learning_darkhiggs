#ifndef PDFREWEIGHT_H
#define PDFREWEIGHT_H

#include <map>
#include <string>
#include <vector>

#include "TRandom3.h"

#include "SimpleAnalysisFramework/Reweight.h"

class AnalysisEvent;

class TruthDecayer : public Reweighter {
public:
  TruthDecayer();
  void init(std::vector<std::string> &options);
  double reweightEvent(AnalysisEvent *event);

private:
  TRandom3 _random;
  int _status;
  std::map<int, float> _decays;
};

#endif
