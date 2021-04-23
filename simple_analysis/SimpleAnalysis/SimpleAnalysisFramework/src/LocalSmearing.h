#ifndef LOCALSMEAR_H
#define LOCALSMEAR_H

#include "SimpleAnalysisFramework/Smearing.h"

class TruthEvent;
class LocalSmear : public Smearer {
public:
  LocalSmear();
  void init(std::vector<std::string> &);
  TruthEvent *smearEvent(AnalysisEvent *);

private:
  bool _apply;
};

#endif
