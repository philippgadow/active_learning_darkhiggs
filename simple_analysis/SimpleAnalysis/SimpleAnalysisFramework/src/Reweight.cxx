#include "SimpleAnalysisFramework/Reweight.h"
#include <vector>

std::vector<Reweighter *> *getReweighterList() {
  static std::vector<Reweighter *> *list = new std::vector<Reweighter *>;
  return list;
}
