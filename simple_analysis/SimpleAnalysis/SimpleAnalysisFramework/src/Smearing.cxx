#include "SimpleAnalysisFramework/Smearing.h"
#include <vector>

std::vector<Smearer *> *getSmearerList() {
  static std::vector<Smearer *> *list = new std::vector<Smearer *>;
  return list;
}
