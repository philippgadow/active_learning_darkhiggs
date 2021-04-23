#include "LocalSmearing.h"
#include "SimpleAnalysisFramework/AnalysisEvent.h"
#include "SimpleAnalysisFramework/TruthEvent.h"

/*
   Does not apply any smearing on its own, but activatives any analysis specific
   smearing/efficiency
 */

DefineSmearer(LocalSmear, "localSmear,L",
              "Apply analysis specific smearing [true|false]");

void LocalSmear::init(std::vector<std::string> &options) {
  _apply = true;
  for (const auto &option : options) {
    if (option.find("help") == 0) {
      std::cout << "Specify 'true' or 'false'" << std::endl;
    }
    if (option.find("true") == 0)
      _apply = true;
    if (option.find("false") == 0)
      _apply = false;
  }
}

TruthEvent *LocalSmear::smearEvent(AnalysisEvent *event) {
  TruthEvent *truth = static_cast<TruthEvent *>(event);
  truth->setSmearing(_apply);
  return truth;
}
