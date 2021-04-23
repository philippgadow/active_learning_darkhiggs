#include "SimpleAnalysisFramework/Reader.h"
#include "SimpleAnalysisFramework/AnalysisClass.h"
#include "SimpleAnalysisFramework/AnalysisRunner.h"
#include "SimpleAnalysisFramework/Reweight.h"
#include "SimpleAnalysisFramework/Smearing.h"

std::vector<Reader *> *getReaderList() {
  static std::vector<Reader *> *list = new std::vector<Reader *>;
  return list;
}

void Reader::SetAnalyses(std::vector<AnalysisClass *> &analysisList) {
  _analysisRunner = new AnalysisRunner(analysisList);
}

void Reader::SetSmearing(Smearer *smear) {
  _analysisRunner->SetSmearing(smear);
}

void Reader::SetMultiRuns(int runs) { _analysisRunner->SetMultiRuns(runs); }

void Reader::AddReweighting(Reweighter *reweighter) {
  _analysisRunner->AddReweighting(reweighter);
}

void Reader::SetMCWeightIndex(int mcwidx) {
  _analysisRunner->SetMCWeightIndex(mcwidx);
}

int Reader::getMCWeightIndex() { return _analysisRunner->getMCWeightIndex(); }

void Reader::processFiles(const std::vector<std::string> &inputNames,
                          unsigned int nevents) {
  _analysisRunner->init();
  processFilesInternal(inputNames, nevents);
  _analysisRunner->final();
}
