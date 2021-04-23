#ifndef READER_H
#define READER_H

#include <boost/program_options.hpp>
#include <string>
#include <vector>
namespace po = boost::program_options;

#include "SimpleAnalysisFramework/AnalysisRunner.h"

class AnalysisClass;
class Smearer;
class Reweighter;
class TFile;

class Reader;
std::vector<Reader *> *
getReaderList(); // for automatically tracking available readers

class Reader {
public:
  Reader() { getReaderList()->push_back(this); };
  virtual ~Reader(){};
  virtual bool isFileType(TFile * /*file*/, po::variables_map & /* vm */) {
    return false;
  };
  virtual void Init(){};
  virtual void AddOptions(po::options_description & /*desc */){};
  virtual void SetAnalyses(std::vector<AnalysisClass *> &analysisList);
  virtual void SetSmearing(Smearer *smear);
  virtual void SetMultiRuns(int runs);
  virtual void AddReweighting(Reweighter *reweighter);
  virtual void SetMCWeightIndex(int mcwidx);

  virtual int getMCWeightIndex();
  virtual void processFiles(const std::vector<std::string> &inputNames,
                            unsigned int nevents = -1);

protected:
  virtual void processFilesInternal(const std::vector<std::string> &inputNames,
                                    unsigned int nevents) = 0;
  AnalysisRunner *_analysisRunner;

private:
};

#define DefineReader(NAME)                                                     \
  static const Reader *NAME_instance __attribute__((used)) = new NAME()

#endif
