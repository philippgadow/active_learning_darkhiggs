#ifndef EXECBASE_H
#define EXECBASE_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

#include <TFile.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "SimpleAnalysisFramework/AnalysisClass.h"
#include "SimpleAnalysisFramework/AnalysisRunner.h"
#include "SimpleAnalysisFramework/OutputHandler.h"
#include "SimpleAnalysisFramework/Reader.h"
#include "SimpleAnalysisFramework/Reweight.h"
#include "SimpleAnalysisFramework/Smearing.h"

class execBase {
public:
  execBase(const char *desc, const char *name = "", bool dumpText = false)
      : _desc(desc), _outputName(name), _mergedOutput(0), _dumpText(dumpText){};

protected:
  static void splitCommaString(const std::string &names,
                               std::vector<std::string> &result) {
    std::stringstream ss(names);
    while (ss.good()) {
      std::string substr;
      getline(ss, substr, ',');
      result.push_back(substr);
    }
  }

  virtual void AddOptions(){};
  virtual void SetupAnalysList() = 0;

public:
  int run(int argc, char **argv) {

    _desc.add_options()("help,h", "print usage and exit");
    if (_outputName.size())
      _desc.add_options()("output,o", po::value<std::string>(&_outputName),
                          ("Output name - default: " + _outputName).c_str());
    else
      _desc.add_options()("output,o", po::value<std::string>(&_outputName),
                          "Output name - if not supplied use analyses names");

    _desc.add_options()("input-files", po::value<std::vector<std::string>>(),
                        "Comma-separated list of input files")("ntuple,n",
                                                               "Fill ntuple")(
        "multiRuns,M", po::value<int>()->default_value(1),
        "Run over each event multiple times - meant for smearing analysis "
        "only")("mcweight,w", po::value<int>()->default_value(0),
                "MC weight index to apply (set to -1 to ignore it, i.e. =1.)")(
        "nevents", po::value<int>()->default_value(-1),
        "number of events to run on (set to -1 to ignore it");
    this->AddOptions();
    for (auto *reader : *getReaderList()) {
      reader->AddOptions(_desc);
    }
    for (auto *reweighter : *getReweighterList()) {
      _desc.add_options()(reweighter->getOption().c_str(),
                          po::value<std::string>(),
                          reweighter->getDesc().c_str());
    }
    for (auto *smearer : *getSmearerList()) {
      _desc.add_options()(smearer->getOption().c_str(),
                          po::value<std::string>(), smearer->getDesc().c_str());
    }
    po::positional_options_description p;
    p.add("input-files", -1);

    po::store(
        po::command_line_parser(argc, argv).options(_desc).positional(p).run(),
        _vm);
    po::notify(_vm);

    _doNtuple = _vm.count("ntuple");

    int mcwindex = 0;
    if (_vm.count("mcweight")) {
      mcwindex = _vm["mcweight"].as<int>();
    }
    std::cout << "MCWeightIndex = " << mcwindex
              << (mcwindex >= 0 ? "" : " . No MC weight will be applied.")
              << std::endl;
    int nevents = -1;
    if (_vm.count("nevents")) {
      nevents = _vm["nevents"].as<int>();
    }

    int multiRuns = 1;
    if (_vm.count("multiRuns")) {
      multiRuns = _vm["multiRuns"].as<int>();
    }
    if (multiRuns != 1) {
      std::cout
          << "Using each event " << multiRuns
          << " times - note statistical errors are not correct at the moment"
          << std::endl;
    }

    _outputName = _outputName.substr(0, _outputName.rfind(".root"));
    _outputName = _outputName.substr(0, _outputName.rfind(".txt"));
    if (_outputName.size() != 0) {
      _mergedOutput = new TFile((_outputName + ".root").c_str(), "RECREATE");
      std::cout << "Output merged into: " << _outputName << ".[txt|root]"
                << std::endl;
    } else
      std::cout << "Output split per analysis" << std::endl;

    SetupAnalysList();

    if (_vm.count("help") || (_vm.count("input-files") == 0)) {
      std::cout << _desc << std::endl;
      return 1;
    }

    std::vector<std::string> inputFileNames;
    for (const auto &fileNames :
         _vm["input-files"].as<std::vector<std::string>>()) {
      splitCommaString(fileNames, inputFileNames);
    }
    std::cout << "Files to analyze: ";
    for (const auto &fileName : inputFileNames)
      std::cout << fileName << " ";
    std::cout << std::endl;

    TFile *fh = TFile::Open(inputFileNames[0].c_str());
    if (fh == 0) {
      std::cerr << "Failed to open the first file: " << inputFileNames[0]
                << std::endl;
      return 2;
    }
    Reader *reader = 0;
    for (auto *readerCand : *getReaderList()) {
      if (readerCand->isFileType(fh, _vm)) {
        reader = readerCand;
        break;
      }
    }
    if (!reader) {
      std::cerr << "Unknown input format in: " << inputFileNames[0]
                << std::endl;
      return 2;
    }
    reader->Init();
    reader->SetAnalyses(_analysisList);
    reader->SetMultiRuns(multiRuns);
    for (auto *smearer : *getSmearerList()) {
      if (_vm.count(smearer->getLongOption())) {
        std::vector<std::string> options;
        splitCommaString(_vm[smearer->getLongOption()].as<std::string>(),
                         options);
        smearer->init(options);
        reader->SetSmearing(smearer);
      }
    }
    for (auto *reweighter : *getReweighterList()) {
      if (_vm.count(reweighter->getLongOption())) {
        std::vector<std::string> options;
        splitCommaString(_vm[reweighter->getLongOption()].as<std::string>(),
                         options);
        reweighter->init(options);
        reader->AddReweighting(reweighter);
      }
    }
    reader->SetMCWeightIndex(mcwindex);
    reader->processFiles(inputFileNames, nevents);
    delete reader;

    std::ostream *oFile;
    if (_dumpText && _mergedOutput) {
      oFile = new std::ofstream(_outputName + ".txt");
    } else {
      oFile = new std::ostringstream;
    }
    bool first = true;
    for (const auto &analysis : _analysisList) {
      if (!_mergedOutput) {
        delete oFile;
        oFile = new std::ofstream(analysis->name() + ".txt");
      }
      analysis->getOutput()->saveRegions(*oFile, first);
      if (_mergedOutput)
        first = false;
    }
    delete oFile;

    return 0;
  }

protected:
  po::options_description _desc;
  std::string _outputName;
  TFile *_mergedOutput;
  bool _dumpText;
  bool _doNtuple;
  std::vector<AnalysisClass *> _analysisList;
  po::variables_map _vm;
};

#endif
