#include "execBase.h"

class analysisExe : public execBase {
  using execBase::execBase;

  void AddOptions() {
    _desc.add_options()("analyses,a", po::value<std::string>(),
                        "Comma-separated list of analyses to run")(
        "listanalyses,l", "List available analyses and exit");
  }

  void SetupAnalysList() {
    if (_vm.count("listanalyses")) {
      std::cout << "Select analyses among:" << std::endl;
      for (auto *analysis : *getAnalysisList()) {
        std::cout << " " << analysis->name() << std::endl;
      }
      exit(0);
    }
    std::vector<std::string> analysisNames;
    if (_vm.count("analyses"))
      splitCommaString(_vm["analyses"].as<std::string>(), analysisNames);

    std::cout << "Analyses to run: ";
    for (const auto &name : analysisNames)
      std::cout << name << " ";
    if (analysisNames.size() == 0)
      std::cout << "all";
    std::cout << std::endl;

    bool selectAnalysis = analysisNames.size() > 0;
    bool selectByYear =
        analysisNames.size() > 0 && isdigit(analysisNames[0][0]) &&
        analysisNames[0].size() == 4; // can specify year 2016, etc.
    for (auto *analysis : *getAnalysisList()) {
      if (selectByYear) {
        bool inYear = false;
        for (const auto &year : analysisNames)
          if (analysis->name().find(year) != std::string::npos)
            inYear = true;
        if (!inYear)
          continue;
      } else if (selectAnalysis) {
        const auto namePtr = std::find(analysisNames.begin(),
                                       analysisNames.end(), analysis->name());
        if (namePtr == analysisNames.end())
          continue;
        analysisNames.erase(namePtr);
      }
      TFile *oRoot = _mergedOutput;
      if (!oRoot) {
        oRoot = new TFile((analysis->name() + ".root").c_str(), "RECREATE");
      }
      OutputHandler *output = new OutputHandler(oRoot, _doNtuple);
      if (_mergedOutput)
        output->title(analysis->name());
      analysis->setOutput(output);
      _analysisList.push_back(analysis);
    }
    if (analysisNames.size() != 0 && !selectByYear) {
      std::cerr << "Unknown analysis requested: ";
      for (const auto &name : analysisNames)
        std::cerr << name << " ";
      std::cerr << std::endl;
      exit(1);
    }

    std::sort(_analysisList.begin(), _analysisList.end(),
              [](AnalysisClass *a, AnalysisClass *b) {
                return b->name() > a->name();
              });
  }
};

int main(int argc, char **argv) {

  analysisExe mainRunner("Run one or more truth-level analyses", "", true);
  return mainRunner.run(argc, argv);
}
