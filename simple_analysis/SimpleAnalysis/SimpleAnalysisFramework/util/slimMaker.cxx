#include "SimpleAnalysisFramework/NtupleMaker.h"
#include "execBase.h"

class slimExe : public execBase {
  using execBase::execBase;

  void AddOptions() {
    minElecPt = 3;
    minMuonPt = 3;
    minTauPt = 15;
    minPhotonPt = 15;
    minJetPt = 15;
    minFatJetPt = 100;
    minTrackJetPt = 10;
    _desc.add_options()("minElectronPt", po::value<double>(&minElecPt),
                        "Minimum electron pt [default: 3 GeV]")(
        "minMuonPt", po::value<double>(&minMuonPt),
        "Minimum muon pt [default: 3 GeV]")("minTauPt",
                                            po::value<double>(&minTauPt),
                                            "Minimum tau pt [default: 15 GeV]")(
        "minPhotonPt", po::value<double>(&minPhotonPt),
        "Minimum photon pt [default: 15 GeV]")(
        "minJetPt", po::value<double>(&minJetPt),
        "Minimum jet pt [default: 15 GeV]")(
        "minFatJetPt", po::value<double>(&minFatJetPt),
        "Minimum fat jet pt [default: 100 GeV]")(
        "minTrackJetPt", po::value<double>(&minTrackJetPt),
        "Minimum track jet pt [default: 10 GeV]");
  }

  void SetupAnalysList() {
    AnalysisClass *writer =
        new NtupleMaker(minElecPt, minMuonPt, minTauPt, minPhotonPt, minJetPt,
                        minFatJetPt, minTrackJetPt);
    OutputHandler *output = new OutputHandler(_mergedOutput, true);
    writer->setOutput(output);
    _analysisList.push_back(writer);
  }

private:
  double minElecPt;
  double minMuonPt;
  double minTauPt;
  double minPhotonPt;
  double minJetPt;
  double minFatJetPt;
  double minTrackJetPt;
};

int main(int argc, char **argv) {

  slimExe mainRunner("Slim xAOD truth to small ntuple", "slimNtuple", false);
  return mainRunner.run(argc, argv);

  return 0;
}
