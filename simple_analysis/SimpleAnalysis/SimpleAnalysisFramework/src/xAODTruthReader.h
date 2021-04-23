#include "MCTruthClassifier/MCTruthClassifier.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "SimpleAnalysisFramework/AnalysisClass.h"
#include "SimpleAnalysisFramework/Reader.h"
#include "xAODRootAccess/TEvent.h"
#include <TFile.h>
#include <TTree.h>
#include <vector>

using std::vector;

class xAODTruthReader : public Reader {

public:
  xAODTruthReader();
  ~xAODTruthReader();
  void AddOptions(po::options_description &desc);
  bool isFileType(TFile *fh, po::variables_map &vm);
  void Init();

protected:
  bool processEvent(xAOD::TEvent *event, xAOD::TStore *store);
  void processFilesInternal(const std::vector<std::string> &inputNames,
                            unsigned int nevents);
  int getTruthOrigin(const xAOD::TruthParticle *part);
  int getTruthType(const xAOD::TruthParticle *part);
  int getMotherID(const xAOD::TruthParticle *part,
                  bool traverseProdChain = false);
  xAOD::TruthParticleContainer *
  findTruthParticles(xAOD::TStore *store,
                     const xAOD::TruthParticleContainer *truthparticles,
                     std::vector<int> pdgIds, int status = 1);
  xAOD::TruthParticleContainer *
  findTruthBSMParticles(xAOD::TStore *store,
                        const xAOD::TruthParticleContainer *truthparticles);
  const xAOD::TruthParticle *
  getFlavourSibling(const xAOD::TruthParticle *particle);

private:
  ST::SUSYObjDef_xAOD *_susytools;
  MCTruthClassifier *_mctool;
  xAOD::TEvent *_event;
  bool _useVisTau;
  bool _useTruthBSM;
};
