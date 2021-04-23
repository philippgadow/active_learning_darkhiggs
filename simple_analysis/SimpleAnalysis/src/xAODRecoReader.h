#include "xAODRootAccess/TEvent.h"
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "SimpleAnalysisFramework/AnalysisClass.h"
#include "SimpleAnalysisFramework/Reader.h"

using std::vector;

class AsgElectronLikelihoodTool;
class AsgPhotonIsEMSelector;
class ElectronPhotonShowerShapeFudgeTool;
class JetCalibrationTool;
class JetVertexTaggerTool;
class JERSmearingTool;
class JetCleaningTool;
class BTaggingSelectionTool;
namespace CP {
  class IEgammaCalibrationAndSmearingTool;
  class IMuonCalibrationAndSmearingTool;
  class IIsolationCorrectionTool;
  class IsolationSelectionTool;
  class MuonSelectionTool;
  class JetJvtEfficiency;
}
namespace TauAnalysisTools {
  class TauOverlappingElectronLLHDecorator;
  class TauSmearingTool;
  class TauSelectionTool;
}
class xAODRecoReader : public Reader {

 public:
  xAODRecoReader();
  ~xAODRecoReader();
  void AddOptions(po::options_description& desc);
  bool isFileType(TFile * fh, po::variables_map&  vm);
  void Init();

 protected:
  bool processEvent(xAOD::TEvent *event,xAOD::TStore *store);
  void processFilesInternal(const std::vector<std::string>& inputNames, unsigned int nevents);
  int getTruthOrigin(const xAOD::TruthParticle *part);
  int getTruthType(const xAOD::TruthParticle *part);
  xAOD::TruthParticleContainer* findTruthParticles(xAOD::TStore *store,
                           const xAOD::TruthParticleContainer* truthparticles,
                           std::vector<int> pdgIds, int status=1);

 private:
  ST::SUSYObjDef_xAOD *_susytools;
  AsgElectronLikelihoodTool* _ElecVeryLooseLH;
  AsgElectronLikelihoodTool* _ElecLooseLH;
  AsgElectronLikelihoodTool* _ElecLooseBLLH;
  AsgElectronLikelihoodTool* _ElecMediumLH;
  AsgElectronLikelihoodTool* _ElecTightLH;
  AsgPhotonIsEMSelector*     _photonLoose;
  AsgPhotonIsEMSelector*     _photonTight;
  ElectronPhotonShowerShapeFudgeTool* _electronPhotonShowerShapeFudgeTool;
  CP::IsolationSelectionTool* _isoTrackLoose;
  CP::IsolationSelectionTool* _isoLoose;
  CP::IsolationSelectionTool* _isoGradient;
  CP::IsolationSelectionTool* _isoGradientLoose;
  CP::IsolationSelectionTool* _isoFixedCutLoose;
  CP::IsolationSelectionTool* _isoFixedCutTight;
  CP::IsolationSelectionTool* _isoFixedCutTightTrackOnly;
  CP::IEgammaCalibrationAndSmearingTool* _ElecCalibTool;
  CP::IIsolationCorrectionTool* _ElecIsoCorrTool; //!
  CP::IMuonCalibrationAndSmearingTool* _MuonCalibTool;
  CP::MuonSelectionTool*      _MuonSelectionTool;
  JetCalibrationTool*         _fatJetCalibrationTool;
  JetCalibrationTool*         _JetCalibrationTool;
  JetVertexTaggerTool*        _JetVertexTaggerTool;
  CP::JetJvtEfficiency*       _JetJvtEfficiencyTool;
  JERSmearingTool*            _JERSmearingTool;
  JetCleaningTool*            _JetCleaningLooseTool;
  JetCleaningTool*            _JetCleaningTightTool;
  BTaggingSelectionTool*     _Btagging85Tool;
  BTaggingSelectionTool*     _Btagging77Tool;
  BTaggingSelectionTool*     _Btagging70Tool;
  BTaggingSelectionTool*     _Btagging60Tool;
  TauAnalysisTools::TauOverlappingElectronLLHDecorator * _tauElORdecorator;
  TauAnalysisTools::TauSmearingTool* _tauSmearingTool;
  TauAnalysisTools::TauSelectionTool* _TauBDTLooseTool;
  TauAnalysisTools::TauSelectionTool* _TauBDTMediumTool;
  TauAnalysisTools::TauSelectionTool* _TauBDTTightTool;
  xAOD::TEvent* _event;
};
