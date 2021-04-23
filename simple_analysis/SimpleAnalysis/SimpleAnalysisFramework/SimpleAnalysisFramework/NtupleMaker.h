#include "SimpleAnalysisFramework/AnalysisClass.h"

class NtupleMaker : public AnalysisClass {
public:
  NtupleMaker(double elecPt, double muonPt, double tauPt, double photonPt,
              double jetPt, double fatJetPt, double trackJetPt)
      : AnalysisClass("NtupleMaker"), _minElecPt(elecPt), _minMuonPt(muonPt),
        _minTauPt(tauPt), _minPhotonPt(photonPt), _minJetPt(jetPt),
        _minFatJetPt(fatJetPt), _minTrackJetPt(trackJetPt){};
  void Init(){};
  void ProcessEvent(AnalysisEvent *event);

private:
  double _minElecPt;
  double _minMuonPt;
  double _minTauPt;
  double _minPhotonPt;
  double _minJetPt;
  double _minFatJetPt;
  double _minTrackJetPt;
};
