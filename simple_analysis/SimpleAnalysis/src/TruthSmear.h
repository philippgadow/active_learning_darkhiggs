#ifndef TRUTHSMEAR_H
#define TRUTHSMEAR_H

#include <vector>
#include <string>

#include "TRandom3.h"

#include "SimpleAnalysisFramework/Smearing.h"

class TruthEvent;
class AnalysisEvent;
namespace Upgrade {
  class UpgradePerformanceFunctions;
}

class TruthSmear : public Smearer
{
 public:
  TruthSmear();
  void init(std::vector<std::string>&);
  TruthEvent *smearEvent(AnalysisEvent*);

 private:
  void setTaggerHist(TFile *ff, int tagger, int operating_point, std::string dirName, bool trackJet = false);
  float getFlavourTagEfficiency(double ptGeV, float eta, int flavour, int tagger, bool trackJet = false);

  bool smearElectrons;
  bool smearMuons;
  bool smearTaus;
  bool smearPhotons;
  bool smearJets;
  bool smearMET;
  bool addPileupJets;
  bool useHGTD0;
  bool useHGTD1;
  bool useMuonHighEta;
  bool optPhotons;
  Upgrade::UpgradePerformanceFunctions *m_upgrade;
  Upgrade::UpgradePerformanceFunctions *m_upgradeMuonTight;
  Upgrade::UpgradePerformanceFunctions *m_upgradeMuonHighPt;
  TRandom3 m_random;
  std::map<int,std::array<TH2D*,4>> m_fEffs;
  std::map<int,std::array<TH2D*,4>> m_fEffs_trackjets;
};


#endif
