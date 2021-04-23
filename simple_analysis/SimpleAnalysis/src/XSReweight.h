#ifndef XSREWEIGHT_H
#define XSREWEIGHT_H

#include <vector>
#include <string>

#include "SimpleAnalysisFramework/Reweight.h"
#include "SUSYTools/SUSYCrossSection.h"

class AnalysisEvent;

class XSReweighter : public Reweighter
{
 public:
  XSReweighter();
  void init(std::vector<std::string>& options);
  double reweightEvent(AnalysisEvent *event);

 private:
  SUSY::CrossSectionDB *xsecDB; 
  float lumi;
  int nEvents;
};


#endif
