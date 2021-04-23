#ifndef LPXREWEIGHT_H
#define LPXREWEIGHT_H

#include <vector>
#include <string>
#include "SimpleAnalysisFramework/AnalysisClass.h"
#ifdef ROOTCORE_PACKAGE_LPXSignalReweightingTool

#include "LPXSignalReweightingTool/ZPrimeSignalModule.h"
#include "LPXKfactorTool/LPXKfactorTool.h"

#include "Reweight.h"

class AnalysisEvent;
class ZPrimeSignalModule;
class LPXKfactorTool;


class LPXReweight : public Reweighter
{
 public:
  LPXReweight();
  void init(std::vector<std::string>& options);
  double reweightEvent(AnalysisEvent *event);

 private:
  ZPrimeSignalModule* _ZPrimeRWTool;
  LPXKfactorTool* _kfactorTool;
  std::string _model;
  float _mass;
};


#endif
#endif
