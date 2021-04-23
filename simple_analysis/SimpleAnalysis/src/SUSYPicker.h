#ifndef SUSYPICKER_H
#define SUSYPICKER_H

#include <vector>
#include <string>

#include "SimpleAnalysisFramework/Reweight.h"

class AnalysisEvent;

class SUSYPicker : public Reweighter
{
 public:
  SUSYPicker();
  void init(std::vector<std::string>& ranges);
  double reweightEvent(AnalysisEvent *event);

 private:
  std::vector<std::pair<int,int>> _ranges;
};


#endif
