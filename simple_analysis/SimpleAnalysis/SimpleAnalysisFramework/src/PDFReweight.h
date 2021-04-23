#ifndef PDFREWEIGHT_H
#define PDFREWEIGHT_H

#include <string>
#include <vector>

#include "SimpleAnalysisFramework/Reweight.h"

class AnalysisEvent;
namespace LHAPDF {
class PDF;
}

class PDFReweighter : public Reweighter {
public:
  PDFReweighter();
  void init(std::vector<std::string> &options);
  double reweightEvent(AnalysisEvent *event);

private:
  double inputEnergy;
  double outputEnergy;
  LHAPDF::PDF *pdf;
};

#endif
