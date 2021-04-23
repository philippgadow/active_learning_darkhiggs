#ifndef PDFVARIATIONS_H
#define PDFVARIATIONS_H

#include <string>
#include <vector>

#include "SimpleAnalysisFramework/Reweight.h"

class AnalysisEvent;
namespace LHAPDF {
class PDF;
}

class PDFVariations : public Reweighter {
private:
  double calcWeight(LHAPDF::PDF *pdf, AnalysisEvent *event, double pdf1,
                    double pdf2);
  std::string output(const std::vector<double> &accVariations);

public:
  PDFVariations();
  void init(std::vector<std::string> &options);
  double reweightEvent(AnalysisEvent *event);

private:
  LHAPDF::PDF *genpdf;
  std::vector<LHAPDF::PDF *> cteq66pdfs;
  std::vector<LHAPDF::PDF *> mstw2008pdfs;
};

#endif
