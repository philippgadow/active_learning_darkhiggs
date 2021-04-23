#include <cmath>
#include <functional>

#include "PDFVariations.h"
#include "SimpleAnalysisFramework/AnalysisEvent.h"
#include "SimpleAnalysisFramework/OutputHandler.h"

#include <LHAPDF/LHAPDF.h>

DefineReweighter(PDFVariations, "pdfVariations,P",
                 "PDF reweight to '<generatedPdfName>'");

// Doing PDF variations according to UL method in arXiv:1206.2892

std::string PDFVariations::output(const std::vector<double> &accVariations) {
  if (accVariations.size() == 0) {
    return ",PDFMean,PDFError";
  }
  double cteq0 = accVariations[0];
  double cteqSum = 0;
  for (unsigned int ii = 1; ii <= cteq66pdfs.size() / 2; ii++)
    cteqSum += pow(accVariations.at(2 * ii - 1) - accVariations.at(2 * ii), 2);
  double ctequpdown = 0.5 * sqrt(cteqSum);

  unsigned int offset = cteq66pdfs.size();
  double mstw0 = accVariations[offset];
  double mstwUpSum = 0;
  double mstwDownSum = 0;
  for (unsigned int ii = 1; ii <= mstw2008pdfs.size() / 2; ii++) {
    double up = std::max(std::max(accVariations.at(offset + 2 * ii - 1) - mstw0,
                                  accVariations.at(offset + 2 * ii) - mstw0),
                         0.);
    mstwUpSum += up * up;
    double down =
        std::max(std::max(mstw0 - accVariations.at(offset + 2 * ii - 1),
                          mstw0 - accVariations.at(offset + 2 * ii)),
                 0.);
    mstwDownSum += down * down;
  }
  double mstwUp = sqrt(mstwUpSum);
  double mstwDown = sqrt(mstwDownSum);
  double U = std::max(cteq0 + ctequpdown, mstw0 + mstwUp);
  double L = std::min(cteq0 - ctequpdown, mstw0 - mstwDown);
  std::string result =
      "," + std::to_string((U + L) / 2.) + "," + std::to_string((U - L) / 2.);
  return result;
}

void PDFVariations::init(std::vector<std::string> &options) {
  genpdf = LHAPDF::mkPDF(options[0], 0);
  cteq66pdfs = LHAPDF::mkPDFs("cteq66");
  mstw2008pdfs = LHAPDF::mkPDFs("MSTW2008nlo90cl");
  std::vector<std::string> pdfNames;
  for (unsigned int ii = 0; ii < cteq66pdfs.size(); ii++)
    pdfNames.push_back("cteq66_" + std::to_string(ii));
  for (unsigned int ii = 0; ii < mstw2008pdfs.size(); ii++)
    pdfNames.push_back("mstw2008_" + std::to_string(ii));

  OutputHandler::addVariationNames(pdfNames);
  std::cout << "Doing " << (cteq66pdfs.size() + mstw2008pdfs.size())
            << " PDF variations assuming " << options[0] << " as input PDF"
            << std::endl;

  using std::placeholders::_1;
  std::function<std::string(const std::vector<double> &)> outputFunc =
      std::bind(&PDFVariations::output, this, _1);
  OutputHandler::setVariationOutput(outputFunc);
}

double PDFVariations::calcWeight(LHAPDF::PDF *pdf, AnalysisEvent *event,
                                 double xpdf1_in, double xpdf2_in) {
  double xpdf1_out =
      pdf->xfxQ(event->getPDF_id1(), event->getPDF_x1(), event->getPDF_scale());
  double xpdf2_out =
      pdf->xfxQ(event->getPDF_id2(), event->getPDF_x2(), event->getPDF_scale());
  double weight = 1;
  if (xpdf1_in > 0)
    weight *= xpdf1_out / xpdf1_in;
  if (xpdf2_in > 0)
    weight *= xpdf2_out / xpdf2_in;
  return weight;
}

double PDFVariations::reweightEvent(AnalysisEvent *event) {
  double xpdf1_in = 1.0;
  double xpdf2_in = 1.0;
  xpdf1_in = genpdf->xfxQ(event->getPDF_id1(), event->getPDF_x1(),
                          event->getPDF_scale());
  xpdf2_in = genpdf->xfxQ(event->getPDF_id2(), event->getPDF_x2(),
                          event->getPDF_scale());
  std::vector<float> pdfWeights;
  for (auto pdf : cteq66pdfs) {
    pdfWeights.push_back(calcWeight(pdf, event, xpdf1_in, xpdf2_in));
  }
  for (auto pdf : mstw2008pdfs) {
    pdfWeights.push_back(calcWeight(pdf, event, xpdf1_in, xpdf2_in));
  }
  OutputHandler::addVariationValues(pdfWeights);
  return 1;
}
