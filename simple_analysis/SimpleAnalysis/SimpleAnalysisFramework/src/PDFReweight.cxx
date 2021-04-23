#include "PDFReweight.h"
#include "SimpleAnalysisFramework/AnalysisEvent.h"

#include <LHAPDF/LHAPDF.h>

DefineReweighter(PDFReweighter, "pdfReweight,p",
                 "PDF reweight to '<pdfName>[,energyIn,energyOut]'");

void PDFReweighter::init(std::vector<std::string> &options) {
  inputEnergy = 8.;
  outputEnergy = 8.;
  pdf = LHAPDF::mkPDF(
      options[0], 0); // FIXME: should one be able to specify iteration of pdf?
  std::cout << "Reweighting to " << options[0] << " PDF" << std::endl;
  if (options.size() == 3) {
    inputEnergy = stod(options[1]);
    outputEnergy = stod(options[2]);
    std::cout << "Reweighting from " << inputEnergy << " TeV to "
              << outputEnergy << " TeV" << std::endl;
  }
}

double PDFReweighter::reweightEvent(
    AnalysisEvent *event) { // FIXME: always reweights with ratio of new PDF
                            // ignoring generated PDF
  double eScale = inputEnergy / outputEnergy;
  double xpdf1_in = 1.0;
  double xpdf1_out = 1.0;
  double xpdf2_in = 1.0;
  double xpdf2_out = 1.0;
  xpdf1_in =
      pdf->xfxQ(event->getPDF_id1(), event->getPDF_x1(), event->getPDF_scale());
  xpdf2_in =
      pdf->xfxQ(event->getPDF_id2(), event->getPDF_x2(), event->getPDF_scale());
  xpdf1_out =
      pdf->xfxQ(event->getPDF_id1(), std::min(eScale * event->getPDF_x1(), 1.),
                event->getPDF_scale());
  xpdf2_out =
      pdf->xfxQ(event->getPDF_id2(), std::min(eScale * event->getPDF_x2(), 1.),
                event->getPDF_scale());

  double weight = 1;
  if (xpdf1_in > 0)
    weight *= xpdf1_out / xpdf1_in;
  if (xpdf2_in > 0)
    weight *= xpdf2_out / xpdf2_in;
  return weight;
}
