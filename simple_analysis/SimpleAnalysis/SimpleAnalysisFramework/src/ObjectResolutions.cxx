#include "ObjectResolutions.h"

#include "objres/JERData16EMTopo.h"
#include "objres/JERMC16EMTopo.h"

#include "objres/jetphi100.h"
#include "objres/jetphi20.h"
#include "objres/jetphi50.h"

#include "objres/r0_MS_MC_BARREL.h"
#include "objres/r1_ID_MC_BARREL.h"
#include "objres/r1_MS_MC_BARREL.h"
#include "objres/r2_ID_MC_BARREL.h"
#include "objres/r2_MS_MC_BARREL.h"

#include "objres/r0_MS_MC_ENDCAP.h"
#include "objres/r1_ID_MC_ENDCAP.h"
#include "objres/r1_MS_MC_ENDCAP.h"
#include "objres/r2_ID_MC_ENDCAP.h"
#include "objres/r2_MS_MC_ENDCAP.h"

#include "objres/electronConst90.h"
#include "objres/electronNoise90.h"
#include "objres/electronSampling90.h"

#include "objres/photonConst90.h"
#include "objres/photonNoise90.h"
#include "objres/photonSampling90.h"

#include "objres/tauRes1p0nBin0.h"
#include "objres/tauRes1p0nBin1.h"
#include "objres/tauRes1p0nBin2.h"
#include "objres/tauRes1p0nBin3.h"
#include "objres/tauRes1p0nBin4.h"

static double getHistValue(double *histdata, double x, double y,
                           bool interpolateX = false) {
  int nBinsX = int(histdata[0]);
  int nBinsY = int(histdata[1]);
  double *xbins = histdata + 2;
  double *ybins = histdata + 2 + nBinsX + 1;
  double *values = histdata + 2 + nBinsX + 1 + nBinsY + 1;
  x = std::min(x, xbins[nBinsX]);
  if (x < xbins[0] || x > xbins[nBinsX]) {
    std::cout << "x=" << x << " out of range [" << xbins[0] << ":"
              << xbins[nBinsX] << "]" << std::endl;
    return 0;
  }
  if (y < ybins[0] || y > ybins[nBinsY]) {
    std::cout << "y=" << y << " out of range [" << ybins[0] << ":"
              << ybins[nBinsY] << "]" << std::endl;
    return 0;
  }
  int xbin = 0;
  for (; xbin < nBinsX; xbin++)
    if (x < xbins[xbin + 1])
      break;
  int ybin = 0;
  for (; ybin < nBinsY; ybin++)
    if (y < ybins[ybin + 1])
      break;
  double value = values[xbin + ybin * nBinsX];
  if (interpolateX) {
    int nextbin = xbin + 1;
    if (nextbin < nBinsX && x < (xbins[nextbin] + xbins[xbin]) /
                                    2.) { // should interpolate to previous bin
      xbin = xbin - 1;
      value = values[xbin + ybin * nBinsX];
    }
    if (xbin >= 0 &&
        xbin < nBinsX - 1) { // simple linear interpolation like TGraph would do
      double value2 = values[xbin + 1 + ybin * nBinsX];
      double x1 = 0.5 * (xbins[xbin] + xbins[xbin + 1]);
      double x2 = 0.5 * (xbins[xbin + 1] + xbins[xbin + 2]);
      double frac = (x - x1) / (x2 - x1);
      value = (1 - frac) * value + frac * value2;
    }
  }
  return value;
}

static double getJetPtResolution(AnalysisObject &obj) {
  double et = std::min(3000., std::max(17., obj.Et())); // constraint 5-50 GeV
  double eta = std::min(4.5, fabs(obj.Eta()));
  double dataRes = getHistValue(hist_JERData16EMTopo, et, eta, true);
  double MCRes = getHistValue(hist_JERMC16EMTopo, et, eta, true);
  return std::max(dataRes, MCRes);
}

static double getJetPhiResolution(AnalysisObject &obj) {
  if (obj.Pt() < 50)
    return getHistValue(hist_jetphi20, obj.Eta(), obj.Phi());
  if (obj.Pt() < 100)
    return getHistValue(hist_jetphi50, obj.Eta(), obj.Phi());
  return getHistValue(hist_jetphi100, obj.Eta(), obj.Phi());
}

static double getMuonPtResolution(AnalysisObject &obj) {
  double *resHistsBarrel[] = {
      hist_r1_ID_MC_BARREL, hist_r2_ID_MC_BARREL, hist_r0_MS_MC_BARREL,
      hist_r1_MS_MC_BARREL, hist_r2_MS_MC_BARREL,
  };
  double *resHistsEndcap[] = {
      hist_r1_ID_MC_ENDCAP, hist_r2_ID_MC_ENDCAP, hist_r0_MS_MC_ENDCAP,
      hist_r1_MS_MC_ENDCAP, hist_r2_MS_MC_ENDCAP,
  };
  double **resHists = resHistsEndcap;
  if (fabs(obj.Eta()) < 1.05)
    resHists = resHistsBarrel;
  double pars[5];
  for (int ii = 0; ii < 5; ii++)
    pars[ii] = getHistValue(resHists[ii], obj.Phi(), obj.Eta());
  // we use same muon pt for everything, i.e. ignore energy loss before MS
  double IDResSq = pow(pars[0], 2) + pow(pars[1] * obj.Pt(), 2);
  double MSResSq =
      pow(pars[2] / obj.Pt(), 2) + pow(pars[3], 2) + pow(pars[4] * obj.Pt(), 2);
  return sqrt(IDResSq * MSResSq / (IDResSq + MSResSq));
}

static double getElectronEtResolution(AnalysisObject &obj) {
  const double rsampling =
      getHistValue(hist_electronSampling90, fabs(obj.Eta()), 0);
  const double rnoise = getHistValue(hist_electronNoise90, fabs(obj.Eta()), 0);
  const double rconst = getHistValue(hist_electronConst90, fabs(obj.Eta()), 0);
  const double energy = obj.Energy();
  double sigma2 = rsampling * rsampling / energy +
                  rnoise * rnoise / energy / energy + rconst * rconst;
  double et = std::min(50., std::max(5., obj.Et())); // constraint 5-50 GeV

  double pileupNoiseMeV = sqrt(32.) * (60. + 40. * log(et / 10.) / log(5.));
  double pileupSigma2 = pow(pileupNoiseMeV / 1000. / obj.Et(),
                            2); // not clear why Egamma uses Et and not E here?
  return sqrt(sigma2 + pileupSigma2);
}

static double getPhotonEtResolution(AnalysisObject &obj) {
  // photon resolution is taken from true unconverted photons
  const double rsampling =
      getHistValue(hist_photonSampling90, fabs(obj.Eta()), 0);
  const double rnoise = getHistValue(hist_photonNoise90, fabs(obj.Eta()), 0);
  const double rconst = getHistValue(hist_photonConst90, fabs(obj.Eta()), 0);
  const double energy = obj.Energy();
  double sigma2 = rsampling * rsampling / energy +
                  rnoise * rnoise / energy / energy + rconst * rconst;
  double et = std::min(50., std::max(5., obj.Et())); // constraint 5-50 GeV

  double pileupNoiseMeV = sqrt(32.) * (60. + 40. * log(et / 10.) / log(5.));
  double pileupSigma2 =
      pow(pileupNoiseMeV / 1000. / obj.Et(),
          2); // not clear why Egamma group uses Et and not E here?
  return sqrt(sigma2 + pileupSigma2);
}

static double getTauPtResolution(AnalysisObject &obj) {
  double eta = fabs(obj.Eta());
  double ptMeV =
      1000. *
      std::min(499., std::max(15., obj.Pt())); // only defined for 15-499 GeV
  double *hist = hist_tauRes1p0nBin0;
  if (eta > 0.3)
    hist = hist_tauRes1p0nBin1;
  if (eta > 0.8)
    hist = hist_tauRes1p0nBin2;
  if (eta > 1.3)
    hist = hist_tauRes1p0nBin3;
  if (eta > 1.6)
    hist = hist_tauRes1p0nBin4;
  double res = getHistValue(hist, ptMeV, 0, true);
  return res;
}

void getObjectResolution(AnalysisObject &obj, double &pt_reso,
                         double &phi_reso) {
  switch (obj.type()) {
  case ELECTRON:
    pt_reso = getElectronEtResolution(obj);
    phi_reso = 0.004;
    break;
  case MUON:
    pt_reso = getMuonPtResolution(obj);
    phi_reso = 0.001;
    break;
  case TAU:
    pt_reso = getTauPtResolution(obj);
    phi_reso = 0.01;
    break;
  case PHOTON:
    pt_reso = getPhotonEtResolution(obj);
    phi_reso = 0.004;
    break;
  case JET:
    pt_reso = getJetPtResolution(obj);
    phi_reso = getJetPhiResolution(obj);
    break;
  default:
    std::cout << "ERROR: getting resolution for unsupported object type: "
              << obj.type() << std::endl;
    pt_reso = 0.1;
    phi_reso = 0.001;
  }
}
