#include "SimpleAnalysisFramework/METSignificanceCalculator.h"
#include "ObjectResolutions.h"

static void rotateXY(TMatrix &mat, TMatrix &mat_new, double phi) {
  double c = cos(phi);
  double s = sin(phi);
  double cc = c * c;
  double ss = s * s;
  double cs = c * s;

  mat_new(0, 0) =
      mat(0, 0) * cc + mat(1, 1) * ss - cs * (mat(1, 0) + mat(0, 1));
  mat_new(0, 1) =
      mat(0, 1) * cc - mat(1, 0) * ss + cs * (mat(0, 0) - mat(1, 1));
  mat_new(1, 0) =
      mat(1, 0) * cc - mat(0, 1) * ss + cs * (mat(0, 0) - mat(1, 1));
  mat_new(1, 1) =
      mat(0, 0) * ss + mat(1, 1) * cc + cs * (mat(1, 0) + mat(0, 1));
}

double calcMETSignificance(AnalysisEvent *event, bool applyOverlapRemoval) {
  auto electrons = event->getElectrons(10, 2.47);
  auto photons = event->getPhotons(10, 2.47);
  auto muons = event->getMuons(10, 2.5);
  auto jets = event->getJets(20., 4.5);
  auto taus = event->getTaus(10, 2.5);
  auto metVec = event->getMET();

  if (applyOverlapRemoval) { // This is a default overlap removal

    auto radiusCalcLepton = [](const AnalysisObject &lepton,
                               const AnalysisObject &) {
      return std::min(0.4, 0.04 + 10 / lepton.Pt());
    };

    auto muJetSpecial = [](const AnalysisObject &jet,
                           const AnalysisObject &muon) {
      if (jet.pass(NOT(BTag77MV2c10)) &&
          (jet.pass(LessThan3Tracks) || muon.Pt() / jet.Pt() > 0.5))
        return 0.2;
      else
        return 0.;
    };

    muons = AnalysisClass::overlapRemoval(muons, electrons, 0.01,
                                          NOT(MuCaloTaggedOnly));
    electrons = AnalysisClass::overlapRemoval(electrons, muons, 0.01);

    jets =
        AnalysisClass::overlapRemoval(jets, electrons, 0.2, NOT(BTag77MV2c10));
    electrons = AnalysisClass::overlapRemoval(electrons, jets, 0.2);

    jets = AnalysisClass::overlapRemoval(jets, muons, muJetSpecial,
                                         NOT(BTag77MV2c10));
    muons = AnalysisClass::overlapRemoval(muons, jets, 0.2);

    muons = AnalysisClass::overlapRemoval(muons, jets, radiusCalcLepton);
    electrons =
        AnalysisClass::overlapRemoval(electrons, jets, radiusCalcLepton);
  }

  return calcMETSignificance(electrons, photons, muons, jets, taus, metVec);
}

double calcMETSignificance(AnalysisObjects &electrons, AnalysisObjects &photons,
                           AnalysisObjects &muons, AnalysisObjects &jets,
                           AnalysisObjects &taus, AnalysisObject &metVec) {
  auto objects = electrons + photons + muons + jets + taus;

  auto softVec = metVec;
  double met = metVec.Et();

  TMatrix cov_sum(2, 2);

  TMatrix particle_u(2, 2), particle_u_rot(2, 2);
  for (auto obj : objects) {
    softVec += obj; // soft term is everything not included in hard objects
    double pt_reso = 0.0, phi_reso = 0.0;
    getObjectResolution(obj, pt_reso, phi_reso);
    particle_u(0, 0) = pow(pt_reso * obj.Pt(), 2);
    particle_u(1, 1) = pow(phi_reso * obj.Pt(), 2);
    rotateXY(particle_u, particle_u_rot, metVec.DeltaPhi(obj));
    cov_sum += particle_u_rot;
  }

  // add soft term resolution (fixed 10 GeV)
  particle_u(0, 0) = 10 * 10;
  particle_u(1, 1) = 10 * 10;
  rotateXY(particle_u, particle_u_rot, metVec.DeltaPhi(softVec));
  cov_sum += particle_u_rot;

  // calculate significance
  double varL = cov_sum(0, 0);
  double varT = cov_sum(1, 1);
  double covLT = cov_sum(0, 1);

  double significance = 0;
  double rho = 0;
  if (varL != 0) {
    rho = covLT / sqrt(varL * varT);
    if (fabs(rho) >= 0.9)
      rho = 0; // too large - ignore it
    significance = met / sqrt((varL * (1 - pow(rho, 2))));
  }
  return significance;
}
