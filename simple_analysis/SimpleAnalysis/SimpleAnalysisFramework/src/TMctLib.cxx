#define tmctlib_cxx
#include "TMctLib.h"
#include <iostream>

#include <TLorentzVector.h>
#include <TVector2.h>

TMctLib::TMctLib() {}

TMctLib::~TMctLib() {}

double TMctLib::mctcorr(const TLorentzVector v1, const TLorentzVector v2,
                        const TLorentzVector vds, const TVector2 ptm,
                        const double ecm, const double mxlo) {
  double v1t[4] = {v1.E(), v1.Px(), v1.Py(), v1.Pz()};
  double v2t[4] = {v2.E(), v2.Px(), v2.Py(), v2.Pz()};
  double vdst[4] = {vds.E(), vds.Px(), vds.Py(), vds.Pz()};
  double ptmt[2] = {ptm.Px(), ptm.Py()};
  return mctcorr(v1t, v2t, vdst, ptmt, ecm, mxlo);
}

double TMctLib::mct(const TLorentzVector v1, const TLorentzVector v2) {
  double v1t[4] = {v1.E(), v1.Px(), v1.Py(), v1.Pz()};
  double v2t[4] = {v2.E(), v2.Px(), v2.Py(), v2.Pz()};
  return mct(v1t, v2t);
}

double TMctLib::mt2(const TLorentzVector v1, const TLorentzVector v2,
                    const TLorentzVector vds, const TVector2 ptm,
                    const double ecm, const double mxlo) {
  double v1t[4] = {v1.E(), v1.Px(), v1.Py(), v1.Pz()};
  double v2t[4] = {v2.E(), v2.Px(), v2.Py(), v2.Pz()};
  double vdst[4] = {vds.E(), vds.Px(), vds.Py(), vds.Pz()};
  double ptmt[2] = {ptm.Px(), ptm.Py()};
  return mt2(v1t, v2t, vdst, ptmt, ecm, mxlo);
}

double TMctLib::mt2neg(const TLorentzVector v1, const TLorentzVector v2,
                       const TVector2 ptm, const double mxlo) {
  double v1t[4] = {v1.E(), v1.Px(), v1.Py(), v1.Pz()};
  double v2t[4] = {v2.E(), v2.Px(), v2.Py(), v2.Pz()};
  double ptmt[2] = {ptm.Px(), ptm.Py()};
  return mt2neg(v1t, v2t, ptmt, mxlo);
}

double TMctLib::mcy(const TLorentzVector v1, const TLorentzVector v2,
                    const TLorentzVector vds, const TVector2 ptm) {
  double v1t[4] = {v1.E(), v1.Px(), v1.Py(), v1.Pz()};
  double v2t[4] = {v2.E(), v2.Px(), v2.Py(), v2.Pz()};
  double vdst[4] = {vds.E(), vds.Px(), vds.Py(), vds.Pz()};
  double ptmt[2] = {ptm.Px(), ptm.Py()};
  return mcy(v1t, v2t, vdst, ptmt);
}

double TMctLib::mcx(const TLorentzVector v1, const TLorentzVector v2,
                    const TLorentzVector vds, const TVector2 ptm) {
  double v1t[4] = {v1.E(), v1.Px(), v1.Py(), v1.Pz()};
  double v2t[4] = {v2.E(), v2.Px(), v2.Py(), v2.Pz()};
  double vdst[4] = {vds.E(), vds.Px(), vds.Py(), vds.Pz()};
  double ptmt[2] = {ptm.Px(), ptm.Py()};
  return mcx(v1t, v2t, vdst, ptmt);
}
