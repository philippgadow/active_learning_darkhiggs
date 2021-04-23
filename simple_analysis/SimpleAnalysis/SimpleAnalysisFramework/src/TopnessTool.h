// use -*- C++ -*- mode
#ifndef TOPNESSTOOL_H
#define TOPNESSTOOL_H

#include <vector>

#include "TLorentzVector.h"
#include "TVector2.h"

class TopnessTool {
public:
  // Initialize object using measured kinematic variables.  two jets are
  // required. j1 is always assumed to be a b-jet. if there is exactly 1
  // b-tagged jet in the event, three jets should be provided. the combinations
  // (j1,j2) and (j1,j3) are evaluated in that case. units are assumed to be
  TopnessTool(const TLorentzVector &p_l, const TVector2 &p_miss_xy,
              const TLorentzVector &p_j1, const TLorentzVector &p_j2);
  TopnessTool(const TLorentzVector &p_l, const TVector2 &p_miss_xy,
              const TLorentzVector &p_j1, const TLorentzVector &p_j2,
              const TLorentzVector &p_j3);

  // Run minimization to find missing W and neutrino kinematic parameters.
  // All possible b-jet assignments are considered.
  // The best topness value found is returned.
  double minimize(int seed = 1, int n_iterations = 1, double tolerance = 4e-7,
                  double step = 20e3);

  // return value of topness for given parameters.
  double log_S(const double *params) const;

  // set delta terms to be extracted by the accessors
  // this is written for the study of the resolution constants
  void setDelta(const std::vector<double> &params);

  // delta terms accessors
  double getDeltaW();
  double getDeltat1l();
  double getDeltatWmiss();
  double getDeltat();

  // parameters are meaningful only after running a minimization.
  double p_nu_x;
  double p_nu_y;
  double p_nu_z;
  double p_W_z;

private:
  TLorentzVector p_l;
  TVector2 p_miss_xy;
  TLorentzVector p_j1;
  TLorentzVector p_j2;
  TLorentzVector p_j3;
  double delta_W;
  double delta_t1l;
  double delta_t_Wmiss;
  double delta_t;
  bool use_first_term;
  bool use_second_term;
  bool use_third_term;
  bool use_fourth_term;
};

#endif
