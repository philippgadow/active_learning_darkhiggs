#include "TopnessTool.h"

#include "TRandom3.h"
#include <algorithm>
#include <cassert>
#include <iostream>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

// top and W mass in MeV
const double m_top = 173.5e3;
const double m_W = 80.4e3;

// resolution parameters, also in MeV
const double a_W = 5e3;
const double a_t = 15e3;
const double a_CM = 1e6;

// squared constants
const double m_top_sqr = m_top * m_top;
const double m_W_sqr = m_W * m_W;
const double a_W_sqr = a_W * a_W;
const double a_t_sqr = a_t * a_t;
const double a_CM_sqr = a_CM * a_CM;

// Initialize object using measured kinematic variables.  two jets are
// required. if there is ==1 b-tagged jet in the event, three jets should be
// provided. the combinations (b,j1) and (b,j2) are evaluated in that case.
TopnessTool::TopnessTool(const TLorentzVector &p_l_, const TVector2 &p_miss_xy_,
                         const TLorentzVector &p_j1_,
                         const TLorentzVector &p_j2_)
    : p_nu_x(0), p_nu_y(0), p_nu_z(0), p_W_z(0), p_l(p_l_ * 1e3),
      p_miss_xy(p_miss_xy_ * 1e3), p_j1(p_j1_ * 1e3), p_j2(p_j2_ * 1e3),
      use_first_term(true), use_second_term(true), use_third_term(true),
      use_fourth_term(true) {}

TopnessTool::TopnessTool(const TLorentzVector &p_l_, const TVector2 &p_miss_xy_,
                         const TLorentzVector &p_j1_,
                         const TLorentzVector &p_j2_,
                         const TLorentzVector &p_j3_)
    : p_nu_x(0), p_nu_y(0), p_nu_z(0), p_W_z(0), p_l(p_l_ * 1e3),
      p_miss_xy(p_miss_xy_ * 1e3), p_j1(p_j1_ * 1e3), p_j2(p_j2_ * 1e3),
      p_j3(p_j3_ * 1e3), use_first_term(true), use_second_term(true),
      use_third_term(true), use_fourth_term(true) {}

double TopnessTool::minimize(int seed, int n_iterations, double tolerance,
                             double step) {
  TRandom3 r(seed);

  // Create ROOT minimizer
  ROOT::Math::Minimizer *minimum =
      ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "");

  // set tolerance , etc...
  minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(10000);       // for GSL
  minimum->SetTolerance(tolerance);
  minimum->SetPrintLevel(-2);

  // create function wrapper for minimizer
  ROOT::Math::Functor f(this, &TopnessTool::log_S, 4);
  minimum->SetFunction(f);

  // for 2 b-jets, there are two possible assignments.  if 3 jets were provided,
  // both combinations (j1,j2) and (j1,j3) will be tested, in both orders.
  int n_assignments = p_j3.E() > 0 ? 4 : 2;

  // run several minimizations, starting from different parameter values.
  // keep overall best result.
  double fmin_best = 999;
  TLorentzVector j1_tmp, j2_tmp;

  for (int i = 0; i < n_iterations; ++i) {
    // initial parameters are chosen random between +-4 TeV (parameters given to
    // amoeba are are in GeV to reduce numerical imprecision). eveywhere else
    // units are MeV.  the same value is used for all 4 parameters.
    std::vector<double> start(4, 4e3 * (2 * r.Rndm() - 1));

    for (int j = 0; j < n_assignments; ++j) {
      // Set the free variables to be minimized !
      minimum->SetVariable(0, "x", start[0], step);
      minimum->SetVariable(1, "y", start[1], step);
      minimum->SetVariable(2, "z", start[2], step);
      minimum->SetVariable(3, "a", start[3], step);
      minimum->Minimize();
      // keep track of overall best result
      if (minimum->MinValue() < fmin_best || i + j == 0) {
        fmin_best = minimum->MinValue();
        // save the best jets to temporary containers for setDelta
        j1_tmp = p_j1;
        j2_tmp = p_j2;
      }

      // next time, run minimization using different b-jet assignment
      std::swap(p_j1, p_j2);
      // if three jets are used, switch j2 and j3 every other time
      if (p_j3.E() > 0 && j % 2 == 1)
        std::swap(p_j2, p_j3);
    }
  }

  delete minimum;

  p_j1 = j1_tmp;
  p_j2 = j2_tmp;

  std::vector<double> params_fordelta(4);
  params_fordelta[0] = p_nu_x;
  params_fordelta[1] = p_nu_y;
  params_fordelta[2] = p_nu_z;
  params_fordelta[3] = p_W_z;

  setDelta(params_fordelta);

  return fmin_best;
}

static inline double sqr(double x) { return x * x; }

double TopnessTool::log_S(const double *params) const {

  //// rewritten without TLorentzVector ////
  double E_nu = sqrt(params[0] * params[0] + params[1] * params[1] +
                     params[2] * params[2]);
  double p1[4] = {p_l[0] + params[0], p_l[1] + params[1], p_l[2] + params[2],
                  p_l[3] + E_nu};
  double delta_m_W1_sqr =
      p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2] - p1[3] * p1[3] + m_W_sqr;
  p1[0] += p_j1[0];
  p1[1] += p_j1[1];
  p1[2] += p_j1[2];
  p1[3] += p_j1[3];

  double p2[4] = {-params[0] + p_miss_xy.Px(), -params[1] + p_miss_xy.Py(),
                  params[3], 0};
  p2[3] =
      p_j2[3] + sqrt(m_W_sqr + p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]);
  p2[0] += p_j2[0];
  p2[1] += p_j2[1];
  p2[2] += p_j2[2];

  double delta_m_side1_sqr =
      p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2] - p1[3] * p1[3] + m_top_sqr;
  double delta_m_side2_sqr =
      p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2] - p2[3] * p2[3] + m_top_sqr;
  double delta_m_ttbar_sqr = sqr(p1[0] + p2[0]) + sqr(p1[1] + p2[1]) +
                             sqr(p1[2] + p2[2]) - sqr(p1[3] + p2[3]) +
                             4 * m_top_sqr;

  double result2 =
      log((use_first_term ? sqr(delta_m_W1_sqr / a_W_sqr) : 0.0) +
          (use_second_term ? sqr(delta_m_side1_sqr / a_t_sqr) : 0.0) +
          (use_third_term ? sqr(delta_m_side2_sqr / a_t_sqr) : 0.0) +
          (use_fourth_term ? sqr(delta_m_ttbar_sqr / a_CM_sqr) : 0.0));

  return result2;
}

void TopnessTool::setDelta(const std::vector<double> &params) {

  //// rewritten without TLorentzVector ////
  double E_nu = sqrt(params[0] * params[0] + params[1] * params[1] +
                     params[2] * params[2]);
  double p1[4] = {p_l[0] + params[0], p_l[1] + params[1], p_l[2] + params[2],
                  p_l[3] + E_nu};
  // double delta_m_W1_sqr = p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] -
  // p1[3]*p1[3] + m_W_sqr;
  delta_W =
      p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2] - p1[3] * p1[3] + m_W_sqr;
  p1[0] += p_j1[0];
  p1[1] += p_j1[1];
  p1[2] += p_j1[2];
  p1[3] += p_j1[3];

  double p2[4] = {-params[0] + p_miss_xy.Px(), -params[1] + p_miss_xy.Py(),
                  params[3], 0};
  p2[3] =
      p_j2[3] + sqrt(m_W_sqr + p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]);
  p2[0] += p_j2[0];
  p2[1] += p_j2[1];
  p2[2] += p_j2[2];

  // double delta_m_side1_sqr = p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] -
  // p1[3]*p1[3] + m_top_sqr; double delta_m_side2_sqr = p2[0]*p2[0] +
  // p2[1]*p2[1] + p2[2]*p2[2] - p2[3]*p2[3] + m_top_sqr; double
  // delta_m_ttbar_sqr = sqr(p1[0] + p2[0]) + sqr(p1[1] + p2[1]) + sqr(p1[2] +
  // p2[2]) - sqr(p1[3] + p2[3]) + 4 * m_top_sqr;

  delta_t1l =
      p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2] - p1[3] * p1[3] + m_top_sqr;
  delta_t_Wmiss =
      p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2] - p2[3] * p2[3] + m_top_sqr;
  delta_t = sqr(p1[0] + p2[0]) + sqr(p1[1] + p2[1]) + sqr(p1[2] + p2[2]) -
            sqr(p1[3] + p2[3]) + 4 * m_top_sqr;

  // double result2 = log( ( use_first_term  ? sqr(delta_m_W1_sqr / a_W_sqr) :
  // 0.0 ) +
  //                       ( use_second_term ? sqr(delta_m_side1_sqr / a_t_sqr)
  //                       : 0.0 ) + ( use_third_term  ? sqr(delta_m_side2_sqr /
  //                       a_t_sqr) : 0.0 ) + ( use_fourth_term ?
  //                       sqr(delta_m_ttbar_sqr / a_CM_sqr): 0.0 ) );

  // return result2;
}

// The delta terms accessors
double TopnessTool::getDeltaW() { return delta_W; }

double TopnessTool::getDeltat1l() { return delta_t1l; }

double TopnessTool::getDeltatWmiss() { return delta_t_Wmiss; }

double TopnessTool::getDeltat() { return delta_t; }
