#include "SimpleAnalysisFramework/AnalysisObject.h"
#include "SimpleAnalysisFramework/AnalysisClass.h"

AnalysisObject operator+(const AnalysisObject &lhs, const AnalysisObject &rhs) {
  const TLorentzVector &tlhs = lhs;
  return AnalysisObject(tlhs + rhs, lhs.charge() + rhs.charge(), 0, COMBINED, 0,
                        0);
}

AnalysisObjects operator+(const AnalysisObjects &lhs,
                          const AnalysisObjects &rhs) {
  AnalysisObjects combined;
  for (const auto &cand : lhs)
    combined.push_back(cand);
  for (const auto &cand : rhs)
    combined.push_back(cand);
  AnalysisClass::sortObjectsByPt(combined);
  return combined;
}
