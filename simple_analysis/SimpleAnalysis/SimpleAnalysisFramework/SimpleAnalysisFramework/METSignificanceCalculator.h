#ifndef METSIGNIFICANCECALCULATOR_H
#define METSIGNIFICANCECALCULATOR_H

class AnalysisEvent;
class AnalysisObject;
class AnalysisObjects;

double calcMETSignificance(AnalysisEvent *event,
                           bool applyOverlapRemoval = false);
double calcMETSignificance(AnalysisObjects &electrons, AnalysisObjects &photons,
                           AnalysisObjects &muons, AnalysisObjects &jets,
                           AnalysisObjects &taus, AnalysisObject &metVec);

#endif
