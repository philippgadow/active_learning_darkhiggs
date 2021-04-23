#include "SimpleAnalysisFramework/AnalysisRunner.h"

void AnalysisRunner::processEvent(TruthEvent *inputEvent, int eventNumber) {
  for (int run = 0; run < _runs; run++) {
    TruthEvent *event;
    if (_smear)
      event = _smear->smearEvent(inputEvent);
    else
      event = inputEvent;
    double weight = 1.;
    if (_mcwindex >= int(event->getMCWeights().size())) {
      throw std::runtime_error(
          "The specified MC weight index is out of range! ");
    }
    if (_mcwindex >= 0)
      weight = event->getMCWeights()[_mcwindex];
    event->sortObjects();
    OutputHandler::resetVariationValues();
    for (auto reweighter : _reweighter)
      weight *= reweighter->reweightEvent(event);
    if (weight != 0) { // skip events with 0 weight
      for (const auto &analysis : _analysisList) {
        analysis->getOutput()->setEventWeight(weight);
        analysis->getOutput()->ntupVar("Event", eventNumber);
        analysis->setEventNumber(eventNumber);
        analysis->ProcessEvent(event);
        analysis->getOutput()->ntupFill();
      }
    }
    if (event != inputEvent)
      delete event;
  }
}
