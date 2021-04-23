#include "TruthDecayer.h"
#include "SimpleAnalysisFramework/AnalysisClass.h"
#include "SimpleAnalysisFramework/TruthEvent.h"

// Class for decaying stable BSM particles *invisibly* with user specified
// lifetimes Note at the moment no actual decay products are added, only the
// decay vertex but the Etmiss is adjusted accordingly
// FIXME: should also allow for LLV neutralino2->chargino->neutralino1

DefineReweighter(TruthDecayer, "decay,D",
                 "Decay HS particles '<pdgId=Lifetime in "
                 "ns>[,seed=number][,status=status-code]'");

void TruthDecayer::init(std::vector<std::string> &options) {
  _random.SetSeed(12345);
  _status = 1;
  for (const auto &option : options) {
    std::string::size_type split = option.find("=");
    if (split == std::string::npos)
      throw std::runtime_error("Invalid option for truth decay");
    if (option.find("status=") == 0) {
      _status = stoi(option.substr(split + 1));
    } else if (option.find("seed=") == 0) {
      _random.SetSeed(stoi(option.substr(5)));
    } else {
      _decays[stoi(option.substr(0, split))] =
          stof(option.substr(split + 1)) * 1e-9;
    }
  }
  std::cout << "Decaying HS truth particles with status:" << _status
            << std::endl;
  for (const auto &decay : _decays)
    std::cout << " " << decay.first << " with lifetime " << decay.second * 1e9
              << " ns" << std::endl;
}

double TruthDecayer::reweightEvent(AnalysisEvent *event) {
  TruthEvent *truthEvent = dynamic_cast<TruthEvent *>(event);
  if (!truthEvent)
    throw std::runtime_error("Could not convert to truthEvent");
  for (auto &obj : truthEvent->_HSTruth) {
    if (obj.status() != _status)
      continue;
    auto decay = _decays.find(abs(obj.pdgId()));
    if (decay == _decays.end())
      continue;
    // std::cout<<"doing decay"<<obj.pdgId()<<" "<<decay->second<<std::endl;
    TVector3 flight =
        obj.Vect().Unit() * _random.Exp(TMath::C() * 1e3 * decay->second *
                                        obj.Gamma() * obj.Beta());
    obj.setDecayVtx(obj.prodVtx() + flight);
    obj.setId((obj.id() & 0x80ffffff) | StablePart);
    if (_status == 1)
      truthEvent->_met += obj; // since it disappears it becomes invisible -
                               // FIXME: this is wrong if it is LSP
  }
  return 1.; // we don't actually reweight the events
}
