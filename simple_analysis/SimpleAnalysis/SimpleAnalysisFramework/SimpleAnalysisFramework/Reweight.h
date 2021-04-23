#ifndef REWEIGHT_H
#define REWEIGHT_H
#include <string>
#include <vector>

class AnalysisEvent;
class Reweighter;
std::vector<Reweighter *> *
getReweighterList(); // for automatically tracking which reweightings have been
                     // defined

class Reweighter {
public:
  Reweighter(const std::string &option, const std::string &desc)
      : _option(option), _desc(desc) {
    getReweighterList()->push_back(this);
  };
  Reweighter(){};
  virtual void init(std::vector<std::string> &options) = 0;
  virtual double reweightEvent(AnalysisEvent *event) = 0;
  virtual ~Reweighter(){};
  const std::string &getOption() { return _option; };
  const std::string getLongOption() {
    return _option.substr(0, _option.find(","));
  };
  const std::string &getDesc() { return _desc; };

private:
  std::string _option;
  std::string _desc;
};

#define DefineReweighter(NAME, OPTION, DESC)                                   \
  NAME::NAME() : Reweighter(OPTION, DESC) {}                                   \
  static const Reweighter *NAME##_instance __attribute__((used)) = new NAME()

#endif
