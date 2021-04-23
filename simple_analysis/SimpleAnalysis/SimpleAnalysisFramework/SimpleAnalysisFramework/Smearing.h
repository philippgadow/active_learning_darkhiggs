#ifndef SMEARING_H
#define SMEARING_H
#include <string>
#include <vector>

class AnalysisEvent;
class TruthEvent;
class Smearer;
std::vector<Smearer *> *getSmearerList(); // for automatically tracking which
                                          // smearings have been defined

class Smearer {
public:
  Smearer(const std::string &option, const std::string &desc)
      : _option(option), _desc(desc) {
    getSmearerList()->push_back(this);
  };
  Smearer(){};
  virtual void init(std::vector<std::string> &options) = 0;
  virtual TruthEvent *smearEvent(AnalysisEvent *event) = 0;
  virtual ~Smearer(){};
  const std::string &getOption() { return _option; };
  const std::string getLongOption() {
    return _option.substr(0, _option.find(","));
  };
  const std::string &getDesc() { return _desc; };

private:
  std::string _option;
  std::string _desc;
};

#define DefineSmearer(NAME, OPTION, DESC)                                      \
  NAME::NAME() : Smearer(OPTION, DESC) {}                                      \
  static const Smearer *NAME_instance __attribute__((used)) = new NAME()
#endif
