#ifndef MVA_H
#define MVA_H

#include "TRandom3.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "MVAUtils/BDT.h"
#include "TFile.h"
#include "TTree.h"

#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/NNLayerConfig.hh"
#include "lwtnn/lightweight_network_config.hh"
#include "lwtnn/parse_json.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <onnxruntime/core/session/onnxruntime_cxx_api.h>
#pragma GCC diagnostic pop

class MVA {
public:
  MVA(const std::string &name, const std::vector<std::string> variableDefs = {}
      // specifies the order of variables in the input vector
      )
      : m_name(name), m_variableDefs(variableDefs){};

  void setEventNumber(int eventNumber) { m_eventNumber = eventNumber; };
  virtual double evaluate(const std::vector<double> &values,
                          const std::string nodeName = "") = 0;
  virtual std::vector<double>
  evaluateMulti(const std::vector<double> & /*values */, int /* numClasses */) {
    throw std::runtime_error("multi output not supported for this type of MVA");
  };

protected:
  std::string m_name;
  std::vector<std::string> m_variableDefs;
  int m_eventNumber;
};

class TMVAReader : public MVA {
public:
  TMVAReader(const std::string &name,
             const std::vector<std::string> &variableDefs,
             const std::string fname1, const std::string fname2 = "");
  virtual double evaluate(const std::vector<double> &values,
                          const std::string nodeName = "");
  ~TMVAReader() {
    delete m_bdt1;
    delete m_bdt2;
  };

private:
  TMVA::Reader *m_bdt1; // for even eventnumber sample
  TMVA::Reader *m_bdt2; // for odd  eventnumber sample
  std::vector<Float_t> m_variables;
};

class MVAUtilsReader : public MVA {
public:
  MVAUtilsReader(const std::string &name, const std::string fname1,
                 const std::string fname2 = "");
  virtual double evaluate(const std::vector<double> &values,
                          const std::string nodeName = "");
  virtual std::vector<double> evaluateMulti(const std::vector<double> &values,
                                            int numClasses);
  ~MVAUtilsReader() {
    delete m_bdt1;
    delete m_bdt2;
  };

private:
  MVAUtils::BDT *m_bdt1;
  MVAUtils::BDT *m_bdt2;
};

class LWTNNReader
    : public MVA { // FIXME: This class is not implemented/supported yet
public:
  LWTNNReader(const std::string &name, const std::string fname1,
              const std::string fname2 = "",
              const std::vector<std::string> variableDefs = {});
  virtual double evaluate(const std::vector<double> &values,
                          const std::string nodeName = "");
  ~LWTNNReader() {
    delete m_NN1;
    delete m_NN2;
  };

private:
  lwt::LightweightGraph *m_NN1;
  lwt::LightweightGraph *m_NN2;
};

class ONNXReader : public MVA {
public:
  ONNXReader(const std::string &name, const std::string fname1,
             const std::string fname2 = "");
  virtual double evaluate(const std::vector<double> &values,
                          const std::string nodeName = "");
  ~ONNXReader() {
    delete m_NN1;
    delete m_NN2;
    delete m_env;
  };

private:
  Ort::Env *m_env;
  Ort::Session *m_NN1;
  Ort::Session *m_NN2;
  std::vector<int64_t> m_input_dimension;
  std::vector<const char *> m_input_node_names;
  std::vector<const char *> m_output_node_names;
};

#endif
