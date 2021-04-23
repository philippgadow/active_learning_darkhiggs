#include "SimpleAnalysisFramework/MVA.h"
#include "SimpleAnalysisFramework/AnalysisClass.h"

TMVAReader::TMVAReader(const std::string &name,
                       const std::vector<std::string> &variableDefs,
                       const std::string fname1, const std::string fname2)
    : MVA(name, variableDefs) {

  // This loads the library
  TMVA::Tools::Instance();

  // Initialize reader(s)
  m_bdt1 = new TMVA::Reader("!Color:Silent");
  if (fname2 != "")
    m_bdt2 = new TMVA::Reader("!Color:Silent");
  else
    m_bdt2 = 0;

  m_variables.resize(variableDefs.size(), 0);

  int idx = 0;
  for (auto label : variableDefs) {
    std::string var = label + " := " + label;
    m_bdt1->AddVariable(var, &m_variables[idx]);
    if (m_bdt2)
      m_bdt2->AddVariable(var, &m_variables[idx]);
    idx++;
  }

  m_bdt1->BookMVA(name, FindFile(fname1));
  if (m_bdt2)
    m_bdt2->BookMVA(name, FindFile(fname2));
}

double TMVAReader::evaluate(const std::vector<double> &values,
                            const std::string /* nodeName */) {
  TMVA::Reader *bdt = m_bdt1;
  if (m_bdt2 && ((m_eventNumber % 2) == 1))
    bdt = m_bdt2;

  if (values.size() != m_variables.size())
    throw std::runtime_error("Wrong number of variables into TMVAReader");
  for (size_t ii = 0; ii < values.size(); ii++)
    m_variables[ii] = values[ii];
  return bdt->EvaluateMVA(m_name);
}

MVAUtilsReader::MVAUtilsReader(const std::string &name,
                               const std::string fname1,
                               const std::string fname2)
    : MVA(name, {}) {
  TFile *f1 = TFile::Open(FindFile(fname1).c_str(), "READ");
  TTree *tree1 = nullptr;
  f1->GetObject(name.c_str(), tree1);
  if (tree1 == nullptr)
    throw std::runtime_error("Did not find MVA tree");
  m_bdt1 = new MVAUtils::BDT(tree1);
  m_bdt2 = nullptr;
  if (fname2 != "") {
    TFile *f2 = TFile::Open(FindFile(fname2).c_str(), "READ");
    TTree *tree2 = nullptr;
    f2->GetObject(name.c_str(), tree1);
    if (tree2 == nullptr)
      throw std::runtime_error("Did not find MVA tree");
    m_bdt2 = new MVAUtils::BDT(tree2);
  }
}

double MVAUtilsReader::evaluate(const std::vector<double> &values,
                                const std::string /* nodeName */) {
  MVAUtils::BDT *bdt = m_bdt1;
  if (m_bdt2 && ((m_eventNumber % 2) == 1))
    bdt = m_bdt2;
  std::vector<float> floatValues(values.begin(), values.end());

  return bdt->GetResponse(floatValues);
}

std::vector<double>
MVAUtilsReader::evaluateMulti(const std::vector<double> &values,
                              int numClasses) {
  MVAUtils::BDT *bdt = m_bdt1;
  if (m_bdt2 && ((m_eventNumber % 2) == 1))
    bdt = m_bdt2;

  std::vector<float> floatValues(values.begin(), values.end());
  auto results = bdt->GetMultiResponse(floatValues, numClasses);

  std::vector<double> doubleResult(results.begin(), results.end());
  return doubleResult;
}

ONNXReader::ONNXReader(const std::string &name, const std::string fname1,
                       const std::string fname2)
    : MVA(name, {}) {
  m_env = new Ort::Env(ORT_LOGGING_LEVEL_WARNING, name.c_str());
  Ort::SessionOptions sessionOptions;
  sessionOptions.SetIntraOpNumThreads(1);

  m_NN1 = new Ort::Session(*m_env, FindFile(fname1).c_str(), sessionOptions);
  m_NN2 = nullptr;
  if (fname2 != "")
    m_NN2 = new Ort::Session(*m_env, FindFile(fname2).c_str(), sessionOptions);

  Ort::AllocatorWithDefaultOptions allocator;
  m_input_node_names.reserve(m_NN1->GetInputCount());
  for (size_t i = 0; i < m_NN1->GetInputCount(); ++i) {
    m_input_node_names.push_back(m_NN1->GetInputName(i, allocator));
  }

  Ort::TypeInfo type_info = m_NN1->GetInputTypeInfo(
      0); // FIXME: Here we are assuming only one input node
  auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
  m_input_dimension = tensor_info.GetShape();

  m_output_node_names.reserve(m_NN1->GetOutputCount());
  for (size_t i = 0; i < m_NN1->GetOutputCount(); ++i)
    m_output_node_names.push_back(m_NN1->GetOutputName(i, allocator));
}

double ONNXReader::evaluate(const std::vector<double> &values,
                            const std::string /* nodeName */) {
  Ort::Session *nn = m_NN1;
  if (m_NN2 && ((m_eventNumber % 2) == 1))
    nn = m_NN2;
  std::vector<float> floatValues(values.begin(), values.end());
  auto memory_info =
      Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
      memory_info, floatValues.data(), floatValues.size(),
      m_input_dimension.data(), m_input_dimension.size());
  auto output_tensor =
      nn->Run(Ort::RunOptions{nullptr}, m_input_node_names.data(),
              &input_tensor, m_input_node_names.size(),
              m_output_node_names.data(), m_output_node_names.size());

  return *output_tensor.front().GetTensorMutableData<float>();
}
