#include "SimpleAnalysisFramework/AnalysisEvent.h"
#include "LPXReweight.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#ifdef ROOTCORE_PACKAGE_LPXSignalReweightingTool

DefineReweighter(LPXReweight,"LPXReweight","arguments - use ""help"" for list");

void LPXReweight::init(std::vector<std::string>& options) {
  std::vector<std::string> ddoptions;
  for(const auto& option : options) { //boost::program_options expects "--..." options
    if (option.find("--")==0) ddoptions.push_back(option);
    else ddoptions.push_back("--"+option);
  }
  po::options_description desc("LPX Reweighting Options ('--' is optional)");
  desc.add_options()
    ("help", "list options and exit")
    ("channel", po::value<std::string>()->default_value("electron"), "dilepton finalstate")
    ("gmZMode", po::value<int>()->default_value(3), "gamma-Z mode")
    ("model", po::value<std::string>()->default_value("SSM"), "model to use")
    ("mass", po::value<float>()->default_value(5000), "Mass")
    ;
  po::variables_map vm;
  po::store( po::command_line_parser(ddoptions).options(desc).run(), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    exit(1);
  }
  std::cout<<"Doing Z' reweighting:"<<std::endl;
  std::cout<<" assuming input channel is: "<<vm["channel"].as<std::string>()<<std::endl;
  std::cout<<" gamma-Z mode is: "<<vm["gmZMode"].as<int>()<<std::endl;
  std::cout<<" reweighting to: "<<vm["model"].as<std::string>()<<" with mass "<<vm["mass"].as<float>()<<" GeV"<<std::endl;

  _mass = vm["mass"].as<float>();
  _model = vm["model"].as<std::string>();
  _ZPrimeRWTool = new ZPrimeSignalModule("ZPrimeRWTool");
  _ZPrimeRWTool->setChannel(vm["channel"].as<std::string>()).isSuccess();
  _ZPrimeRWTool->setgmZMode(vm["gmZMode"].as<int>()).isSuccess(); // e.g. pure Z' contribution only
  _ZPrimeRWTool->initialize().isSuccess();

  _kfactorTool = new LPXKfactorTool("LPXKfactorTool");
  _kfactorTool->setProperty("isMC15",true).isSuccess();
  _kfactorTool->initialize().isSuccess();

}



double LPXReweight::reweightEvent(AnalysisEvent */*event*/) { 
  if (!_ZPrimeRWTool->execute().isSuccess()) {
    std::cout<<"Error calculating LPX reweighting"<<std::endl;
    exit(1);
  }
  _ZPrimeRWTool->setZPrimeMass(_mass).isSuccess();
  if (_model=="SSM")
    _ZPrimeRWTool->setModelParametersSSM(_model).isSuccess();
  else
    _ZPrimeRWTool->setModelParametersE6(_model).isSuccess();

  double rw = _ZPrimeRWTool->getRWFactor();
  double kFactorWeight = 1.;
  if (!_kfactorTool->getKFactor(kFactorWeight).isSuccess()) {
    std::cout<<"Error calculating LPX kFactor"<<std::endl;
    exit(1);
  }
  return rw*kFactorWeight;
}
#endif
