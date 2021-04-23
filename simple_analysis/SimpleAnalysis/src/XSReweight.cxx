#include "XSReweight.h"
#include "SimpleAnalysisFramework/AnalysisClass.h"


DefineReweighter(XSReweighter,"xsReweight,x","XS reweight nEvents to LUMI '<LUMI (ifb)>,<nEvents>[,<XS DB file>]'"); 

void XSReweighter::init(std::vector<std::string>& options) {
  std::string file_db = "";
  lumi = stof(options[0]);
  nEvents = stoi(options[1]);
  if ( options.size()<3 || options[2]=="")
    file_db = "dev/PMGTools/PMGxsecDB_mc16.txt";
  else file_db = options[2];

  std::cout << "Loading XS DB from " << PathResolverFindCalibFile(file_db) << std::endl;
  xsecDB = new SUSY::CrossSectionDB(PathResolverFindCalibFile(file_db), true);

  if (!xsecDB)
  {
    std::cout << "Couldn't load cross section database from " << file_db << "!" << std::endl;
  } 
}

double XSReweighter::reweightEvent(AnalysisEvent *event) { 
  // *1000 to convert the xSec into units of fb
  return lumi * xsecDB->xsectTimesEff(event->getMCNumber()) * 1000 / nEvents;
}
