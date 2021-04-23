/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/*
 *  Author: Shunsuke Adachi, Univeristy of Tokyo
 *
 *  Calculate BDT score from weight file.
 */

#ifndef BDT_H
#define BDT_H

#include "TRandom3.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// for TMVA
#include "PathResolver/PathResolver.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

namespace BDT {

class BDTReader {

private:
  // MVA method name
  std::string m_methodname;

  // Variable definitions
  std::vector<std::string> m_variableDefs;

  // XML weight files
  std::string m_xmlfilename1; // for even eventnumber sample
  std::string m_xmlfilename2; // for odd  eventnumber sample
  bool m_useDoubleXML = false;
  // Reader
  TMVA::Reader *m_reader1; // for even eventnumber sample
  TMVA::Reader *m_reader2; // for odd  eventnumber sample
  TMVA::Reader *m_reader;  // reader used in each event

  // Event number to separate samples
  ULong64_t m_eventNumber = 0;
  TRandom3 m_rand;
  // Variable containers to be filled with value of variable
  std::map<std::string, Float_t *> *m_variables;

public:
  // Constructor
  BDTReader(std::string methodname, std::vector<std::string> variableDefs,
            std::string xmlfilename1, std::string xmlfilename2 = "")
      : m_methodname(methodname), m_variableDefs(variableDefs),
        m_xmlfilename1(xmlfilename1), m_xmlfilename2(xmlfilename2),
        m_reader1(0), m_reader2(0), m_reader(0) {

    // This loads the library
    TMVA::Tools::Instance();

    if (m_xmlfilename2 != "")
      m_useDoubleXML = true;
    m_rand = TRandom3();

    // Initialize reader
    m_reader1 = new TMVA::Reader("!Color:Silent");
    if (m_useDoubleXML)
      m_reader2 = new TMVA::Reader("!Color:Silent");

    // Initialize variable containers
    m_variables = new std::map<std::string, Float_t *>;

    // Set variables
    setVariables();

    // Book MVA method
    TString methodName = m_methodname + TString(" method");
    m_reader1->BookMVA(
        methodName,
        PathResolverFindCalibFile("SimpleAnalysisCodes/" + m_xmlfilename1));
    if (m_useDoubleXML) {
      m_reader2->BookMVA(
          methodName,
          PathResolverFindCalibFile("SimpleAnalysisCodes/" + m_xmlfilename2));
    }
  };

private:
  // Set variable's address to reader
  void setVariables(void) {
    for (auto label : m_variableDefs) {
      std::string var = label + " := " + label;
      Float_t *tmp = new Float_t();
      m_reader1->AddVariable(var, tmp);
      if (m_useDoubleXML)
        m_reader2->AddVariable(var, tmp);
      m_variables->insert(std::make_pair(label, tmp));
    }
    return;
  };

public:
  // Initililzation for each event (Need to be called before setValue() &
  // getBDT() )
  void initInEventLoop(ULong64_t eventNumber = 0) {
    m_eventNumber = eventNumber;
    // Select reader
    if (m_eventNumber == 0)
      m_eventNumber = m_rand.Integer(99999999);
    if (m_useDoubleXML) {
      if (m_eventNumber % 2 == 0) {
        m_reader = m_reader1;
      } else {
        m_reader = m_reader2;
      }
    } else
      m_reader = m_reader1;
    // Clean values in variable containers
    for (auto var : *m_variables) {
      *(var.second) = 0.;
    }
  };

  // Set value of variables
  void setValue(std::string label, Float_t value) {
    std::map<std::string, Float_t *>::iterator itr = m_variables->find(label);
    std::map<std::string, Float_t *>::iterator end = m_variables->end();
    if (itr == end) {
      std::cout << "  BDT::BDTReader::setValue() : Error!! Cannot find "
                   "variable label \""
                << label << "\"!" << std::endl;
      std::cout << "  -->> Skip set \"" << label << "\" = " << value << " !!"
                << std::endl;
      return;
    } else {
      *((*itr).second) = value;
    }
    return;
  };

  // Get BDT score
  double getBDT(void) {
    TString methodName = m_methodname + TString(" method");
    return m_reader->EvaluateMVA(methodName);
  };

}; // end of class BDTReader

} // end of namespace BDT

#endif
