#include "SUSYPicker.h"
#include "SimpleAnalysisFramework/AnalysisEvent.h"
#include <iostream>
#include <sstream>

static void splitDashString(const std::string& names,std::vector<int>& result) {
  std::stringstream ss(names);
  while( ss.good() ) {
    std::string substr;
    getline( ss, substr, '-' );
    result.push_back( stoi(substr) );
  }
}

DefineReweighter(SUSYPicker,"susyProcess,S",
		 "select only SUSY events with process number in certain ranges (example: '2-4,6,51-62'"); 

void SUSYPicker::init(std::vector<std::string>& ranges) {
  std::cout<<"Only considering events with susy points in the ranges:";
  for(const auto& range : ranges) {
    std::vector<int> points;
    splitDashString(range,points);
    if (points.size()==1) points.push_back(points[0]);
    if (points.size()!=2) {
      std::cerr<<"Not understanding range of SUSY processes: "<<range<<std::endl;
      exit(1);
    }
    _ranges.emplace_back(points[0],points[1]);
    std::cout<<" "<<points[0]<<"-"<<points[1];
  }
  std::cout<<std::endl;
}

double SUSYPicker::reweightEvent(AnalysisEvent *event) { 
  int process=event->getSUSYChannel();
  for( const auto& range : _ranges) {
    if (process>=range.first && process<=range.second) return 1.;
  }
  return 0;
}
