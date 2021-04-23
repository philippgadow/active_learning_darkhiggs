#ifndef OUTPUTHANDLER_H
#define OUTPUTHANDLER_H
#include <TH1.h>
#include <TTree.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

class OutputHandler {
public:
  OutputHandler(TDirectory *oFile, bool doNtuple = false)
      : _eventWeight(0), _ntuple(0), _oFile(oFile), _doNtuple(doNtuple){};
  void init();
  int addEntry(const std::string &label);
  void addEntries(const std::vector<std::string> &labels);
  void addHistogram(const std::string &label, int bins, float min, float max);
  void addHistogram(const std::string &label, int bins, float *edges);
  void addHistogram(const std::string &label, std::vector<float> &edges);
  void addHistogram(const std::string &label,
                    std::vector<std::string> charBins);
  void addHistogram(const std::string &label, int binsX, float minX, float maxX,
                    int binsY, float minY, float maxY);

  void setEventWeight(double weight); // should only be called once per event
                                      // and before the pass method is called
  void
  adjustEventWeight(double weight); // for further adjustment of event weight
  void pass(const std::string &name, double weight = 1);
  void fillHistogram(const std::string &label, double x);
  void fillHistogram(const std::string &label, const char *bin);
  void fillHistogram(const std::string &label, double x, double y);
  void fillHistogramWeighted(const std::string &label, double x, double y,
                             double val);
  void ntupVar(const std::string &label, int value);
  void ntupVar(const std::string &label, float value);
  void ntupVar(const std::string &label, double value) {
    ntupVar(label, (float)value);
  };
  void ntupVar(const std::string &label, std::vector<int> &values);
  void ntupVar(const std::string &label, std::vector<float> &values);
  void ntupFill();
  void saveRegions(std::ostream &filehandle, bool header = false);
  void title(std::string title) { _title = title + "__"; };
  // the following are common for all instances of output handler
  static void addVariationNames(std::vector<std::string> &names);
  static void addVariationValues(std::vector<float> &values);
  static void resetVariationValues();
  static void setVariationOutput(
      std::function<std::string(const std::vector<double> &)> outFunc);

protected:
  void createNtuple();
  void pass(int num, double weight = 1);

private:
  double _eventWeight;
  std::map<std::string, int> label2idx;
  std::map<int, std::string> idx2label;
  std::vector<int> eventCounts;
  std::vector<double> weightedSums;
  std::vector<double> weightedSquaredSums;
  std::vector<std::vector<double> *> weightedVariationSums;
  std::string _title;
  std::unordered_map<std::string, TH1 *> histograms;
  std::unordered_map<std::string, int> ntupInt;
  std::unordered_map<std::string, float> ntupFloat;
  std::unordered_map<std::string, std::vector<int> *> ntupVectorInt;
  std::unordered_map<std::string, std::vector<float> *> ntupVectorFloat;
  std::vector<int *> ntupInts;
  std::vector<float *> ntupFloats;
  std::vector<std::vector<int> *> ntupVectorInts;
  std::vector<std::vector<float> *> ntupVectorFloats;
  TTree *_ntuple;
  TDirectory *_oFile;
  bool _doNtuple;

  static std::vector<std::string> variationNames;
  static std::vector<float> variationValues;
  static std::function<std::string(const std::vector<double> &)> _outFunc;
};
#endif
