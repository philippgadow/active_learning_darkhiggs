#include <math.h>
#include <stdlib.h>

#include "SimpleAnalysisFramework/OutputHandler.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

void OutputHandler::init() {
  if (eventCounts.size() == 0)
    addEntry("All");
}

int OutputHandler::addEntry(const std::string &name) {
  int num = eventCounts.size();
  if (label2idx.find(name) != label2idx.end()) {
    std::cerr << "Duplicate signal region label: " << name << std::endl;
    exit(1);
  }
  if (name.find_first_not_of(
          "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890_") !=
      std::string::npos) {
    std::cerr << "Illegal signal region name: " << name << std::endl;
    std::cerr << "Signal region names should only have alphanumeric characters "
                 "and '_'\n"
              << std::endl;
    exit(1);
  }
  label2idx[name] = num;
  idx2label[num] = name;
  eventCounts.push_back(0);
  weightedSums.push_back(0);
  weightedSquaredSums.push_back(0);
  if (variationNames.size())
    weightedVariationSums.push_back(
        new std::vector<double>(variationNames.size()));

  return num;
}

void OutputHandler::addEntries(const std::vector<std::string> &labels) {
  for (const auto &label : labels)
    addEntry(label);
}

void OutputHandler::addHistogram(const std::string &label, int bins, float min,
                                 float max) {
  TH1 *hist = new TH1D((_title + label).c_str(), label.c_str(), bins, min, max);
  hist->Sumw2();
  hist->SetDirectory(_oFile);
  histograms[label] = hist;
}

void OutputHandler::addHistogram(const std::string &label, int bins,
                                 float *edges) {
  TH1 *hist = new TH1D((_title + label).c_str(), label.c_str(), bins, edges);
  hist->Sumw2();
  hist->SetDirectory(_oFile);
  histograms[label] = hist;
}

void OutputHandler::addHistogram(const std::string &label,
                                 std::vector<float> &edges) {
  float *my_edges = new float[edges.size()];
  for (unsigned int e = 0; e < edges.size(); ++e)
    my_edges[e] = edges[e];
  TH1 *hist = new TH1D((_title + label).c_str(), label.c_str(),
                       edges.size() - 1, my_edges);
  hist->Sumw2();
  hist->SetDirectory(_oFile);
  histograms[label] = hist;
}

void OutputHandler::addHistogram(const std::string &label,
                                 std::vector<std::string> charBins) {
  TH1 *hist = new TH1D((_title + label).c_str(), label.c_str(), charBins.size(),
                       0, charBins.size());
  for (unsigned int c = 0; c < charBins.size(); c++)
    hist->GetXaxis()->SetBinLabel(c + 1, charBins[c].c_str());
  hist->Sumw2();
  hist->SetDirectory(_oFile);
  histograms[label] = hist;
}

void OutputHandler::addHistogram(const std::string &label, int binsX,
                                 float minX, float maxX, int binsY, float minY,
                                 float maxY) {
  TH1 *hist = new TH2D((_title + label).c_str(), label.c_str(), binsX, minX,
                       maxX, binsY, minY, maxY);
  hist->Sumw2();
  hist->SetDirectory(_oFile);
  histograms[label] = hist;
}

void OutputHandler::pass(int num, double weight) {
  weight *= _eventWeight;
  if (idx2label.find(num) == idx2label.end()) {
    std::cerr << "Unknown signal region number: " << num << std::endl;
    exit(1);
  }
  if (num == 0) {
    ntupVar("eventWeight", weight);
    if (variationNames.size() != variationValues.size())
      throw(std::runtime_error(
          "Variation names and values not equal size - should not happen"));
    for (unsigned int ii = 0; ii < variationNames.size(); ii++) {
      ntupVar(variationNames[ii], variationValues[ii]);
    }
  } else {
    ntupVar(idx2label[num], weight);
  }
  eventCounts[num]++;
  weightedSums[num] += weight;
  weightedSquaredSums[num] += weight * weight;
  if (variationValues.size()) {
    for (unsigned int ii = 0; ii < variationValues.size(); ii++)
      weightedVariationSums[num]->at(ii) += weight * variationValues[ii];
  }
}

void OutputHandler::pass(const std::string &name, double weight) {
  if (label2idx.find(name) == label2idx.end()) {
    std::cerr << "Unknown signal region label: " << name << std::endl;
    exit(1);
  }
  pass(label2idx[name], weight);
}

void OutputHandler::fillHistogram(const std::string &label, double x) {
  if (histograms.find(label) == histograms.end()) {
    std::cerr << "Unknown histogram label: " << label << std::endl;
    exit(1);
  }
  histograms[label]->Fill(x, _eventWeight);
}

void OutputHandler::fillHistogram(const std::string &label, const char *bin) {
  if (histograms.find(label) == histograms.end()) {
    std::cerr << "Unknown histogram label: " << label << std::endl;
    exit(1);
  }
  histograms[label]->Fill(bin, _eventWeight);
}

void OutputHandler::fillHistogram(const std::string &label, double x,
                                  double y) {
  if (histograms.find(label) == histograms.end()) {
    std::cerr << "Unknown histogram label: " << label << std::endl;
    exit(1);
  }
  ((TH2 *)histograms[label])->Fill(x, y, _eventWeight);
}

void OutputHandler::fillHistogramWeighted(const std::string &label, double x,
                                          double y, double val) {
  if (histograms.find(label) == histograms.end()) {
    std::cerr << "Unknown histogram label: " << label << std::endl;
    exit(1);
  }
  ((TH2 *)histograms[label])->Fill(x, y, val);
}

void OutputHandler::setEventWeight(
    double weight) { // called once before analysis code is run
  for (auto &ntup : ntupInts)
    *ntup = 0;
  for (auto &ntup : ntupFloats)
    *ntup = 0;
  for (auto &ntupVec : ntupVectorFloats)
    ntupVec->clear();
  for (auto &ntupVec : ntupVectorInts)
    ntupVec->clear();

  _eventWeight = weight;
}

void OutputHandler::adjustEventWeight(double weight) { _eventWeight *= weight; }

void OutputHandler::createNtuple() {
  _ntuple = new TTree((_title + "ntuple").c_str(), "Simple Analysis ntuple");
  _ntuple->SetDirectory(_oFile);
}

void OutputHandler::ntupVar(const std::string &label, int value) {
  if (!_doNtuple)
    return;
  if (!_ntuple)
    createNtuple();
  if (ntupInt.find(label) == ntupInt.end()) {
    ntupInt[label] = 0;
    TBranch *branch = _ntuple->Branch(label.c_str(), &(ntupInt[label]),
                                      (label + "/I").c_str());
    ntupInts.push_back(&(ntupInt[label]));
    for (int ii = 0; ii < _ntuple->GetEntries(); ++ii)
      branch->Fill(); // backfill
  }
  ntupInt[label] = value;
}

void OutputHandler::ntupVar(const std::string &label, float value) {
  if (!_doNtuple)
    return;
  if (!_ntuple)
    createNtuple();
  if (ntupFloat.find(label) == ntupFloat.end()) {
    ntupFloat[label] = 0;
    TBranch *branch = _ntuple->Branch(label.c_str(), &(ntupFloat[label]),
                                      (label + "/F").c_str());
    ntupFloats.push_back(&(ntupFloat[label]));
    for (int ii = 0; ii < _ntuple->GetEntries(); ++ii)
      branch->Fill(); // backfill
  }
  ntupFloat[label] = value;
}

void OutputHandler::ntupVar(const std::string &label,
                            std::vector<float> &values) {
  if (!_doNtuple)
    return;
  if (!_ntuple)
    createNtuple();
  if (ntupVectorFloat.find(label) == ntupVectorFloat.end()) {
    ntupVectorFloat[label] = new std::vector<float>;
    TBranch *branch = _ntuple->Branch(label.c_str(), ntupVectorFloat[label]);
    ntupVectorFloats.push_back(ntupVectorFloat[label]);
    for (int ii = 0; ii < _ntuple->GetEntries(); ++ii)
      branch->Fill(); // backfill
  }
  std::vector<float> *dest = ntupVectorFloat[label];
  for (float value : values)
    dest->push_back(value);
}

void OutputHandler::ntupVar(const std::string &label,
                            std::vector<int> &values) {
  if (!_doNtuple)
    return;
  if (!_ntuple)
    createNtuple();
  if (ntupVectorInt.find(label) == ntupVectorInt.end()) {
    ntupVectorInt[label] = new std::vector<int>;
    TBranch *branch = _ntuple->Branch(label.c_str(), ntupVectorInt[label]);
    ntupVectorInts.push_back(ntupVectorInt[label]);
    for (int ii = 0; ii < _ntuple->GetEntries(); ++ii)
      branch->Fill(); // backfill
  }
  std::vector<int> *dest = ntupVectorInt[label];
  for (int value : values)
    dest->push_back(value);
}

void OutputHandler::ntupFill() { // called once after analysis code is run
  pass(0);                       // count all events
  if (_ntuple)
    _ntuple->Fill();
}

void OutputHandler::saveRegions(std::ostream &filehandle, bool header) {
  //  csv_ofile << "SR,Acceptance,Sum,WeightedSum,WeightedSquareSum"<<std::endl;
  std::vector<double> variationEff;
  if (header) {
    filehandle << "SR,events,acceptance,err";
    if (variationValues.size())
      filehandle << _outFunc(variationEff);
    filehandle << std::endl;
    // Not so obvious, but we only want to save to ROOT file once
    _oFile->Write();
    _oFile->Close();
  }
  for (const auto &idx : idx2label) {
    filehandle << _title + idx.second << "," << eventCounts[idx.first];
    if (idx.first == 0) {
      filehandle << "," << weightedSums[0] << "," << weightedSquaredSums[0];
    } else {
      double acc = weightedSums[idx.first] / weightedSums[0];
      double err = sqrt(weightedSquaredSums[idx.first]) / weightedSums[0];
      filehandle << "," << acc << "," << err;
    }
    if (variationValues.size()) {
      variationEff.clear();
      for (unsigned int ii = 0; ii < variationValues.size(); ii++)
        variationEff.push_back(weightedVariationSums[idx.first]->at(ii) /
                               weightedVariationSums[0]->at(ii));
      filehandle << _outFunc(variationEff);
    }
    filehandle << std::endl;
  }
}

std::vector<std::string> OutputHandler::variationNames;
std::vector<float> OutputHandler::variationValues;
std::function<std::string(const std::vector<double> &)>
    OutputHandler::_outFunc = 0;

void OutputHandler::addVariationNames(std::vector<std::string> &names) {
  for (const auto &name : names)
    variationNames.push_back(name);
}

void OutputHandler::addVariationValues(std::vector<float> &values) {
  for (const auto &value : values)
    variationValues.push_back(value);
}

void OutputHandler::resetVariationValues() { variationValues.clear(); }

void OutputHandler::setVariationOutput(
    std::function<std::string(const std::vector<double> &)> outFunc) {
  if (_outFunc)
    throw(std::runtime_error("Setting more than one output changer is not "
                             "supported - file bug report!"));
  _outFunc = outFunc;
}
