---
title: Implementing ML methods
---

## Overview

It may be useful or essential to implement a ML method developed for your analysis in your SimpleAnalysis instance. The [ATLAS Machine-Learning forum](https://atlasml.web.cern.ch/atlasml/) gives some [recommendations](https://atlasml.web.cern.ch/atlasml/#porting-to-analysisreconstruction) on implementation in analysis code. Neural Networks are recommended through the ltwnn package and BDT are recommended via [MVAUtils](https://gitlab.cern.ch/atlas/athena/tree/21.2/Reconstruction/MVAUtils). These methods are lightweight and are easily integrated into ATLAS analysis code. The workflow proposed:

1. Train the model in your preferred method (python, notebooks).
2. Convert the trained model into MVAUtils/ltwnn friendly format.
3. Implement in analysis code via MVAUtils/ltwnn.


TMVA usage is also directly implemented in SimpleAnalysis, however in this section we'll illustrate the implementation of a BDT trained using XGBoost via MVAUtils as it's a more general use-case (it's also possible to convert TMVA models to MVAUtils). The model we'll implement is a simple binary model trained to distinguish a SUSY signal and \\(t \bar t \\) (multi-classification is also supported in MVAUtils).

## Converting the XGBoost model to ROOT file format (optional)

(Feel free to skip this step, the converted file `SATutorial_xgboost.root` can be found in this [cernbox link](https://cernbox.cern.ch/index.php/s/C2PjxKdxMssOjtr))

MVAUtils has some nice scripts for converting common model formats to ROOT tree format in  the `utils/` directory. Converting a binary XGBoost model with the script ```convertXGBoostToRootTree.py```
Running is as follows:

```sh
python convertXGBoostToRootTree.py <model_file> <output(.root)> --objective <objective> --test-file <test-file>
```

```<model_file>``` is the trained model, saved from training, ```<output>``` is the name of the output root file (must end in `.root`),  `<test-file>` is a file containing a numpy table of test input variables. `objective` is the objective function used in training. The test file is used for validating the output of the converted model is consistent with native python implementation.

The following lines are an example python script for producing a test file containing a numpy table:

```python
#!/usr/bin/python
import numpy as np
variables=["met","mbb","mct","mt"]
data_input = np.random.uniform(20, 500, size=(100, len(variables)))
print(data_input)
np.save('test-file.npy', data_input)
```

To run these scripts, you'll need to setup making the relevant packages. Many of the ML python packages are available on cvmfs using ```lsetup views``` (xgboost, pytorch etc):
```sh
setupATLAS
asetup 21.2.115,AnalysisBase #Just use the same AB for convenience
lsetup "views LCG_96b x86_64-centos7-gcc8-opt"#Replace with relevant system
```

You can now copy the model file `xgboost_model.model` can be found at this [cernbox link](https://cernbox.cern.ch/index.php/s/5ktpUsOqcq6CUgc)) to the current directory and convert:
```sh
cp /cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBase/21.2.115/InstallArea/x86_64-centos7-gcc8-opt/src/Reconstruction/MVAUtils/util/convertXGBoostToRootTree.py .
python convertXGBoostToRootTree.py xgboost_model.model SATutorial_xgboost.root --objective binary:logistic --test-file test-file.npy
```

## Implementation in SimpleAnalysis

If you skipped the previous step, go ahead and download the pre-converted BDT model ```SATutorial_xgboost.root``` via this [cernbox link](https://cernbox.cern.ch/index.php/s/toMh8zzVCZdkECH).

Move across the converted model into the `data` directory of your SimpleAnalysis clone:

```sh
mv SATutorial_xgboost.root $TUTORIAL_DIR/SimpleAnalysis/SimpleAnalysis/data/
```


Inside the `Init` function we'll add a new histogram to fill with BDT output and we'll create the MVA instance:

```cpp
void SATutorialCode::Init()
{
  // Define signal regions
  addRegions({"SR_h_Low_bin1", "SR_h_Low_bin2", "SR_h_Low_bin3"});
  addRegions({"SR_h_Med_bin1", "SR_h_Med_bin2", "SR_h_Med_bin3"});
  addRegions({"SR_h_High_bin1", "SR_h_High_bin2", "SR_h_High_bin3"});

  // Preselection
  addRegions({"preselection"});

  // Initialise some exemplary histograms
  addHistogram("hist_met",100,0,1000);
  addHistogram("hist_mct",100,0,1000);
//New code------------------------------------------------
  addHistogram("hist_BDT",100,0,1);

  //Initialise the BDT from the ROOT file in data
  TString input_file = PathResolverFindCalibFile("SimpleAnalysis/SATutorial_xgboost.root");
  TFile* f = TFile::Open(input_file);
  TTree* tree = nullptr;
  f->GetObject("xgboost", tree);//'xgboost' is the name of the tree in the ROOT file
  m_MVAUtilsBDT = new MVAUtils::BDT(tree);
//~New code-----------------------------------------------
}
```
The BDT will be initialised along with the histograms and regions specified. Now just to get some values out of it! In `ProcessEvent`:

```cpp
// Exemplary ntuple branches
ntupVar("met", met);
ntupVar("mbb", mbb);
//New code-------------------------------------------------
std::vector<float> BDT_input{met,mbb,mct,mt};//Create the input vector
float output = m_MVAUtilsBDT->GetClassification(BDT_input);//GetClassification for binary
fill("hist_BDT", 1-output);//1-since XGBoost uses opposite to normal convention
//~New code------------------------------------------------
```

Go ahead and re-compile the package and run the signal file input that we last ran:

```sh
cd $TUTORIAL_DIR/build
cmake ../SimpleAnalysis #Need to re-run cmake to update the data/ file
make
cd $TUTORIAL_DIR/run
simpleAnalysis -a <analysis_name> inputs/C1N2_Wh_hbb_700p0_0p0_lep_DAOD_TRUTH3.root -n
```

Now we can open up the output file and see the distribution of BDT scores for the signal, drawing the histogram `hist_BDT`:

![BDT output](images/BDT_output.png)

## Large statistics samples (optional)

If you like you can run the signal and background that this BDT is trained on. These are a SUSY signal sample and a \\(t \bar t \\) sample located the file ```DAOD_TRUTH3/```  in the  [cernbox link](https://cernbox.cern.ch/index.php/s/nKsOy7ZjipmdCAF). Be aware these are very large high-statistics samples. The inputs consist of many subfiles which can be run in SimpleAnalysis using wildcards:

```sh
simpleAnalysis -a SATutorialCode C1N2_Wh_hbb_300p0_150p0_lep_DAOD_TRUTH3/* -n -s layout=run2 --nevents 10000
mv SATutorialCode.root signal.root ##This will be overwritten when we run the next file
simpleAnalysis -a SATutorialCode mc15_13TeV.407346.PhPy8EG_A14_ttbarMET300_400_hdamp258p75_nonallhad.deriv.DAOD_TRUTH3.e6414_e5984_p3655/* -n -s layout=run2 --nevents 10000
```
The output scores are now available to see in the histogram (after having run with full statistics):

![signal ttbar](images/signal_ttbar.png)
