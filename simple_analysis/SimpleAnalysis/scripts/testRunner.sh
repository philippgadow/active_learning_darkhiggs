#!/bin/bash
echo "Testing on Gtt sample"
simpleAnalysis -o testGtt --nevents 5000 root://eosatlas.cern.ch//eos/atlas/user/a/aagaard/SimpleAnalysisTestSamples/slimGtt370179.root
rm -f Gtt.ref
xrdcp root://eosatlas.cern.ch//eos/atlas/user/a/aagaard/SimpleAnalysisTestSamples/slimGtt370179.reference Gtt.ref
echo "differences for Gtt"
diff Gtt.ref testGtt.txt
echo "Testing on slepton sample"
simpleAnalysis -o testSlepton root://eosatlas.cern.ch//eos/atlas/user/a/aagaard/SimpleAnalysisTestSamples/slimSlepton393162.root
rm -f slepton.ref
xrdcp root://eosatlas.cern.ch//eos/atlas/user/a/aagaard/SimpleAnalysisTestSamples/slimSlepton393162.reference slepton.ref
echo "differences for sleptons"
diff slepton.ref testSlepton.txt
