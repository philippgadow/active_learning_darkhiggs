# Simplified ATLAS SUSY analysis framework

Holds a collections of SUSY analyses. These can be run over samples in different
input formats:

 * DAOD_TRUTH      (TRUTH1 and TRUTH3 tested)
 * xAOD            (either truth-level and or reco-level - the latter with some constraints)
 * slimmed ntuples (Reduced ntuples produced from above input)

It provides the analysis acceptance per control and signal region as well as
optionally histograms or ntuples with event level objects. Smearing of the truth-level objects
can optionally be done using the latest upgrade performance smearing functions.

## Docs

Documentation is found at [https://simpleanalysis.docs.cern.ch](https://simpleanalysis.docs.cern.ch).

## Help

For more help, please contact ATLAS SUSY Working Group Conveners at atlas-phys-susy-conveners (AT) cern.ch.

# Acknowledgements

The ATLAS Collaboration thanks Christopher Rogan for granting explicit permission for [RestFrames](http://restframes.com/) to be used in this code.
