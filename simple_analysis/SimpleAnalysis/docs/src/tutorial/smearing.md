---
title: Smearing in SimpleAnalysis
---

## What is smearing?

In ATLAS, fully reconstructed MC samples undergo a detector simulation using either a full ATLAS simulation based on Geant4 or an ATLAS fast simulation. In many cases, however, even an ATLAS fast simulation is computationally too expensive, e.g. because the number of distinct signal models is simply too high for it to be computationally feasible. Smearing is a procedure that approximates a full detector simulation by using parameterisations for object reconstruction and identification efficiencies as well as object resolutions.

Smearing is one of the main features of SimpleAnalysis and is used to mimic the impact that detector effects have on truth distributions. Dedicated smearing is desirable in situations where one cannot afford fast detector simulations like AtlfastII (because they are not CPU-efficient enough). It uses the [UpgradePerformanceFunctions](https://acode-browser1.usatlas.bnl.gov/lxr/source/athena/PhysicsAnalysis/UpgradePhys/SmearingFunctions/UpgradePerformanceFunctions/?v=21.2) which rely on object efficiencies and resolutions parameterised in &eta; and &phi;.

!!! warning "Status of smearing functions"
    Most of the object resolutions and efficiencies in the UpgradePerformanceFunctions stem from early Run-2 conditions and need to be updated. However, the smearing in its current state has proven to be sufficiently good for many use-cases, even though it is currently a bit outdated.

## What objects are being smeared


### Electrons
For truth electrons, the identification efficiencies taken into account depend on &eta; and p<sub>T</sub> as well as on the reconstruction working point used. The probability of finding a fake electron in a truth jet is also estimated using a map in &eta; and p<sub>T</sub>. The range of the p<sub>T</sub> interpolation for the identification efficiencies and fake rates extends from 7 GeV to 120 GeV. If the p<sub>T</sub> of an electron is outside of that range, the efficiency and fake rate value from the respective bound of the corresponding &eta;-bin are used. The probability of an electron being falsely identified as a photon is approximated with separate fixed values for the barrel and end-cap regions. The transverse energy of the truth electron is subsequently smeared through a Gaussian term with a mean corresponding to the truth value and a standard deviation equal to the &eta;- and p<sub>T</sub>-dependent energy resolution.

### Photons
The identification efficiencies for photons are only parameterised in p<sub>T</sub> for the tight working point. Fake rates are emulated for both jets coming from the primary vertex and pileup jets. The p<sub>T</sub> of the photon is smeared with a Gaussian corresponding to the energy resolution.


### Muons
The identification efficiencies for truth muons are only parameterised in &eta; and the reconstruction working point used. As for truth electrons, the transverse energy of truth muons is smeared with a Gaussian term. The mean of the Gaussian is equal to the true transverse energy value and the resolution is approximated separately for the barrel and the end-cap regions.

### Taus
Tau identification efficiencies are parameterised in &eta; and p<sub>T</sub> for 1-prong and 3-prong decays consdiering taus with p<sub>T</sub> >20 GeV and &eta; <4, for loose, medium and tight working points. Fake tau identification rates (jets identified as taus) are also emulated considering jets with p<sub>T</sub> >20 GeV and &eta; <4.

### Jets
The energy of truth jets is smeared with a Gaussian term using a resolution approximation based on five &eta; bins ranging from &#124;&eta;&#124;=0 to &#124;&eta;&#124;=4.5. Only jets with a truth p<sub>T</sub> from 10 GeV to 1500 GeV are smeared. For truth jets with p<sub>T</sub> > 20 GeV, the flavour-tagging efficiencies are approximated using measured values from fully reconstructed ttbar MC samples. The approximated efficiencies are parameterised in &eta;, p<sub>T</sub> and the flavour of the jet but also in the specific flavour tagger and working point employed.

### Missing transverse energy
Finally, the smeared E<sub>T</sub><sup>miss</sup> is computed using the transverse momenta of the smeared truth objects and an approximation for the track soft terms (TST). The TST is approximated using results from Z&rarr;e<sup>+</sup>+e<sup>-</sup> events, allowing to derive a distribution of the mean TST projected in the direction longitudinal to the total transverse momentum of all hard objects p<sub>T</sub><sup>hard</sup> in an event. The measured resolution parallel and perpendicular to p<sub>T</sub><sup>hard</sup> is used to smear the nominal TST value.

## How to enable smearing

Enabling the smearing in SimpleAnalysis is rather straightforward. Just specify the smearing configuration using `-s layout=run2` (you can also set this to `run4` if you are doing HL-LHC studies) in the command line when running over input files:
```sh
simpleAnalysis -s layout=run2 -a MyAnalysisName C1N2_Wh_hbb_700p0_0p0_lep_DAOD_TRUTH3.root
```
