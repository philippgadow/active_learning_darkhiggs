import uproot
from os.path import join


def relative_efficiency(cutflow, i):
    return float(cutflow[i]) / float(cutflow[i-1])

def base_efficiency(cutflow, i):
    return float(cutflow[i]) / float(cutflow[0])

def main():
    # get input files
    f_simpleanalysis = join('data', 'xampp_312363_simpleana_100022_monoSbb_zp2000_dm200_dh130', 'histograms.root')
    f_xampp_mc16a = join('data', 'xampp_312363_simpleana_100022_monoSbb_zp2000_dm200_dh130', 'user.changqia.22265276.XAMPP._000001.root')
    dsid = '312363'

    # get cutflow histograms
    with uproot.open(f_simpleanalysis) as f:
        cutflow_simpleanalysis_merged = f['Cutflow_merged']
        cutflow_simpleanalysis_resolved = f['Cutflow_resolved']
 
    with uproot.open(f_xampp_mc16a) as f:
        cutflow_xampp_merged = f['Histos_0L_SR_Merged_Nominal/InfoHistograms/DSID_312363_CutFlow']
        cutflow_xampp_resolved = f['Histos_0L_SR_Resolved_Nominal/InfoHistograms/DSID_312363_CutFlow']
    

    cuts_merged = [
        "Initial",
        "PassGRL",
        "passLArTile",
        "Trigger",
        "HasVtx",
        "BadJet",
        "CosmicMuon",
        "BadMuon",
        "PFlow Electron veto",
        "IsMETTrigPassed",
        "0 baseline electrons",
        "0 baseline muons",
        "Tau Veto",
        "Extended Tau Veto",
        "MetTST>500",
        ">=1 fat-jets",
        ">= 2 b-tagged track jets",
        "(mJ > 40 || mJ_corr > 40)",
        "(mJ> 50 && mJ < 270)",
        "N_associated_trkJets>=2",
        "(TrackJet_1passOR = true && TrackJet_2passOR = true)",
        "|DeltaPhiMin3|>20deg",
    ]

    cuts_resolved = [
        "Initial",
        "PassGRL",
        "passLArTile",
        "Trigger",
        "HasVtx",
        "BadJet",
        "CosmicMuon",
        "BadMuon",
        "PFlow Electron veto",
        "IsMETTrigPassed",
        "0 baseline electrons",
        "0 baseline muons",
        "Tau Veto",
        "Extended Tau Veto",
        "MetTST>150",
        "MetTST<=500",
        ">=2 jets",
        ">=2 b-tags",
        "(mjj > 40 || mjj_corr > 40)",
        "(mjj > 50 && mjj < 280)",
        "|DeltaPhiMin3|>20deg",
        "METSig>12",
        "(pt_jj > 100000 || pt_jj_corr > 100000)",
        "mT_METclosestBJet > 170000",
        "mT_METfurthestBJet > 200000",
        "<=4 jets",
    ]

    # convert histogram values to lists
    v_cutflow_simpleanalysis_merged = list(cutflow_simpleanalysis_merged.values())
    v_cutflow_xampp_merged = list(cutflow_xampp_merged.values())

    v_cutflow_simpleanalysis_resolved = list(cutflow_simpleanalysis_resolved.values())
    v_cutflow_xampp_resolved = list(cutflow_xampp_resolved.values())
   
    header = 'cut,events_xampp,events_simpleanalysis,efficiency_xampp,efficiency_simpleanalysis,relativeefficiency_xampp,relativeefficiency_simpleanalysis'

    output_merged = header + '\n'
    for i, cut in enumerate(v_cutflow_simpleanalysis_merged):
        if i == 0: continue
        cut_xp = v_cutflow_xampp_merged[i]
        cut_sa = v_cutflow_simpleanalysis_merged[i]
        output_merged += cuts_merged[i] + ',' + \
                         str(cut_xp) + ',' + str(cut_sa) + ',' + \
                         str(base_efficiency(v_cutflow_xampp_merged, i)) + ',' + str(base_efficiency(v_cutflow_simpleanalysis_merged, i)) + ',' + \
                         str(relative_efficiency(v_cutflow_xampp_merged, i)) + ',' + str(relative_efficiency(v_cutflow_simpleanalysis_merged, i)) + \
                         '\n'
    with open('cutflow_{dsid}_merged.csv'.format(dsid=dsid), 'w') as f:
        f.write(output_merged)

    output_resolved = header + '\n'
    for i, cut in enumerate(v_cutflow_simpleanalysis_resolved):
        if i == 0: continue
        cut_xp = v_cutflow_xampp_resolved[i]
        cut_sa = v_cutflow_simpleanalysis_resolved[i]
        output_resolved += cuts_resolved[i] + ',' + \
                           str(cut_xp) + ',' + str(cut_sa) + ',' + \
                           str(base_efficiency(v_cutflow_xampp_resolved, i)) + ',' + str(base_efficiency(v_cutflow_simpleanalysis_resolved, i)) + ',' + \
                           str(relative_efficiency(v_cutflow_xampp_resolved, i)) + ',' + str(relative_efficiency(v_cutflow_simpleanalysis_resolved, i)) + \
                           '\n'
    with open('cutflow_{dsid}_resolved.csv'.format(dsid=dsid), 'w') as f:
        f.write(output_resolved)

if __name__ == '__main__':
    main()
