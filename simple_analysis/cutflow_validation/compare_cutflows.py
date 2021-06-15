import uproot
from os.path import join
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt


def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('inputFile_xamppanalysis')
    parser.add_argument('inputFile_simpleanalysis')
    parser.add_argument('dsid')
    return parser

def relative_efficiency(cutflow, i):
    return float(cutflow[i]) / float(cutflow[i-1])

def base_efficiency(cutflow, i):
    return float(cutflow[i]) / float(cutflow[0])

def plotCutFlow(h_xampp, h_simpleanalysis, labels, output):
    fig, ax = plt.subplots()
    ax.plot(h_xampp, '+')
    ax.plot(h_simpleanalysis, '+')
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
        rotation_mode="anchor")
    ax.set_ylabel('Cut efficiency')
    fig.tight_layout()
    fig.savefig('{output}.png'.format(output=output))

def main():
    # get input files
    args = getArgumentParser().parse_args()
    f_simpleanalysis = args.inputFile_simpleanalysis
    f_xampp_mc16a = args.inputFile_xamppanalysis
    dsid = args.dsid

    # get cutflow histograms
    with uproot.open(f_simpleanalysis) as f:
        cutflow_simpleanalysis_merged = f['Cutflow_merged']
        cutflow_simpleanalysis_resolved = f['Cutflow_resolved']
 
    with uproot.open(f_xampp_mc16a) as f:
        cutflow_xampp_merged = f['Histos_0L_SR_Merged_Nominal/InfoHistograms/DSID_{dsid}_CutFlow'.format(dsid=dsid)]
        cutflow_xampp_resolved = f['Histos_0L_SR_Resolved_Nominal/InfoHistograms/DSID_{dsid}_CutFlow'.format(dsid=dsid)]
    

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
    efficiency_merged_xp = []
    efficiency_merged_sa = []
    for i, cut in enumerate(v_cutflow_simpleanalysis_merged):
        if i == 0: continue
        cut_xp = v_cutflow_xampp_merged[i]
        cut_sa = v_cutflow_simpleanalysis_merged[i]
        eff_xp = base_efficiency(v_cutflow_xampp_merged, i)
        eff_sa = base_efficiency(v_cutflow_simpleanalysis_merged, i)
        output_merged += cuts_merged[i] + ',' + \
                         str(cut_xp) + ',' + str(cut_sa) + ',' + \
                         str(eff_xp) + ',' + str(eff_sa) + ',' + \
                         str(relative_efficiency(v_cutflow_xampp_merged, i)) + ',' + str(relative_efficiency(v_cutflow_simpleanalysis_merged, i)) + \
                         '\n'
        efficiency_merged_xp.append(eff_xp)
        efficiency_merged_sa.append(eff_sa)
    with open('cutflow_{dsid}_merged.csv'.format(dsid=dsid), 'w') as f:
        f.write(output_merged)
    plotCutFlow(efficiency_merged_xp, efficiency_merged_sa, cuts_merged, 'cutflow_{dsid}_merged'.format(dsid=dsid))

    output_resolved = header + '\n'
    efficiency_resolved_xp = []
    efficiency_resolved_sa = []
    for i, cut in enumerate(v_cutflow_simpleanalysis_resolved):
        if i == 0: continue
        cut_xp = v_cutflow_xampp_resolved[i]
        cut_sa = v_cutflow_simpleanalysis_resolved[i]
        eff_xp = base_efficiency(v_cutflow_xampp_resolved, i)
        eff_sa = base_efficiency(v_cutflow_simpleanalysis_resolved, i)
        output_resolved += cuts_resolved[i] + ',' + \
                           str(cut_xp) + ',' + str(cut_sa) + ',' + \
                           str(base_efficiency(v_cutflow_xampp_resolved, i)) + ',' + str(base_efficiency(v_cutflow_simpleanalysis_resolved, i)) + ',' + \
                           str(relative_efficiency(v_cutflow_xampp_resolved, i)) + ',' + str(relative_efficiency(v_cutflow_simpleanalysis_resolved, i)) + \
                           '\n'
        efficiency_resolved_xp.append(eff_xp)
        efficiency_resolved_sa.append(eff_sa)
    with open('cutflow_{dsid}_resolved.csv'.format(dsid=dsid), 'w') as f:
        f.write(output_resolved)
    plotCutFlow(efficiency_resolved_xp, efficiency_resolved_sa, cuts_resolved, 'cutflow_{dsid}_resolved'.format(dsid=dsid))

if __name__ == '__main__':
    main()
