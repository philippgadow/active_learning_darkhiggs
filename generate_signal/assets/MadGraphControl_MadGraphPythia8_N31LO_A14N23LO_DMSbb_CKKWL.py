import MadGraphControl.MadGraphUtils
from MadGraphControl.MadGraphUtils import *


#####################
# Settings
#####################
# safe factor to ensure a sufficient number of event has been generated
safefactor=2.5

# determine whether MadGraph reweight module should be run
try:
    reweight
except NameError:
    from AthenaCommon import Logging
    locallog = Logging.logging.getLogger('monoSbb')
    locallog.info("reweight not set, defaulting to 'False'")
    reweight = False


evgenConfig.contact = ["Paul Philipp Gadow <pgadow@cern.ch>"]
evgenConfig.generators = ["MadGraph", "Pythia8", "EvtGen"]

# Get parameters from physics short name
from MadGraphControl.MadGraphUtilsHelpers import get_physics_short
phys_short = get_physics_short()
mzp = int(phys_short.split('_zp')[1].split('_')[0])
mdm = int(phys_short.split('_dm')[1].split('_')[0])
mhs = int(phys_short.split('_dh')[1].split('_')[0])


###########
# Process
###########
# form full process string and set up directory
process = """
import model DarkHiggs2MDM
generate p p > zp > n1 n1 hs QED<=2, (hs > b b~) @0
add process p p > zp > n1 n1 hs j QED<=2, (hs > b b~) @1
output -f
"""
process_dir = new_process(process)


########
# PDF
########
MadGraphControl.MadGraphUtils.MADGRAPH_PDFSETTINGS={
    'central_pdf':315200,     # NNPDF31_lo_as_0130
    'pdf_variations':315200,  # NNPDF31_lo_as_0130
    'alternative_pdfs':[13200,25000,27000],  # CT14lo, MMHT2014lo68cl, MSHT20lo_as130
    'scale_variations':[0.5,1.,2.],
}


###########
# Run card
###########
# determine ktdurham cut from dark Higgs mass
# (ktdurham cut sets scale at which event description is split between parton shower and matrix element) 
try:
    ktdurham = int(mhs / 4)
    assert ktdurham > 40
except AssertionError:
    ktdurham = 40

# set settings for run card
nevents = runArgs.maxEvents*safefactor if runArgs.maxEvents>0 else safefactor*evgenConfig.nEventsPerJob
settings = {'lhe_version':'3.0',
            'cut_decays': 'F',
            'event_norm': 'sum',
            'drjj': "0.0",         # required for CKKW-L jet matching
            'ickkw': 0,            # required for CKKW-L jet matching
            'ktdurham': ktdurham,  # required for CKKW-L jet matching
            'dparameter': "0.4",   # required for CKKW-L jet matching
            'xqcut': "0.0",        # required for CKKW-L jet matching
            'nevents': nevents
          }
modify_run_card(process_dir=process_dir,runArgs=runArgs,settings=settings)


##################
# Parameter card
##################
# set parameters for parameter card
params = {}
# mass
params['mass'] = {'54': mhs, '55': mzp, '1000022': mdm}
# couplings
params['frblock'] = {'1': gq , '2': gx , '3': th}
# decay width
params['decay'] = {'54':"AUTO" , '55':"AUTO"}
modify_param_card(process_dir=process_dir, params=params)


##################
# Reweight card
##################
if reweight:
    gx_scan = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
               1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
               2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
               3.1, 3.2, 3.3, 3.4, 3.5] # gchi < (4 pi)^0.5 \approx 3.55 (perturbativity bound)
    reweightCommand=""
    for i_gx in gx_scan:
        reweightCommand += "launch --rwgt_name=rwgt_gx_{gx_str}\n".format(gx_str=str(i_gx).replace('.', 'p'))
        reweightCommand += "set frblock 2 {gx}\n\n".format(gx=i_gx)
    rcard = open(os.path.join(process_dir,'Cards', 'reweight_card.dat'), 'w')
    rcard.write(reweightCommand)
    rcard.close()


###################
# Event generation
###################
generate(runArgs=runArgs, process_dir=process_dir)

# multi-core capability
check_reset_proc_number(opts)

# put output into the appropriate place for the transform
arrange_output(process_dir=process_dir, runArgs=runArgs, lhe_version=3, saveProcDir=False)

# option: disable TestHepMC
# if hasattr(testSeq, "TestHepMC"):
#     testSeq.remove(TestHepMC())


##################
# Parton shower
##################
# showering with Pythia 8
evgenConfig.description = "Dark Higgs (bb) Dark Matter from 2MDM UFO"
evgenConfig.keywords = ["exotic", "BSM"]
evgenConfig.process = "generate p p > zp > n1 n1 hs, (hs > b b~)"

include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
include("Pythia8_i/Pythia8_MadGraph.py")

# Pythia settings: make the dark matter invisible
# syntax: particle data = name antiname spin=2s+1 3xcharge colour mass width (left out, so set to 0: mMin mMax tau0)
genSeq.Pythia8.Commands += ["SLHA:allowUserOverride = on",
                            "1000022:all = chi chi 2 0 0 %d 0.0 0.0 0.0 0.0" %(mdm),
                            "1000022:isVisible = false"]

# CKKW-L jet matching
PYTHIA8_nJetMax=1
PYTHIA8_Dparameter=float(settings['dparameter'])
PYTHIA8_Process="guess"
PYTHIA8_TMS=float(settings['ktdurham'])
PYTHIA8_nQuarksMerge=4
include("Pythia8_i/Pythia8_CKKWL_kTMerge.py")
genSeq.Pythia8.Commands+=["Merging:mayRemoveDecayProducts=on"]
# modification of merging to allow pythia to guess the hard process with "guess" syntax
if "UserHooks" in genSeq.Pythia8.__slots__.keys():
    genSeq.Pythia8.UserHooks += ['JetMergingaMCatNLO']
else:
    genSeq.Pythia8.UserHook = 'JetMergingaMCatNLO'
genSeq.Pythia8.CKKWLAcceptance = False

