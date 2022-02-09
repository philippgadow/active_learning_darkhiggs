import MadGraphControl.MadGraph_NNPDF30NLO_Base_Fragment
from MadGraphControl.MadGraphUtils import *

safefactor=2.5
gridpack_mode = False

evgenConfig.contact = ["Changqiao Li <changqiao.li@cern.ch>", "Qibin Liu <qibin.liu@cern.ch>", "Paul Philipp Gadow <pgadow@cern.ch>"]
evgenConfig.generators = ["MadGraph", "Pythia8", "EvtGen"]

# define process
from MadGraphControl.MadGraphUtilsHelpers import get_physics_short
phys_short = get_physics_short()

# Get other mass parameters from the phys short
mzp = int(phys_short.split('_zp')[1].split('_')[0])
mdm = int(phys_short.split('_dm')[1].split('_')[0])
mhs = int(phys_short.split('_dh')[1].split('_')[0])

process_string = """
generate p p > zp > n1 n1 hs QED<=2, (hs > b b~) @0
add process p p > zp > n1 n1 hs j QED<=2, (hs > b b~) @1
"""

# Form full process string and set up directory
process = """
import model DarkHiggs2MDM
{process}
output -f
""".format(process=process_string)
process_dir = new_process(process)

# determine ktdurham cut from dark Higgs mass
# (ktdurham cut sets scale at which event description is split between parton shower and matrix element)
ktdurham = 40

# fetch default LO run_card.dat and set parameters
nevents = runArgs.maxEvents*safefactor if runArgs.maxEvents>0 else safefactor*evgenConfig.nEventsPerJob
settings = {'lhe_version':'3.0',
            #'pdlabel'    : "'lhapdf'",
            #'lhaid'      : 247000,
            'cut_decays': 'F',
            'event_norm': 'sum',
            'drjj': "0.0",         # required for CKKW-L jet matching
            'ickkw': 0,            # required for CKKW-L jet matching
            'ktdurham': ktdurham,  # required for CKKW-L jet matching
            'dparameter': "0.4",    # required for CKKW-L jet matching
            'xqcut': "0.0",    # required for CKKW-L jet matching
            'nevents': nevents,
            #'bwcutoff': 25
          }
modify_run_card(process_dir=process_dir,runArgs=runArgs,settings=settings)

# write parameter card
params = {}
# Set masses
params['MASS'] = {'54':mhs,'55':mzp,'1000022':mdm}
# Set couplings
params['frblock'] = { '1':gq , '2':gx , '3':th }
# Set decay width
params['DECAY'] = { '54':"AUTO" , '55':"AUTO" }
# Get the card set up
modify_param_card(process_dir=process_dir,params=params)

# Perform the real event generation
generate(runArgs=runArgs,grid_pack=gridpack_mode,process_dir=process_dir)

# multi-core capability
check_reset_proc_number(opts)

# Put output into the appropriate place for the transform
arrange_output(process_dir=process_dir,runArgs=runArgs,lhe_version=3,saveProcDir=False)


# showering with Pythia 8
evgenConfig.description = "Dark Higgs (bb~) Dark Matter from 2MDM UFO"
evgenConfig.keywords = ["exotic","BSM"]
evgenConfig.process = "generate p p > zp > n1 n1 hs QED<=2 (+1 jet), (hs > b b~)"

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
