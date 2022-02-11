# coupling settings
gq = 0.25
gx = 2.09
th = 0.01

# mass settings are determined by file name structure:
# mzp <- int(filename.split('_zp')[1].split('_')[0])
# mdm <- int(filename.split('_dm')[1].split('_')[0])
# mhs <- int(filename.split('_dh')[1].split('_')[0])

# default number of events
evgenConfig.nEventsPerJob = 10000
include("MadGraphControl_MadGraphPythia8_N30NLO_A14N23LO_monoSbb_CKKWL.py")

