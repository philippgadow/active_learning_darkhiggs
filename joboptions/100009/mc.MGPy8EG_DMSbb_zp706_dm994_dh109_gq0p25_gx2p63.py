# coupling settings
gq = 0.25
gx = 2.63
th = 0.01

# mass settings are determined by file name structure:
# mzp <- int(filename.split('_zp')[1].split('_')[0])
# mdm <- int(filename.split('_dm')[1].split('_')[0])
# mhs <- int(filename.split('_dh')[1].split('_')[0])

# default number of events
evgenConfig.nEventsPerJob = 10000

# activate MadGraph reweight in g(x)
reweight = False

include("MadGraphControl_MadGraphPythia8_N31LO_A14N23LO_DMSbb_CKKWL.py")
