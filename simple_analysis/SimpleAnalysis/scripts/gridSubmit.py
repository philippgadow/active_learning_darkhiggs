#!/usr/bin/env python

import commands
import os
import re
import sys

def main(args):
    if len(args)!=1: 
        print "Usage: gridSubmit.py <fileList>"
        sys.exit(1)
    listName=args[0]
    if not os.access(listName,os.R_OK):
        print "Cannot read input list:",listName
        sys.exit(1)

    if os.system("which prun >& /dev/null"):
        print "Panda not setup"
        sys.exit(1)

    if os.system("voms-proxy-info --valid 1:00 >& /dev/null"):
        print "Voms proxy not valid for at least one hour"
        print "call: voms-proxy-init --voms atlas"
        sys.exit(1)

    output = commands.getoutput('voms-proxy-info -all')
    for line in output.split('\n'):
        if line.startswith('attribute'):
            match = re.search('nickname =\s*([^\s]+)\s*\(atlas\)',line)
            if match != None:
                nickName = match.group(1)
                break
    # check        
    if nickName == '':
        print "Could not figure out grid nickname??"
        sys.exit(1)

    os.system("rc clean")
    print "Generating tar file"
    os.system("prun --exec 'simpleAnalysis -o results %IN' --outputs='*.txt,*.root' --useRootCore --outTarBall source.tar  --maxFileSize=3000000 --noSubmit --outDS user."+nickName+".test1234" )
    
    for line in open(listName).readlines():
        dsName=line.strip()
        print "input:",dsName
        dsParts=dsName.split(".")
        outName="user.%s.%s.%s.simple.v1" % (nickName,dsParts[1],dsParts[2]) # assumes standard derivation name
        print "output:",outName
        os.system("prun --exec 'simpleAnalysis -o results %IN' --outputs='*.txt,*.root' --useRootCore --inTarBall source.tar  --maxFileSize=3000000 --inDS "+dsName+" --outDS "+outName )
    os.system("rm source.tar")

if __name__ == "__main__":
   sys.exit(main(sys.argv[1:]))
