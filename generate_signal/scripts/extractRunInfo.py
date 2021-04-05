from argparse import ArgumentParser
from os.path import join
import re
import json


def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument("-i", "--inputFile", default=join("workdir_evgen", "log.generate"))
    parser.add_argument("-o", "--outputFile", default='xsec.json')
    return parser


def main():
    args = getArgumentParser().parse_args()

    data = {}

    with open(args.inputFile) as f:
        for l in f.readlines():
            if not 'MetaData:' in l: continue
            
            if 'cross-section' in l:
                data['xsec_nb'] = float(l.split('=')[1].replace('\n', ''))
                data['xsec_fb'] = float(data['xsec_nb']) * 1e6

            if 'GenFiltEff' in l:
                data['filter_eff'] = l.split('=')[1].replace('\n', '')

    with open(args.outputFile, 'w') as fp:
        json.dump(data, fp, sort_keys=True)


if __name__ == "__main__":
    main()
