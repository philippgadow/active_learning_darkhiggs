from argparse import ArgumentParser
from itertools import product
import numpy as np


def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('--mzp_min', default=500)
    parser.add_argument('--mzp_max', default=3500)
    parser.add_argument('--mzp_step', default=500)

    parser.add_argument('--mdh_min', default=50)
    parser.add_argument('--mdh_max', default=170)
    parser.add_argument('--mdh_step', default=20)

    parser.add_argument('--mdm_min', default=50)
    parser.add_argument('--mdm_max', default=500)
    parser.add_argument('--mdm_step', default=50)

    parser.add_argument('--gx', nargs='+', default=[0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0])

    parser.add_argument('--dsid_start', default=100000)

    return parser

def main():
    args = getArgumentParser().parse_args()
    dsid = args.dsid_start

    output = ''

    for mzp, mdh, mdm, gx in product(
        np.arange(args.mzp_min, args.mzp_max + args.mzp_step, args.mzp_step),
        np.arange(args.mdh_min, args.mdh_max + args.mdh_step, args.mdh_step),
        np.arange(args.mdm_min, args.mdm_max + args.mdm_step, args.mdm_step),
        np.array(args.gx)
    ):
        line = "{dsid},{mzp},{mdh},{mdm},0.25,{gx:.2f}".format(
            dsid=dsid, mzp=mzp, mdh=mdh, mdm=mdm, gx=gx
        )
        dsid += 1
        output += line + '\n'

    with open('grid.csv', 'w') as f:
        f.write(output)


if __name__ == "__main__":
    main()
