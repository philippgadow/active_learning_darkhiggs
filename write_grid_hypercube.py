from argparse import ArgumentParser
from itertools import product
import numpy as np
np.random.seed(123)
import pandas as pd

import skopt
from skopt.sampler import Lhs
from skopt.sampler import Grid
from scipy.spatial.distance import pdist



def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('--mzp_min', default=500)
    parser.add_argument('--mzp_max', default=5000)

    parser.add_argument('--mdh_min', default=50)
    parser.add_argument('--mdh_max', default=170)

    parser.add_argument('--mdm_min', default=100)
    parser.add_argument('--mdm_max', default=1210)

    parser.add_argument('--gx_min', default=0.1)
    parser.add_argument('--gx_max', default=3.5)

    parser.add_argument('--dsid_start', default=100000)
    parser.add_argument('--n_samples', default=25000)

    return parser


def main():
    args = getArgumentParser().parse_args()

    space = [
        skopt.space.Real(args.mzp_min, args.mzp_max, name='mzp', prior='uniform'),
        skopt.space.Real(args.mdh_min, args.mdh_max, name='mdh', prior='uniform'),
        skopt.space.Real(args.mdm_min, args.mdm_max, name='mdm', prior='uniform'),
        skopt.space.Real(args.gx_min, args.gx_max, name='g', prior='log-uniform')
    ]

    # sample data
    lhs = Lhs(lhs_type="centered", criterion=None)
    x = np.array(lhs.generate(space, args.n_samples, random_state=42))

    # set up dataframe for writing to file
    df = pd.DataFrame(data=x, columns=["mzp", "mdh", "mdm", "g"])
    df.index.rename('dsid', inplace=True)
    df.index = np.arange(args.dsid_start, args.dsid_start + len(df))

    # write data to file
    df.to_csv('./grid_hypercube.csv', header=False)


if __name__ == "__main__":
    main()
