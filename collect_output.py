from argparse import ArgumentParser
from os import walk
from os.path import join
import json
import numpy as np
import pandas as pd
from copy import deepcopy

def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input_dir', default='data')
    parser.add_argument('-o', '--output', default='summary.csv')
    parser.add_argument('--grid', default='grid.csv')
    return parser


def getGridAsDataFrame(grid_file):
    df = pd.read_csv(
        grid_file,
        names=['dsid', 'mzp', 'mdh', 'mdm', 'gq', 'gx'],
        dtype={'dsid': np.int32, 'mzp': np.int32, 'mdh': np.int32, 'mdm': np.int32, 'gq': np.float64, 'gx': np.float64}
    )
    return df


def getAcceptances(file_data):
    acc = {}
    for row in file_data:
        # ignore all lines which do not contain acceptances
        if not row.strip() or not row.startswith('MET'): continue
        try:
            name, events, acceptance, err = row.strip().split(',')
            acc[name] = float(acceptance)
        except ValueError:
            pass
    # calculate total acceptance
    acc['total'] = np.sum(np.array([float(i) for i in acc.values()]))
    return acc


def getData(input_dir, df_grid):
    data = {}
    for root, dirs, files in walk(input_dir):
        # d: e.g. 100000 (DSID)
        for d in sorted(dirs):
            tempdict = {}
            try:
                # look up DSID in grid file to get mzp, mdh, mdm, gq, gx 
                # signal parameters
                row = df_grid[df_grid['dsid'] == int(d)]
                tempdict['mzp'] = int(row['mzp'])
                tempdict['mdh'] = int(row['mdh'])
                tempdict['mdm'] = int(row['mdm'])
                tempdict['gq'] = float(row['gq'])
                tempdict['gx'] = float(row['gx'])

                # get cross-section from json dictionary
                with open(join(input_dir, d, 'xsec.json')) as f:
                    xsec_fb = json.loads(f.read())['xsec_fb']
                # get acceptances estimated with simple analysis
                with open(join(input_dir, d, 'acceptances.txt')) as f:
                    acc = getAcceptances(f)
                    for entry, value in acc.items():
                        tempdict['acc_' + entry] = value

                # get sensitivity estimate based on model-independent limits
                with open(join(input_dir, d, 'sensitivity.txt')) as f:
                    tempdict['sensitivity'] = f.read().strip()

                # add to final dict
                data[d] = deepcopy(tempdict)

            # skip file in case something was missing
            except (FileNotFoundError, KeyError):
                pass

    df = pd.DataFrame.from_dict(data, orient='index')
    df.index.name = 'dsid'
    return df


def main():
    args = getArgumentParser().parse_args()
    df_grid = getGridAsDataFrame(args.grid)
    data = getData(args.input_dir, df_grid)
    data.to_csv(args.output)
    print(data)
    print('output: ' + args.output)

if __name__ == "__main__":
    main()
