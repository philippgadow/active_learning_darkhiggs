from argparse import ArgumentParser
from os.path import join
import json


def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('-a', '--acceptances')
    parser.add_argument('-x', '--cross_section')
    parser.add_argument('-l', '--limits', default=join('data', 'limits_EXOT-2018-46.csv'))
    parser.add_argument('-o', '--output_file')
    return parser


def getAcceptances(fName):
    acceptances = {}
    with open(fName) as f:
        # skip first two header lines
        for l in f.readlines()[2:]:
            # line content: region, events, acceptance, uncertainty
            l = l.strip().split(',')
            acceptances[l[0]] = float(l[2])
    return acceptances


def getCrossSection(fName):
    with open(fName) as f:
        return float(json.load(f)['xsec_fb'])


def getLimits(fName):
    limits = {}
    with open(fName) as f:
        # skip header line
        for l in f.readlines()[1:]:
            # line content: region, observed limit, expected limit
            l = l.strip().split(',')
            limits[l[0]] = float(l[1])
    return limits


def main():
    args = getArgumentParser().parse_args()

    acceptances = getAcceptances(args.acceptances)
    cross_section = getCrossSection(args.cross_section)
    limits = getLimits(args.limits)

    sensitivity = 0.
    sensitivities = {}
    for region in acceptances.keys():
        sensitivity_region = (cross_section * acceptances[region]) / limits[region]
        sensitivities[region] = sensitivity_region
        sensitivity += sensitivity_region

    with open(args.output_file, 'w') as f:
        f.write(str(sensitivity) + '\n')

    print(sensitivity)
    if (sensitivity > 1):
        print('Signal point is excluded.')

if __name__ == '__main__':
    main()
