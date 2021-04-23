import uproot
import json
from argparse import ArgumentParser
from os.path import join


def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('inputFile')
    parser.add_argument('-e', '--events', type=float, help='Number of generated events')
    parser.add_argument('-x', '--cross_section', help='Path to json file containing cross-section (in fb)')
    parser.add_argument('-l', '--luminosity', type=float, default=139., help='Luminosity (in fb)')
    parser.add_argument('-t', '--template', default=join('data', 'patch_template.json'), help='Path to json file template')
    parser.add_argument('-o', '--output_file', default='patch.json', help='Output json file')
    parser.add_argument('--signalname', default='signal')
    parser.add_argument('--mzp', default='0')
    parser.add_argument('--mdh', default='0')
    parser.add_argument('--mdm', default='0')
    parser.add_argument('--gq', default='0.25')
    parser.add_argument('--gx', default='1.00')
    return parser


def getCrossSection(fName):
    with open(fName) as f:
        return float(json.load(f)['xsec_fb'])


def formatYieldArray(yield_array):
    result = ""
    for i in yield_array:
        result += "                            {i:f},\n".format(i=i)
    return '\n' + result[:-2] + '\n                        '


def dumpYields(yields, args):
    with open(args.template) as f:
        output = f.read()

    output = output.replace('#SIGNALNAME', args.signalname.replace('.', 'p'))
    output = output.replace('#MZP', args.mzp)
    output = output.replace('#MDH', args.mdh)
    output = output.replace('#MDM', args.mdm)
    output = output.replace('#GQ', args.gq)
    output = output.replace('#GX', args.gx)

    for k, yield_array in yields.items():
        output = output.replace('#REGION_' + k, formatYieldArray(yield_array))
    return output


def main():
    args = getArgumentParser().parse_args()
    cross_section = getCrossSection(args.cross_section)

    yields = {}
    with uproot.open(args.inputFile) as f:
        yields['MET150200_2b'] = f['MET150200_2b'].values()
        yields['MET200350_2b'] = f['MET200350_2b'].values()
        yields['MET350500_2b'] = f['MET350500_2b'].values()
        yields['MET500750_2b'] = f['MET500750_2b'].values()
        yields['MET750_2b'] = f['MET750_2b'].values()
        yields['MET150200_3b'] = f['MET150200_3b'].values()
        yields['MET200350_3b'] = f['MET200350_3b'].values()
        yields['MET350500_3b'] = f['MET350500_3b'].values()
        yields['MET500_3b'] = f['MET500_3b'].values()

    for k in yields.keys():
        yields[k] *= (args.luminosity * cross_section / args.events)

    output_string = dumpYields(yields, args)
    with open(args.output_file, 'w') as f:
        f.write(output_string)


if __name__ == "__main__":
    main()
