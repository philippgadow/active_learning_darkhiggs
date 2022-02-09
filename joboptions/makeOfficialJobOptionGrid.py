import argparse
from makeOfficialJobOption import makeJobOption

def getArgumentParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('grid', help='Path to signal grid to be generated')
    return parser


def main():
    args = getArgumentParser().parse_args()
    # parameters (eventually to be outsourced to a config file)
    # dsid, mzp,mdh,mdm,gq,gx
    job_parameters = []
    with open(args.grid) as f:
        for l in f: job_parameters.append(l.strip())

    outputDir = '100xxx'
    template = 'mc.MGPy8EG_DMSbb_template.py'
    for data in job_parameters:
        if not data.split(): continue
        dsid, mzp, mdh, mdm, gq, gx = data.split(',')
        mzp = round(float(mzp))
        mdh = round(float(mdh))
        mdm = round(float(mdm))

        makeJobOption(
            mzp,
            mdh,
            mdm,
            gq,
            gx,
            dsid,
            template,
            outputDir,
            copyLog=True
        )


if __name__ == '__main__':
    main()
