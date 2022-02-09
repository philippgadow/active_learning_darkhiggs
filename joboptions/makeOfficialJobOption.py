import logging
from shutil import copy
from argparse import ArgumentParser
from os import makedirs, symlink, chdir
from os.path import join, basename


def getArguments():
    parser = ArgumentParser()
    parser.add_argument("--mzp", help="mass of Z' resonance", required=True)
    parser.add_argument("--mdh", help="mass of dark Higgs boson", required=True)
    parser.add_argument("--mdm", help="mass of dark matter particle", required=True)
    parser.add_argument("--gq", help="coupling of resonance to quarks", required=True)
    parser.add_argument("--gx", help="dark sector coupling", required=True)
    parser.add_argument("--dsid", help="DSID of job", required=True)
    parser.add_argument("--outputDir", help="Path to output directory", default='100xxx')
    parser.add_argument(
        "--template",
        default="mc.MGPy8EG_DMSbb_template.py",
        help="Path to job option template",
    )
    return parser


def makeJobOption(mzp, mdh, mdm, gq, gx, dsid, template, outputDir="", copyLog=False):
    logging.info("Creating job option with the following parameters:")
    logging.info("- mzp {mzp}".format(mzp=mzp))
    logging.info("- mdh {mdh}".format(mdh=mdh))
    logging.info("- mdm {mdm}".format(mdm=mdm))
    logging.info("- gq {gq}".format(gq=gq))
    logging.info("- gx {gx}".format(gx=gx))

    # make output directory
    makedirs(join(outputDir, dsid))

    # create JobOption from template
    with open(template) as fin:
        outFile = fin.read().replace("#GQ", "{0:.2f}".format(float(gq))).replace("#GX", "{0:.2f}".format(float(gx)))
        outFileName = "mc.MGPy8EG_DMSbb_zp{mzp}_dm{mdm}_dh{mdh}_gq{gq}_gx{gx}.py".format(
            mzp=mzp, mdh=mdh, mdm=mdm,
            gq="{0:.2f}".format(float(gq)).replace('.', 'p'),
            gx="{0:.2f}".format(float(gx)).replace('.', 'p')
        )
    with open(join(outputDir, dsid, outFileName), "w") as fout:
        fout.write(outFile)

    # link top level job option
    chdir(join(outputDir, dsid))
    symlink('../../508xxx/508820/MadGraphControl_MadGraphPythia8_N30NLO_A14N23LO_monoSbb_CKKWL.py', 'MadGraphControl_MadGraphPythia8_N30NLO_A14N23LO_monoSbb_CKKWL.py')
    chdir('../..')

    # copy log
    if copyLog: copy('../data/{dsid}/log.generate'.format(dsid=dsid), join(outputDir, dsid))

def main():
    args = getArguments().parse_args()
    makeJobOption(
        args.mzp,
        args.mdh,
        args.mdm,
        args.gq,
        args.gx,
        args.dsid,
        args.template,
        args.outputDir
    )


if __name__ == "__main__":
    main()
