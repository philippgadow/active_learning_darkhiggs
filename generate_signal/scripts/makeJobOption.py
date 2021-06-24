import logging
from shutil import copy
from argparse import ArgumentParser
from os import makedirs
from os.path import join, basename


def getArguments():
    parser = ArgumentParser()
    parser.add_argument("--mzp", help="mass of Z' resonance", required=True)
    parser.add_argument("--mdh", help="mass of dark Higgs boson", required=True)
    parser.add_argument("--mdm", help="mass of dark matter particle", required=True)
    parser.add_argument("--gq", help="coupling of resonance to quarks", required=True)
    parser.add_argument("--gx", help="dark sector coupling", required=True)
    parser.add_argument("--dsid", help="DSID of job", required=True)
    parser.add_argument("--outputDir", help="Path to output directory", default='jobOptions')
    parser.add_argument(
        "--template",
        default=join("assets", "mc.MGPy8EG_DMSbb_template.py"),
        help="Path to job option template",
    )
    parser.add_argument(
        "--topJobOption",
        default=join(
            "assets", "MadGraphControl_MadGraphPythia8_N31LO_A14N23LO_DMSbb_CKKWL.py"
        ),
        help="Path to top level job option",
    )
    return parser


def makeJobOption(mzp, mdh, mdm, gq, gx, dsid, template, topJobOption, outputDir=""):
    logging.info("Creating job option with the following parameters:")
    logging.info("- mzp {mzp}".format(mzp=mzp))
    logging.info("- mdh {mdh}".format(mdh=mdh))
    logging.info("- mdm {mdm}".format(mdm=mdm))
    logging.info("- gq {gq}".format(gq=gq))
    logging.info("- gx {gx}".format(gx=gx))

    # make output directory
    try:
        makedirs(join(outputDir, dsid))
    except OSError:
        pass

    # create JobOption from template
    with open(template) as fin:
        outFile = fin.read().replace("#GQ", gq).replace("#GX", gx)
        outFileName = "mc.MGPy8EG_DMSbb_zp{mzp}_dm{mdm}_dh{mdh}_gq{gq}_gx{gx}.py".format(
            mzp=mzp, mdh=mdh, mdm=mdm,
            gq="{0:.2f}".format(float(gq)).replace('.', 'p'),
            gx="{0:.2f}".format(float(gx)).replace('.', 'p')
        )
    with open(join(outputDir, dsid, outFileName), "w") as fout:
        fout.write(outFile)

    # copy top level job option
    copy(topJobOption, join(outputDir, dsid, basename(topJobOption)))


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
        args.topJobOption,
        args.outputDir
    )


if __name__ == "__main__":
    main()
