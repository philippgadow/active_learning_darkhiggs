from argparse import ArgumentParser
from os import walk
from os.path import join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input_dir', default='data')
    return parser


def getSensitivities(input_dir):
    sensitivities = {}

    for root, dirs, files in walk(input_dir):
        for d in sorted(dirs):
            try:
                with open(join(input_dir, d, 'sensitivity.txt')) as f:
                    sensitivities[d] = f.read().strip()
            except FileNotFoundError:
                pass

    return sensitivities


def prepareDataFrame(data, sensitivities):
    # prepare dataframe with sensitivities
    df = pd.DataFrame(columns=['mzp', 'mdm', 'ms', 'gx', 'gq', 'sensitivity'])

    for dsid, sensitivity in sensitivities.items():
        mzp = data[dsid][0]
        ms = data[dsid][1]
        mdm = data[dsid][2]
        gq = data[dsid][3]
        gx = data[dsid][4]
        df = df.append({'mzp': mzp, 'ms': ms, 'mdm': mdm,
                        'gq': gq, 'gx': gx, 'sensitivity': sensitivity},
                       ignore_index=True)

    # prepare dataframe
    df = df.replace(r'^\s*$', np.nan, regex=True)
    df = df.astype(float).dropna()
    df = df[df['mdm'] == 200.]
    df = df.pivot('ms', 'mzp', 'sensitivity')

    return df


def plotHeatmap(df):
    # plot limits as heatmap
    f, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(df, ax=ax, annot=True, cmap="coolwarm_r", center=1., fmt='.1f', cbar_kws={'label': 'sensitivity'}, vmin=0.00, vmax=1.2)

    # style plot
    plt.xlabel("Mediator ($Z'$) mass [GeV]", position=(1., 0.), va='bottom', ha='right')
    plt.ylabel('Dark Higgs ($s$) mass [GeV]', position=(0., 1.), va='top', ha='right')
    ax.xaxis.set_label_coords(1., -0.10)
    ax.yaxis.set_label_coords(-0.18, 1.)
    plt.gca().invert_yaxis()
    outName = 'sensitivity_heatmap.png'
    plt.savefig(outName)


def formatForPlot(data):
    X=data.index.values
    Y=data.columns.values
    Z=data.values
    Xi,Yi = np.meshgrid(Y, X)
    return Xi, Yi, Z


def plotContour(df):
    # plot limits as contour
    f, ax = plt.subplots(figsize=(10, 8))
    relicdensity_plt = plt.contour(*formatForPlot(df), levels=[1.0], colors='black', linewidths=[2.], linestyles='dotted');

     # style plot
    plt.xlabel("Mediator ($Z'$) mass [GeV]", position=(1., 0.), va='bottom', ha='right')
    plt.ylabel('Dark Higgs ($s$) mass [GeV]', position=(0., 1.), va='top', ha='right')
    ax.xaxis.set_label_coords(1., -0.10)
    ax.yaxis.set_label_coords(-0.18, 1.)
    outName = 'sensitivity_contour.png'
    plt.savefig(outName)


def main():
    args = getArgumentParser().parse_args()

    # style
    plt.style.use('https://raw.githubusercontent.com/beojan/atlas-mpl/master/atlas_mpl_style/stylesheets/atlas.mplstyle')

    # get sensitivities from output folders
    sensitivities = getSensitivities(args.input_dir)

    # parameters (eventually to be outsourced to a config file)
    # dsid, mzp,mdh,mdm,gq,gx
    job_parameters = [
        "100000,500,50,200,0.25,1.0",
        "100001,500,70,200,0.25,1.0",
        "100002,500,90,200,0.25,1.0",
        "100003,500,110,200,0.25,1.0",
        "100004,500,130,200,0.25,1.0",
        "100005,500,150,200,0.25,1.0",

        "100006,1000,50,200,0.25,1.0",
        "100007,1000,70,200,0.25,1.0",
        "100008,1000,90,200,0.25,1.0",
        "100009,1000,110,200,0.25,1.0",
        "100010,1000,130,200,0.25,1.0",
        "100011,1000,150,200,0.25,1.0",

        "100012,1500,50,200,0.25,1.0",
        "100013,1500,70,200,0.25,1.0",
        "100014,1500,90,200,0.25,1.0",
        "100015,1500,110,200,0.25,1.0",
        "100016,1500,130,200,0.25,1.0",
        "100017,1500,150,200,0.25,1.0",

        "100018,2000,50,200,0.25,1.0",
        "100019,2000,70,200,0.25,1.0",
        "100020,2000,90,200,0.25,1.0",
        "100021,2000,110,200,0.25,1.0",
        "100022,2000,130,200,0.25,1.0",
        "100023,2000,150,200,0.25,1.0",

        "100024,2500,50,200,0.25,1.0",
        "100025,2500,70,200,0.25,1.0",
        "100026,2500,90,200,0.25,1.0",
        "100027,2500,110,200,0.25,1.0",
        "100028,2500,130,200,0.25,1.0",
        "100029,2500,150,200,0.25,1.0",

        "100030,3000,50,200,0.25,1.0",
        "100031,3000,70,200,0.25,1.0",
        "100032,3000,90,200,0.25,1.0",
        "100033,3000,110,200,0.25,1.0",
        "100034,3000,130,200,0.25,1.0",
        "100035,3000,150,200,0.25,1.0",
        
        "100036,3500,50,200,0.25,1.0",
        "100037,3500,70,200,0.25,1.0",
        "100038,3500,90,200,0.25,1.0",
        "100039,3500,110,200,0.25,1.0",
        "100040,3500,130,200,0.25,1.0",
        "100041,3500,150,200,0.25,1.0"
    ]

    data = {l.split(',')[0]: l.split(',')[1:] for l in job_parameters}
    df = prepareDataFrame(data, sensitivities)
    
    plotHeatmap(df)
    plotContour(df)
    

if __name__ == '__main__':
    main()
