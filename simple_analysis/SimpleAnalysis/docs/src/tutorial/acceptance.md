---
title: Acceptance and efficiency plots
---

## Acceptance plots

Acceptance plots are typically part of the material your analysis will need to provide in their HEPdata. Since SimpleAnalysis outputs acceptance values for each region you defined, producing acceptance plots is rather easy. In order to keep things a little bit more interesting, we will not be using `ROOT` in the following plotting example, but rather `matplotlib`. You can, however, of course use any other plotting tool of your liking.

### Preliminaries

Make sure you are starting in a new shell where you haven't set up an analysis release. We won't need that here anyway and setting one up tends to unnecessarily complicate things.

First, start by creating the previous `TUTORIAL_DIR`environment variable again. In your base directory, do:
```sh
export TUTORIAL_DIR=$(pwd)
```

We will need to install a two python packages, `matplotlib`, `pandas` and `atlasify`. While `matplotlib` is probably already installed, `pandas` and `atlasify` typically are not if you haven't used them previously:
```sh
python3 -m pip install --user matplotlib pandas atlasify
```

!!! warning "Python 3"
    Notice how we are using `python3` now, instead of `python 2.7`. We will stick with that choice in this example.

!!! info "Virtual environment"
    You can always install these packages in a [virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) if you want to. Just do `python3 -m venv env` to create a virtual environment called `env` and activate it using `source env/bin/activate`.

Let's go ahead and create a new directory and an empty plotting script:
```sh
mkdir $TUTORIAL_DIR/acceptance
cd $TUTORIAL_DIR/acceptance
touch plot.py
```

### Imports

Using your favorite editor, open the `plot.py` file and start editing. We will first do some necessary imports. We'll need `glob` and `os` for conveniently loading the SimpleAnalysis output into our script. We'll also use `pandas` dataframes to store our acceptances in memory. Then, we will of course need `matplotlib` for plotting and `numpy` for tweaking some plotting stuff.
```python
#!/usr/bin/python3

import glob
import os

import pandas as pd

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
```

### Loading the SimpleAnalysis output

Next, let's define the basepath at which your copied inputs live as well as the signal region you want to plot.
```python
basepath = os.path.join(os.environ['TUTORIAL_DIR'],'inputs/acceptances/')
sr = 'EwkOneLeptonTwoBjets2018__SR_h_High_bin3'
```

The next step is to actually read in all the inputs. The following block of code will just loop over all files and read the text files into a `pandas` dataframe.
```python
# gather all the results from SimpleAnalysis
print("Getting inputs...")

acceptances = pd.DataFrame(columns=['m_x', 'm_y', 'SR', 'events', 'acceptance', 'err'])
for d in glob.glob(os.path.join(basepath,'*.txt')):
    with open(d, 'r') as file:
        file.readline() # get rid of header
        for l in file:
            replacements = [
              ('MadGraphPythia8EvtGen_A14N23LO_C1N2_Wh_hbb_',''),
              ('_lep_DAOD_TRUTH3.txt',''),
              ('p0','.0'),
              ('p5','.5'),
            ]
            masses = os.path.basename(d).split('.')[2]
            for (check, rep) in replacements:
              masses = masses.replace(check, rep)
            masses = masses.split('_')
            columns = l.split(',')
            acceptances = acceptances.append({
                'm_x' : float(masses[0]), # x-axis mass
                'm_y' : float(masses[1]), # y-axis mass
                'SR' : columns[0],
                'events' : float(columns[1]),
                'acceptance' : float(columns[2]),
                'err' : float(columns[3])}, ignore_index=True)
```

!!! warning "Filter efficiencies"
    If your input samples have been produced using a filter (e.g. filtering only events with one lepton in the final state), then you might want to consider to already include this filter efficiency in the acceptance. That way, you will be able to give the analysis acceptance considering the unfiltered process. In any case, you should clearly state whether or not any filter efficiency has already been considered.

### Plotting

Now, we start plotting. Most of the following code is standard `matplotlib` plotting. First, get the x-axis masses, y-axis masses as well as the acceptance values for the signal region we want to plot:
```python
print("Plotting...")

x_list = acceptances.loc[acceptances['SR'] == sr, 'm_x']
y_list = acceptances.loc[acceptances['SR'] == sr, 'm_y']
z_list = acceptances.loc[acceptances['SR'] == sr, 'acceptance']*100 # in percent
```

Now, we set the axis ranges and actually plot 2D contours using 100 contour levels with the `coolwarm` colourmap:
```python
plt.axis([150, 1100, 0, 500]) # axis ranges
cnt = plt.tricontourf(x_list,y_list,z_list,levels=np.linspace(0.0,6.0,100), cmap='coolwarm')
```

The 100 contour levels are used to get smooth transitions. This can sometimes cause weird optical effects in non-rasterized figures (e.g. in PDF format) because of the edgecolor of these contours. So we set all of these edge colors to `"face"` to prevent that from happening.
```python
# fix for the white lines between contour levels
# https://stackoverflow.com/a/32911283
for c in cnt.collections:
    c.set_edgecolor("face")
```

Next, we create a good looking colorbar on the right side of the plot. You could just call `plt.colorbar(cnt)`, and let matplotlib figure out the colorbar ticks and format, but if you want to control it yourself, this is how you would do it.
```python
def myfmt(x, pos):
    return '{0:.1f}'.format(x)

cbar = plt.colorbar(cnt, format=ticker.FuncFormatter(myfmt),ticks=[0.0,1.0,2.0,3.0,4.0,5.0,6.0])
cbar.set_label(r"Acceptance [%]")
```
And then we also want to see a scatter plot of the signal grid as well as the actual acceptance values:
```python
plt.scatter(x_list,y_list,s=5,facecolors='none', edgecolors='k')
for i,_ in enumerate(z_list):
    plt.annotate('{0:.2f}'.format(z_list.iloc[i]), (x_list.iloc[i],y_list.iloc[i]), size=6)
```
Last, but not least, we set the axis labels and add a text label to the plot. If you wish to switch to ATLAS style, you can do so using the `atlasify` package. You might need to install that with `pip3 install --user atlasify`.
```python
plt.xlabel(r"$m(\tilde{\chi}^{\pm}_{1}$/$\tilde{\chi}^{0}_{2})$ [GeV]")
plt.ylabel(r"$m(\tilde{\chi}^{0}_{1})$ [GeV]")
plt.text(0.04,0.86,r"Region: SRHM high $m_{CT}$",transform = plt.gca().transAxes, alpha=1.0, fontsize=10)

# do pip3 install --user atlasify before this
from atlasify import atlasify
atlasify("Internal", enlarge=1.0)
```

Before saving the plot to a PDF file, we tighten the layout as well.
```python
plt.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)
plt.savefig('plot.pdf')

print("Done...")
```
### Running and output

We live post April 2020, so the sun definitely has set on `python 2.x`. Python is dead, long live python, so let's use `python 3.x` to run this
```sh
python3 plot.py # is python2 really dead if we still have to write this?
```

This may give you a couple of deprecation warnings if you use `atlasify`, but that is nothing to worry about for now. Your output should look pretty much like this, a typical acceptance plot:

![Acceptance plot](images/acceptance_plot.png)

## Efficiency plots

!!! abstract "TL;DR"
    Creating efficiency plots is very similar to creating acceptance plots with the difference that you need some additional reco-level inputs from your analysis.

Efficiency plots encapsulate detector effects and show how efficient your object reconstruction is. Multiplying the efficiency with the acceptance (and any generator filter efficiency) gives you the full analysis selection efficiency using reconstructed objects, including any detector effects.

Generating efficiency plots using the SimpleAnalysis acceptance output is also quite easy. The efficiency in any given region is defined as the ratio between number of events using reconstructed objects and number of events using truth objects in that region:

$$
\epsilon = N_\mathrm{reco}/N_\mathrm{truth}
$$

The number of events at full-reconstruction level is taken from your usual analysis MC inputs (e.g. your analysis ntuples). The number of events at truth level can be easily computed using:

$$
N_\mathrm{truth} = a\cdot\sigma\cdot\epsilon_\mathrm{filter}\cdot\mathcal{L}_\mathrm{int}
$$

where $a$ is the acceptance from SimpleAnalysis, $\sigma$ is the cross-section of the process, $\epsilon_\mathrm{filter}$ the MC filter efficiency and $\mathcal{L}_\mathrm{int}$ the integrated luminosity.
