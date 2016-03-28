#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from scipy.optimize import bisect

mpl.use('pgf')

def figsize(scale):
    fig_width_pt = 550.0                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause scatters to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "text.fontsize": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # scatters will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

# I make my own newfig and savefig functions
def newfig(width):
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename):
    plt.savefig('{}.pdf'.format(filename))


fig, ax  = newfig(0.6)

# Cox-DeBoor formula for basis functions
def N(i,p,x,knot):
    if(p==0):
        if(x>=knot[i] and x<=knot[i+1]):
            return 1.0
        else:
            return 0.0
    else:
        a = x-knot[i]
        b = knot[i+p]-knot[i]
        if(b==0.0):
            a=0.0
            b=1.0
        c = knot[i+p+1] -x
        d = knot[i+p+1]-knot[i+1]
        if(d==0.0):
            c=0.0
            d=1.0
        return a /  b * N(i,p-1,x,knot) + c / d * N(i+1,p-1,x,knot)

vectorN = np.vectorize(N,excluded={0,1,3})

knotvect = np.array([0.0, 0.0, 0.0, 0.5, 0.8, 1.0, 1.0, 1.0])
domain = np.linspace(0.0,1.0,100)
N1 = vectorN(0,2,domain,knotvect)
N2 = vectorN(1,2,domain,knotvect)
N3 = vectorN(2,2,domain,knotvect)
N4 = vectorN(3,2,domain,knotvect)
N5 = vectorN(4,2,domain,knotvect)

ax.plot(domain, N1)
ax.plot(domain, N2)
ax.plot(domain, N3)
ax.plot(domain, N4)
ax.plot(domain, N5)
ax.set_xlabel('x')
ax.set_ylabel('N(x)')
ax.set_title('B-splines degree 2')
savefig('problem_1_1_1')
ax.clear()
