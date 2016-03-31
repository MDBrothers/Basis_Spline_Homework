#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from tabulate import tabulate
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
        if(x>=knot[i] and x<knot[i+1]):
            return 1.0
        elif( x==1.0 and knot[i+1] == 1.0): #exception to allow evaluating at 1.0
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

# Tell numpy how to apply the function to an array
vectorN = np.vectorize(N,excluded={0,1,3})

# For problem 1.1
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

#For problem 1.2
domain_u = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
knotvect_orig = knotvect
knotvect_a = np.array([0,0,0,.5,.5,.8,1,1,1])
knotvect_b = np.array([0,0,0,0,.5,.5,.8,.8,1,1,1,1])
knotvect_c = np.array([0,0,0,.3,.5,.8,1,1,1])

control_points_orig = np.array([[50, 50], [150,150], [250,50], [350,250], [450,250]])
control_points_a = np.array([[50, 50], [150,150], [435.0/2.0, 175.0/2.0], [250,50], [350,250], [450,250]])
control_points_b = np.array([[50,50], [117,117], [171,129], [238,62], [270,90], [337,223], [383,250], [450,250]])
control_points_c = np.array([[50, 50], [110,110], [375.0/2.0, 225.0/2.0], [250,50], [350,250], [450,250]])

N1_orig = vectorN(0,2,domain_u,knotvect)
N2_orig = vectorN(1,2,domain_u,knotvect)
N3_orig = vectorN(2,2,domain_u,knotvect)
N4_orig = vectorN(3,2,domain_u,knotvect)
N5_orig = vectorN(4,2,domain_u,knotvect)

N1_a = vectorN(0,2,domain_u,knotvect_a)
N2_a = vectorN(1,2,domain_u,knotvect_a)
N3_a = vectorN(2,2,domain_u,knotvect_a)
N4_a = vectorN(3,2,domain_u,knotvect_a)
N5_a = vectorN(4,2,domain_u,knotvect_a)
N6_a = vectorN(5,2,domain_u,knotvect_a)

N1_b = vectorN(0,3,domain_u,knotvect_b)
N2_b = vectorN(1,3,domain_u,knotvect_b)
N3_b = vectorN(2,3,domain_u,knotvect_b)
N4_b = vectorN(3,3,domain_u,knotvect_b)
N5_b = vectorN(4,3,domain_u,knotvect_b)
N6_b = vectorN(5,3,domain_u,knotvect_b)
N7_b = vectorN(6,3,domain_u,knotvect_b)
N8_b = vectorN(7,3,domain_u,knotvect_b)

N1_c = vectorN(0,2,domain_u,knotvect_c)
N2_c = vectorN(1,2,domain_u,knotvect_c)
N3_c = vectorN(2,2,domain_u,knotvect_c)
N4_c = vectorN(3,2,domain_u,knotvect_c)
N5_c = vectorN(4,2,domain_u,knotvect_c)
N6_c = vectorN(5,2,domain_u,knotvect_c)

N_orig = np.array([N1_orig, N2_orig, N3_orig, N4_orig, N5_orig])
N_a = np.array([N1_a, N2_a, N3_a, N4_a, N5_a, N6_a])
N_b = np.array([N1_b, N2_b, N3_b, N4_b, N5_b, N6_b, N7_b, N8_b])
N_c = np.array([N1_c, N2_c, N3_c, N4_c, N5_c, N6_c])
curve_orig = []
curve_a = []
curve_b = []
curve_c = []
for point in enumerate(domain_u):
    curve_orig.append( np.sum(control_points_orig.transpose()*N_orig[:,point[0]],axis=1) )
    curve_a.append( np.sum(control_points_a.transpose()*N_a[:,point[0]],axis=1) )
    curve_b.append( np.sum(control_points_b.transpose()*N_b[:,point[0]],axis=1) )
    curve_c.append( np.sum(control_points_c.transpose()*N_c[:,point[0]],axis=1) )
headers = ["u=0.0", "u=0.25", "u=0.5", "u=0.75", "u=1.0"]
print tabulate(N_orig, headers, tablefmt="latex")
print tabulate(N_a, headers, tablefmt="latex")
print tabulate(N_b, headers, tablefmt="latex")
print tabulate(N_c, headers, tablefmt="latex")
print tabulate(np.array(curve_orig).transpose(), headers, tablefmt="latex")
print tabulate(np.array(curve_a).transpose(), headers, tablefmt="latex")
print tabulate(np.array(curve_b).transpose(), headers, tablefmt="latex")
print tabulate(np.array(curve_c).transpose(), headers, tablefmt="latex")



# For problem 2.2
knotvec_U = np.array([0,0,0,0,2.4,2.4,2.4,2.4])
knotvec_V = np.array([0,0,0,1.2,1.7,1.7,1.7])
x_space = np.linspace(0.0,2.4,40)
y_space = np.linspace(0.0,1.7,40)
X, Y = np.meshgrid(x_space, y_space)
Z = []

from mpl_toolkits.mplot3d import Axes3D
newfig = plt.figure()
ax3d = newfig.add_subplot(111, projection='3d')

for i in range(4):
    for j in range(4):
        Z.append(  (np.array([ vectorN(i,3,x, knotvec_U)*vectorN(j,2,y,knotvec_V) for x,y in zip(np.ravel(X), np.ravel(Y))])).reshape(X.shape))
        ax3d.set_xlabel('X')
        ax3d.set_ylabel('Y')
        ax3d.set_title('Tensor product basis (i,j)=('+ str(i) +', '+ str(j) +')')
        ax3d.scatter(X,Y,Z[-1])
        savefig('problem_2_2_'+str(i)+'_'+str(j))
        ax3d.clear()
