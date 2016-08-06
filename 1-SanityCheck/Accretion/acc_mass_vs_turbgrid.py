#!/usr/bin/env python
import matplotlib.pyplot as pl
import numpy as np
import sys

def load_data(filename):
    all_lines = []
    with open(filename) as f:
        for line in f:
            if (line[0] == '#'):
                continue
            #for char in line:
            try:
                line_data = np.array(line.split(),dtype=float)
            except ValueError:
                print("Line:")
                print(line)
                print("will give error.")
                print("Check the input file:")
                print(filename)
                print(" for possible corruption.")
            all_lines.append(line_data)
    return all_lines

def prepare_figure():
    #Set relevant dimensions:
    FigSize = 10
    FontSize = 25
    BorderWidth = 3
    pl.rcParams.update({'font.size': FontSize})
    pl.rcParams.update({'axes.linewidth': BorderWidth})
    pl.rcParams.update({'lines.linewidth': 0.3*FontSize})
    pl.rcParams.update({'lines.markersize': 0.7*FontSize})

    #Open figure:
    pl.figure( figsize=( FigSize, FigSize ) )
    pl.tick_params(width=BorderWidth, length=FontSize, which='major')
    pl.tick_params(width=BorderWidth, length=0.3*FontSize, which='minor')

    pl.xlabel('$N_{\\rm grid}$')
    pl.ylabel('$M_{\\rm acc}/M_{\\rm shell}$')
    pl.xlim([64,1024])
    pl.ylim([1e-6,5e-2])


def get_acc_mass(datafile):
    #Read datafile
    DATA = np.array(load_data(datafile))
    Time = DATA[:,0] 
    Nfrac = (DATA[:,1]) / float(DATA[0,1])

    #The total accreted mass:
    Tot_acc_mass = Nfrac[0] - Nfrac[-1]
    return Tot_acc_mass

def get_fit(gridsize,M_measured):
    point_b = -1
    point_a = 0
    k = np.log(M_measured[point_b]/M_measured[point_a]) / np.log(gridsize[point_a]/gridsize[point_b])
    C0 = (M_measured[1] * gridsize[1] **(k))
    x = np.linspace(gridsize[0],gridsize[-1],1e5)

    y =  C0 * x**(-k)
    return x,y

def trace_the_line(X,Y,col):
    pl.loglog( X, Y, col+'o', label = "", basex=2 )
    pl.legend(loc=3)
    x,y = get_fit(X,Y)
    pl.loglog(x,y, col+'--',basex=2)
    #pl.loglog(x,y+Y[-1], col+'--')


prepare_figure()

X = []
Y = []
for gridsize in ([128,256,512]):
    datafile = 'TurbGrid_'+str(gridsize)+'.txt'
    Macc = get_acc_mass(datafile)
    X.append(gridsize)
    Y.append(Macc)
trace_the_line(X,Y,'b')


figfile = 'plot_acc_mass_vs_turbgrid.png'
pl.savefig(figfile)
print('Figure saved to '+figfile)
