#!/usr/bin/env python
import matplotlib.pyplot as pl
import numpy as np

def load_data(filename):
	all_lines = []
	with open(filename) as f:
		for line in f:
			if (line[0] == '#'):
				continue
			#for char in line:
			line_data = np.array(line.split(),dtype=float)
			all_lines.append(line_data)
	return all_lines

#Set relevant dimensions:
FigSize = 10
FontSize = 20
BorderWidth = 3
pl.rcParams.update({'font.size': FontSize})
pl.rcParams.update({'axes.linewidth': BorderWidth})
#Open figure:
pl.figure( figsize=( FigSize, FigSize ) )
pl.tick_params(width=BorderWidth, length=FontSize, which='major')
pl.tick_params(width=BorderWidth, length=0.3*FontSize, which='minor')

#plot data for Racc = 1e-3
datafile = "../Racc_1e-3/M_inside_r/Npart.txt"
X = np.array(load_data(datafile))
Resolution = X[:,0]
M_accreted = X[:,4]
r_acc = 1e-3
label_this = '$r_{\\rm acc} = 10^{-3}$'
imgplot = pl.loglog(Resolution, M_accreted,'bo', markersize = 0.7*FontSize, linewidth = 0.7*FontSize, label = label_this)

#plot data for Racc = 5e-3
datafile = "../Racc_5e-3/M_inside_r/Npart.txt"
X = np.array(load_data(datafile))
Resolution = X[:,0]
M_accreted = X[:,4]
r_acc = 1e-3
label_this = '$r_{\\rm acc} = 5\\times10^{-3}$'
imgplot = pl.loglog(Resolution, M_accreted,'ro', markersize = 0.7*FontSize, linewidth = 0.7*FontSize, label = label_this)

#plot fit
#x = np.linspace(1e5,1e6,1025)
#y = 1.5e-3*(x/1e5)**-1.5
#imgplot = pl.loglog(x,y,'k--',label = '$y \\sim x^{-3/2} $')

pl.legend(loc=1)
pl.xlabel("$N_{\\rm part}$",fontsize=1.5*FontSize)
pl.ylabel("$N_{\\rm acc}/N_{\\rm part}$",fontsize=1.5*FontSize)
#pl.xlim([2e-3,2e-2])
#pl.ylim([1e-3,1e0])
#pl.title("Circularization radius is where the peak is..." )
figfile = "M_accreted.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
