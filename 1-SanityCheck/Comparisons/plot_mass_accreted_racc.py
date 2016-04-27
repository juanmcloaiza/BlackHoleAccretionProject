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
datafile = "./macc_racc_res.txt"
X = np.array(load_data(datafile))
R_acc = X[:,0]
Res = ['dummy','100k','250k','500k','750k','1M']
for i in range(1,6):
	label_this = Res[i]
	imgplot = pl.loglog(R_acc, X[:,i],'o', markersize = 0.7*FontSize, linewidth = 0.7*FontSize, label = label_this)

#plot fit
x = np.linspace(0.8e-3,6e-3,1025)
y = 2*x
imgplot = pl.loglog(x,y,'k--',label = '$ N_{\\rm acc}/N_{\\rm part} \\sim r_{\\rm acc} $')
for i in range(1,10):
	imgplot = pl.loglog(x,x/i**2,'k--')

pl.legend(loc=2)
pl.xlabel("$r_{\\rm acc}$",fontsize=1.5*FontSize)
pl.ylabel("$N_{\\rm acc}/N_{\\rm part}$",fontsize=1.5*FontSize)
#pl.xlim([2e-3,2e-2])
#pl.ylim([1e-3,1e0])
#pl.title("Circularization radius is where the peak is..." )
figfile = "macc_racc.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
