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

#plot shaded region inside r_acc
#pl.axvspan(1e-3, 5e-3, facecolor='k', alpha=0.5)

snap = 5

#Plot the simulations without accretion radius:
for res, color in zip(['001M','002M','004M','008M','010M'],['b-','g-','r-','c-','m-']):
	datafile = "./res_"+res+"_t"+str(snap)+"00.txt"
	X = np.array(load_data(datafile))
	M_cumulative = np.array(X[:,1])
	for i in range(1,len(M_cumulative)):
		M_cumulative[i] += M_cumulative[i-1]
	imgplot = pl.loglog(X[:,0], M_cumulative,color, linewidth = 0.2*FontSize)


pl.legend(loc=2)
pl.xlabel("$r$",fontsize=1.5*FontSize)
pl.ylabel("$N_{\\rm frac}$",fontsize=1.5*FontSize)
#pl.xlim([2e-3,2e-2])
#pl.ylim([1e-3,1e0])
#pl.title("Circularization radius is where the peak is..." )
figfile = "M_inside_r_all_res.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()
