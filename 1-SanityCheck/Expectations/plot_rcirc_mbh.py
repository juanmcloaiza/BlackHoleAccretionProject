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

datafile = ".//fix_l_018.txt"
print("reading " + datafile)
X = np.array(load_data(datafile))
label_this="Analytical: ${l_0} = 0.018$"
imgplot = pl.loglog(X[:,0], X[:,1], label = label_this)

datafile = "../Fix_v/RcircCalc/Particles2Bins/rcirc_vs_mbh.txt"
print("reading " + datafile)
X = np.array(load_data(datafile))
label_this="Simulations: ${l_0} = 0.018$"
imgplot = pl.loglog(X[:,0], X[:,1], 'ro' ,label = label_this)

pl.legend(loc=1)
pl.xlabel("${\\rm M_{BH}}$")
pl.ylabel("$r_{\\rm circ}$")
#pl.xlim([0,1])
#pl.ylim([0,1])
pl.title("Circularization radius as function of black hole mass")
figfile = "rcirc_mbh.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()
