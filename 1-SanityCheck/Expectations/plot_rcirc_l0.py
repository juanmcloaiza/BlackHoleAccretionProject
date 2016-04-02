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

mbhfiles = [0,-1,-2,-3,-4]
colors =   ["r","g","k","m","c"]
for mbhfile in mbhfiles:
	datafile = ".//fix_mbh_1e"+str(mbhfile)+".txt"
	print("reading " + datafile)
	X = np.array(load_data(datafile))
	label_this="${\\rm M_{BH}} = 10^{"+str(mbhfile)+"}$"
	imgplot = pl.loglog(X[:,0], X[:,1], label = label_this)
pl.legend(loc=2)
pl.xlabel("$l_{0}$")
pl.ylabel("$r_{\\rm circ}$")
pl.xlim([1e-4,0.018])
pl.ylim([1e-6,1e-1])
pl.title("Circularization radius as function of initial angular momentum")
figfile = "rcirc_l0.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()
