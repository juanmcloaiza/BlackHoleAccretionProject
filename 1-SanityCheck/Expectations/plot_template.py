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

datafile = "./lacc_mbh.txt"
print("reading " + datafile)
X = np.array(load_data(datafile))

label_this="$r_{\\rm acc} = 10^{-3}$"
imgplot = pl.loglog(X[:,1], X[:,0], label = label_this)

pl.legend(loc=2)
pl.xlabel("${\\rm M_{BH}}$")
pl.ylabel("$l_{\\rm acc}$")
#pl.xlim([0,1])
#pl.ylim([0,1])
pl.title("Minimum angular momentum to avoid accretion, $l_{\\rm acc}$")
figfile = "lacc_mbh.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()
