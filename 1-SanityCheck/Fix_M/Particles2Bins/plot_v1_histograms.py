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

for i in [100,200,300,400,488]:
	datafile = "./v1_r_l_histogram_t"+str(i)+".txt"
	X = np.array(load_data(datafile))
	label_this="$v_{\\rm rot} = 0.1, t =$ "+str(i)
	imgplot = pl.plot(X[:,0], X[:,1],  label = label_this)

pl.legend(loc=1)
pl.xlabel("$r$")
pl.ylabel("$N_{\\rm part}$")
pl.xlim([0,0.01])
#pl.ylim([0,0.30])
pl.title("Circularization radius is where the peak is..." )
figfile = "v1_histograms.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()
