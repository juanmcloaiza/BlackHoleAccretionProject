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

snap=500
for i in [6,7,8,9]:
	datafile = "./M"+str(i)+"_r_l_histogram_t"+str(snap)+".txt"
	X = np.array(load_data(datafile))
	label_this="$M_{\\rm BH} = 10^"+str(i)+" M_{\\odot}, t = "+str(snap)+"$"
	imgplot = pl.plot(X[:,0], X[:,1],  label = label_this)

pl.legend(loc=2)
pl.xlabel("$r$")
pl.ylabel("$N_{\\rm part}$")
#pl.xlim([0,0.03])
#pl.ylim([0,0.30])
pl.title("Circularization radius is where the peak is..." )
figfile = "Rcirc.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()
