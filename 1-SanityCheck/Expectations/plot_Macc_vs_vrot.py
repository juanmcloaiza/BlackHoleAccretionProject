#!/usr/bin/env python
import matplotlib.pyplot as pl
import numpy as np

def load_data(filename):
	all_lines = []
	nfields=2
	with open(filename) as f:
		for line in f:
			if (line[0] == '#'):
				continue
			#for char in line:
			line_data = np.array(line.split(),dtype=float)
			all_lines.append(line_data)
	return all_lines

#datafile = ".//output_04threads.txt"
#X = np.array(load_data(datafile))
#label_this="$N_{\\rm threads} = 4$"
#imgplot = pl.loglog(X[:,0], X[:,1], 'ro-', label = label_this)

racc = 1e-3
rcore = 20 * racc
rshell = 100 * racc

Mcore = 0.02
Mbh = lambda x: 0.01 * x

vrot = np.linspace(0.1,1.0,10)
print(vrot)
Nacc = (Mcore * racc * rshell / vrot**2) * ( (Mbh(10.0)/Mcore) + (racc/rcore)**3 )
pl.loglog(vrot,Nacc,'r')

Nacc = (Mcore * racc * rshell / vrot**2) * ( (Mbh(1.0)/Mcore) + (racc/rcore)**3 )
pl.loglog(vrot,Nacc,'g')

Nacc = (Mcore * racc * rshell / vrot**2) * ( (Mbh(0.1)/Mcore) + (racc/rcore)**3 )
pl.loglog(vrot,Nacc,'b')

Nacc = (Mcore * racc * rshell / vrot**2) * ( (Mbh(0.01)/Mcore) + (racc/rcore)**3 )
pl.loglog(vrot,Nacc,'m')
#pl.legend(loc=2)
pl.xlabel("$v_{\\rm rot}$")
pl.ylabel("$t \\, \\, [s]$")
pl.xlim([0,1])
#pl.ylim([0,1])
pl.title("Accreted masses as function of initial $v_{\\rm rot}$")
figfile = "Macc_vs_vrot.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()

