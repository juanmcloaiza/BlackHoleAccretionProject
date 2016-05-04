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

def get_rcirc(radius,nparts):
	mean_r = 0
	for r,n in zip(radius,nparts):
		mean_r += r*n
	return mean_r

def correct_mass(x,error):
	if (x == 0):
		x = error
	else:
		x = x * (1.0-error) + error
	return x
	

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


snap=5

#Plot the simulations with accretion radius correcting for erros
#in mass calculation:
for res, error in zip(['100k','250k','500k','750k','001M'],[0.00147026, 0.000535889, 0.000154021, 7.73302e-05, 5.09935e-05]):
	datafile = "./res_"+res+"_t"+str(snap)+"00.txt"
	X = np.array(load_data(datafile))
	M_cumulative = np.array(X[:,1])
#	print('#'+res)
	for i in range(1,len(M_cumulative)):
		M_cumulative[i] += M_cumulative[i-1]
	for i in range(1,len(M_cumulative)):
		M_cumulative[i] = correct_mass(M_cumulative[i],error)
#		print(X[i,0], M_cumulative[i])
	label_this= res#+" particles, $t = "+str(snap)+"$"
	imgplot = pl.loglog(X[:,0], M_cumulative,'--', linewidth = 0.2*FontSize, label = label_this)

#Plot the same simulations without correction for mass calculation:
for res, color in zip(['100k','250k','500k','750k','001M'],['b.','g.','r.','c.','m.']):
	datafile = "./res_"+res+"_t"+str(snap)+"00.txt"
	X = np.array(load_data(datafile))
	M_cumulative = np.array(X[:,1])
	for i in range(1,len(M_cumulative)):
		M_cumulative[i] += M_cumulative[i-1]
	imgplot = pl.loglog(X[:,0], M_cumulative,color, linewidth = 0.2*FontSize)

#Plot the simulations without accretion radius:
for res, color in zip(['100k','250k','500k','750k','001M'],['b-','g-','r-','c-','m-']):
	datafile = "../../No_racc/M_inside_r/res_"+res+"_t"+str(snap)+"00.txt"
	X = np.array(load_data(datafile))
	M_cumulative = np.array(X[:,1])
	for i in range(1,len(M_cumulative)):
		M_cumulative[i] += M_cumulative[i-1]
	imgplot = pl.loglog(X[:,0], M_cumulative,color, linewidth = 0.2*FontSize)


pl.legend(loc=2)
pl.xlabel("$r$",fontsize=1.5*FontSize)
pl.ylabel("$N(r)/N_{\\rm part}$",fontsize=1.5*FontSize)
pl.xlim([1e-4,3e-2])
pl.ylim([1e-6,1e0])
#pl.title("Circularization radius is where the peak is..." )
figfile = "M_inside_r_all_res.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()
