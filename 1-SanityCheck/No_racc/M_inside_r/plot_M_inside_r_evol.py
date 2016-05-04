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


res='001M'
#print("#M_bh r_circ")
for snap in ['000','010','050','100','200','300','500']:
	datafile = "./res_"+res+"_t"+snap+".txt"
	X = np.array(load_data(datafile))
	M_cumulative = np.array(X[:,1])
	for i in range(1,len(M_cumulative)):
		M_cumulative[i] += M_cumulative[i-1]
	label_this= snap#+" particles, $t = "+str(snap)+"$"
	imgplot = pl.loglog(X[:,0], M_cumulative,'-', linewidth = 0.5*FontSize, label = label_this)
	#imgplot2 = pl.loglog(X[:,0], X[:,1],'--', linewidth = 0.5*FontSize, label = label_this)
	#print("1e"+str(i-10)+" "+str(get_rcirc(X[:,0],X[:,1])))

#plot expectation: From eq 12:
#M(<r*) / Mshell ~ Mbh r* / rshell^2 / vrot^2

rad = np.linspace(1e-3,1e-1)
m_in = 1e-2 * rad / 0.1 / 0.3**2
pl.plot(rad,m_in,'k--',label='expected slope')


pl.legend(loc=2)
pl.xlabel("$r$",fontsize=1.5*FontSize)
pl.ylabel("$N_{\\rm part}$",fontsize=1.5*FontSize)
pl.xlim([2e-3,1e-1])
pl.ylim([1e-3,1e0])
#pl.title("Circularization radius is where the peak is..." )
figfile = "M_inside_r_evol.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
#pl.show()
