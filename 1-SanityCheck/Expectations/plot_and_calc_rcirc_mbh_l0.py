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

def mass_inside(rad,mbh,M_core,M_bulge,r_core,r_bulge):
	y = []
	for r in rad:
		if(r < r_core):
			y.append(mbh + M_core * (r / r_core)**(3))
		else:
			y.append(mbh + M_bulge * (r / r_bulge))
	return np.array(y)



Mbh = [1,1e-1,1e-2,1e-3,1e-4,0]
colors = ['r-', 'g-', 'b-', 'm-', 'c-', 'k.'] 
M_bulge = 1.0
M_core = 0.02
r_core = 0.02
r_bulge = 1.0
r_shell = 0.1
r_hole = 0.03
r_acc = 1e-3
rad = np.linspace(1e-6,1e-1,1025)
for col, mbh in zip(colors,Mbh):
	m_in = mass_inside(rad,mbh,M_core,M_bulge,r_core,r_bulge)
	v_eq = np.sqrt(m_in/rad)
	l_eq = v_eq * rad
	label_this="${\\rm M_{BH}} =$ "+str(mbh)
	imgplot = pl.loglog( l_eq, rad, col, linewidth = 0.5*FontSize, label = label_this)
#pl.axvspan(r_hole, r_shell, facecolor='g', alpha=0.5)

pl.legend(loc=4)
pl.xlabel("${l_0}$",fontsize=30)
pl.ylabel("$r_{\\rm circ}$",fontsize=30)
pl.xlim([1e-4,2e-2])
pl.ylim([1e-6,1e-1])
#pl.title("Circularization radius as function of black hole mass")
figfile = "rcirc_vs_l0.png"
pl.savefig(figfile)
print("Figure saved to "+figfile)
pl.close()
#pl.show()
