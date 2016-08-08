#!/usr/bin/env python
import matplotlib.pyplot as pl
import numpy as np
import sys

def load_data(filename):
	all_lines = []
	with open(filename) as f:
		for line in f:
			if (line[0] == '#'):
				continue
			#for char in line:
			try:
				line_data = np.array(line.split(),dtype=float)
			except ValueError:
				print("Line:")
				print(line)
				print("will give error.")
				print("Check the input file:")
				print(filename)
				print(" for possible corruption.")
			all_lines.append(line_data)
	return all_lines

def mag(string):
	if(string == 'M'):
		return 1e6
	elif(string == 'k'):
		return 1e3
	else:
		print("Error, string "+string+" cannot be converted to number.")
		sys.exit()
def prepare_figure():
	#Set relevant dimensions:
	FigSize = 10
	FontSize = 25
	BorderWidth = 3
	pl.rcParams.update({'font.size': FontSize})
	pl.rcParams.update({'axes.linewidth': BorderWidth})
	pl.rcParams.update({'lines.linewidth': 0.3*FontSize})
	pl.rcParams.update({'lines.markersize': 0.7*FontSize})

	#Open figure:
	pl.figure( figsize=( FigSize, FigSize ) )
	pl.tick_params(width=BorderWidth, length=FontSize, which='major')
	pl.tick_params(width=BorderWidth, length=0.3*FontSize, which='minor')

	pl.xlabel('$N_{\\rm resol}$')
	pl.ylabel('$M_{\\rm acc}/M_{\\rm shell}$')
	pl.xlim([0.9e5,2e7])
	pl.ylim([1e-6,5e-2])


def compute_acc_rate(Nfrac,Time):
	A = np.zeros((len(Nfrac),int(len(Nfrac)/2)))
	delta_t = np.zeros((len(Nfrac),int(len(Nfrac)/2)))
	Amean = np.zeros(len(Nfrac))
	dtmean = np.zeros(len(Nfrac))
	Astd = np.zeros(len(Nfrac))

#	count = 0
#	for delta_i in range(3,int(len(A)/2)):
#		for i in range(len(A)):
#			if(i+delta_i >= len(A)):
#				A[i,delta_i] = - ( Nfrac[-1] - Nfrac[i-delta_i] ) / ( Time[-1] - Time[i-delta_i] )
#				delta_t[i,delta_i] = ( Time[-1] - Time[i-delta_i] )
#			elif(i-delta_i < 0):
#				A[i,delta_i] = - ( Nfrac[i+delta_i] - Nfrac[0] ) / (Time[i+delta_i] - Time[0])
#				delta_t[i,delta_i] = (Time[i+delta_i] - Time[0])
#			else:
#				A[i,delta_i] = - ( Nfrac[i+delta_i] - Nfrac[i-delta_i] ) / (Time[i+delta_i] - Time[i-delta_i])
#				delta_t[i,delta_i] = (Time[i+delta_i] - Time[i-delta_i])
#
	delta_i = 1
	for i in range(len(A)):
		if( ( i - delta_i >= 0 ) and ( i + delta_i < len(A) )):
				Amean[i] = - ( Nfrac[i+delta_i] - Nfrac[i-delta_i] ) / (Time[i+delta_i] - Time[i-delta_i])
				dtmean[i] = (Time[i+delta_i] - Time[i-delta_i]) / (2.0*delta_i)
#		Amean[i] = np.mean(A[i,:])
#		dtmean[i] = np.mean(delta_t[i,:])
#		Astd[i] = np.std(A[i,:])

	return Amean, dtmean

def get_acc_info(datafile):
	#Read datafile
	DATA = np.array(load_data(datafile))
	Time = DATA[:,0] 
	Nfrac = (DATA[:,1]) / float(DATA[0,1])

	#Obtain resolution and accretion radius
	#from the filename:
	r_acc = float(datafile[5:9])
	resol = float(datafile[-8:-5]) * mag(datafile[-5])

	#Calculate Accretion Rate
	M_dot, dt = compute_acc_rate(Nfrac,Time)

	#Check that the integrated accretion 
	#rates are consistent with 
	#the total accreted mass:
	Tot_acc_mass = Nfrac[0] - Nfrac[-1]
	Int_acc_mass = sum(M_dot * dt)
	error = abs(Tot_acc_mass - Int_acc_mass)/Tot_acc_mass
	print("resol = "+str(resol)+" r_acc = "+str(r_acc)+", accreted mass: "+str(Tot_acc_mass)+", "+str(Int_acc_mass)+", error:"+str(error) )
	return resol, Tot_acc_mass, r_acc

def get_fit(Npart,M_measured):
	point_b = 2
	point_a = 1
	k = np.log(M_measured[point_b]/M_measured[point_a]) / np.log(Npart[point_a]/Npart[point_b])
	M_real = M_measured[-1]
	C0 = (M_measured[0] * Npart[0] **(k))
	x = np.linspace(Npart[0],Npart[-1],1e5)

	y = M_real + C0 * x**(-k)
	return x,y

def trace_the_line(X,Y,col):
	pl.loglog( X, Y, col+'o', label = "$r_{\\rm acc} ="+str(racc)+"$" )
	pl.legend(loc=3)
	x,y = get_fit(X,Y)
	pl.loglog(x,y, col+'--')
	#pl.loglog(x,y+Y[-1], col+'--')

	#for i in range(1,len(X)):
	#	i1 = i - 1
	#	annotation = ( (Y[i]/Y[i1]) ) / ( (X[i]/X[i1]) )
	#	pl.annotate( '$ {:3.2f} $'.format(annotation), xy=(X[i],Y[i]), xytext=(10,0), textcoords='offset points')
	#	print( '{:3.2f}'.format(annotation))


prepare_figure()

X = []
Y = []
for datafile in (
[
'Racc_5e-3_res_100k.txt',
'Racc_5e-3_res_250k.txt',
'Racc_5e-3_res_500k.txt',
'Racc_5e-3_res_750k.txt',
'Racc_5e-3_res_001M.txt',
'Racc_5e-3_res_004M.txt',
]):
	Nres, Macc, racc = get_acc_info(datafile)
	X.append(Nres)
	Y.append(Macc)
trace_the_line(X,Y,'b')


X = []
Y = []
for datafile in (
[
'Racc_1e-3_res_100k.txt',
'Racc_1e-3_res_250k.txt',
'Racc_1e-3_res_500k.txt',
'Racc_1e-3_res_750k.txt',
'Racc_1e-3_res_001M.txt',
'Racc_1e-3_res_002M.txt',
'Racc_1e-3_res_004M.txt',
]):
	Nres, Macc, racc = get_acc_info(datafile)
	X.append(Nres)
	Y.append(Macc)
trace_the_line(X,Y,'r')

X = []
Y = []
for datafile in (
[
'Racc_1e-4_res_100k.txt',
'Racc_1e-4_res_250k.txt',
'Racc_1e-4_res_500k.txt',
'Racc_1e-4_res_750k.txt',
'Racc_1e-4_res_001M.txt',
'Racc_1e-4_res_002M.txt',
'Racc_1e-4_res_004M.txt',
'Racc_1e-4_res_008M.txt',
'Racc_1e-4_res_010M.txt',
]):
	Nres, Macc, racc = get_acc_info(datafile)
	X.append(Nres)
	Y.append(Macc)
trace_the_line(X,Y,'g')



figfile = 'plot_acc_mass_vs_resol.png'
pl.savefig(figfile)
print('Figure saved to '+figfile)
