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

def add_point(datafile,style,labelswitch=True):
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

	pl.loglog( resol, Tot_acc_mass, style, label = "$r_{\\rm acc} ="+str(r_acc)+"$" )
	if(labelswitch):
		pl.legend()


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
	add_point(datafile,'ro',False)

for datafile in (
[
'Racc_1e-4_res_100k.txt',
'Racc_1e-4_res_250k.txt',
'Racc_1e-4_res_500k.txt',
'Racc_1e-4_res_750k.txt',
'Racc_1e-4_res_001M.txt',
'Racc_1e-4_res_004M.txt',
]):
	add_point(datafile,'bo',False)

for datafile in (
[
'Racc_5e-3_res_100k.txt',
'Racc_5e-3_res_250k.txt',
'Racc_5e-3_res_500k.txt',
'Racc_5e-3_res_750k.txt',
'Racc_5e-3_res_001M.txt',
]):
	add_point(datafile,'go',False)


pl.show()
