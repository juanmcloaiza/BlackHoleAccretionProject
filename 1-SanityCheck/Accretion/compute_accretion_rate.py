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

for datafile in (
[
'Racc_1e-3_res_100k.txt',
#'Racc_1e-3_res_250k.txt',
#'Racc_1e-3_res_500k.txt',
#'Racc_1e-3_res_750k.txt',
'Racc_1e-3_res_001M.txt',
'Racc_1e-4_res_100k.txt',
#'Racc_1e-4_res_250k.txt',
#'Racc_1e-4_res_500k.txt',
#'Racc_1e-4_res_750k.txt',
'Racc_1e-4_res_001M.txt',
#'Racc_1e-4_res_004M.txt',
'Racc_5e-3_res_100k.txt',
#'Racc_5e-3_res_250k.txt',
#'Racc_5e-3_res_500k.txt',
#'Racc_5e-3_res_750k.txt',
'Racc_5e-3_res_001M.txt',
]):


	DATA = np.array(load_data(datafile))
	Time = DATA[:,0] 
	Npart = (DATA[:,1]) / float(DATA[0,1])
	A = np.zeros((len(Npart),int(len(Npart)/2)))
	Amean = np.zeros(len(Npart))
	Astd = np.zeros(len(Npart))

	count = 0
	for delta_i in range(3,int(len(A)/2)):
		for i in range(len(A)):
			if(i+delta_i >= len(A)):
				A[i,delta_i] = - ( Npart[-1] - Npart[i-delta_i] ) / ( Time[-1] - Time[i-delta_i] )
			elif(i-delta_i < 0):
				A[i,delta_i] = - ( Npart[i+delta_i] - Npart[-1] ) / (Time[i+delta_i] - Time[-1])
			else:
				A[i,delta_i] += - ( Npart[i+delta_i] - Npart[i-delta_i] ) / (Time[i+delta_i] - Time[i-delta_i])

	for i in range(len(A)):
		Amean[i] = np.mean(A[i,:])
		Astd[i] = np.std(A[i,:])

	#pl.plot(Time,Npart)
	pl.semilogy(Time,Amean,label = datafile)
	pl.legend()
	#pl.plot(Time,Astd,'b-')
	#pl.plot(Time,Amean-Astd,'r--')
	#pl.plot(Time,Amean+Astd,'r--')
	#pl.ylim([Npart[-1],Npart[0]])
pl.show()
