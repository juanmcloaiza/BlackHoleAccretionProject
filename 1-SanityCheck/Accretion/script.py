#!/usr/bin/env python
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

for RACC in ['1e-3']:
	for RES in ['002M','004M']:
		datafile = '/scratch/jcarmona/Racc_'+RACC+'_res_'+RES+'/Output/Accretion.txt'
		print("loading "+datafile)
		DATA = np.array(load_data(datafile))
		
		outfilename = 'Racc_'+RACC+'_res_'+RES+'.txt'
		outfile = open(outfilename,"w")
		for i in range(1,len(DATA)):
			if( DATA[i,1] != DATA[i-1,1] ):
				outfile.write(str(DATA[i,0])+" "+str(DATA[i,1])+" "+str(DATA[i,2])+"\n")
		outfile.close()
