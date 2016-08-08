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

def add_fit(datax,datay):
    a = 15
    b = -15
    A = ( datay[b] - datay[a] ) / ( np.log( datax[b] ) - np.log( datax[a] ) )
    K = datay[1] - A * np.log(datax[1])
    print("accreted mass = {}, fit = {} log(t) + {}".format(datay[-1],A,K))
   
    x = np.linspace(min(datax),max(datax),100)
    y = A * np.log(x) + K
    #pl.plot(x,y,'b--')
    pl.semilogx(x,y,'b--')


def add_plot(datafile,style,labelswitch=True):
    #Read datafile
    DATA = np.array(load_data(datafile))
    Time = DATA[:,0] 
    Npart = (DATA[:,1])
    
    x = Time
    y = Npart[0] - Npart

    #pl.plot( x, y, style, label=".")
    pl.semilogx( x, y, style, label=".")
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
    add_plot(datafile,'ro',False)

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
    add_plot(datafile,'go',False)


#for datafile in (
#[
#'Racc_5e-3_res_100k.txt',
#'Racc_5e-3_res_250k.txt',
#'Racc_5e-3_res_500k.txt',
#'Racc_5e-3_res_750k.txt',
#'Racc_5e-3_res_001M.txt',
#]):
#    add_plot(datafile,'bo',False)

for datafile in (
[
'TurbGrid_256.txt',
'TurbGrid_512.txt',
]):
    add_plot(datafile,'go',False)

pl.show()
