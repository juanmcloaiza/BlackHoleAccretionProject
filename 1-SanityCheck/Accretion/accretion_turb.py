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

    pl.xlabel('$t$')
    pl.ylabel('$M_{\\rm acc}/M_{\\rm shell}( \\times 100)$')
    #pl.xlim([0.9e5,2e7])
    #pl.ylim([1e-6,5e-2])



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
    pl.plot(x,y,'b--')
    #pl.semilogx(x,y,'b--')


def add_plot(datafile,style,label="."):
    #Read datafile
    DATA = np.array(load_data(datafile))
    Time = DATA[:,0] 
    Nfrac = (DATA[:,1]) / float(DATA[0,1])
    dt = np.diff(Time)
    dN = np.diff(1.0 - Nfrac)
    
    x = Time
    y = 100*(1 - Nfrac)

    pl.plot( x, y, style, label=label)
    pl.plot( x[-1], y[-1], style[0]+'*', label=None, markersize='50')
    #pl.loglog( x, y, style, label=".")
    #pl.semilogx( x, y, style, label=label)
    #pl.semilogy( x, y, style, label=".")
    #add_fit(x,y)

    annotation = y[-1]
    pl.annotate( '${:3.2f}$'.format(annotation), xy=(x[-1],y[-1]), xytext=(-60,10), textcoords='offset points')
    pl.legend(loc=2,framealpha=0.7)

prepare_figure()

#for datafile in (
#[
#'Racc_1e-3_res_100k.txt',
#'Racc_1e-3_res_250k.txt',
#'Racc_1e-3_res_500k.txt',
#'Racc_1e-3_res_750k.txt',
#'Racc_1e-3_res_001M.txt',
#'Racc_1e-3_res_002M.txt',
#'Racc_1e-3_res_004M.txt',
#]):
#    add_plot(datafile,'k-',False)

#for datafile in (
#[
#'Racc_1e-4_res_100k.txt',
#'Racc_1e-4_res_250k.txt',
#'Racc_1e-4_res_500k.txt',
#'Racc_1e-4_res_750k.txt',
#'Racc_1e-4_res_001M.txt',
#'Racc_1e-4_res_002M.txt',
#'Racc_1e-4_res_004M.txt',
#'Racc_1e-4_res_008M.txt',
#'Racc_1e-4_res_010M.txt',
#]):
#    add_plot(datafile,'g-',False)

#for datafile in (
#[
#'Racc_5e-3_res_100k.txt',
#'Racc_5e-3_res_250k.txt',
#'Racc_5e-3_res_500k.txt',
#'Racc_5e-3_res_750k.txt',
#'Racc_5e-3_res_001M.txt',
#]):
#    add_plot(datafile,'b-',False)
for datafile in ([
    'TurbGrid_128.txt',
]):
    add_plot(datafile,'b--','$128$')


for datafile in ([
    'TurbGrid_256.txt',
]):
    add_plot(datafile,'r--','$256$')

for datafile in ([
    'TurbGrid_512.txt',
]):
    add_plot(datafile,'g--','$512$')

#pl.show()
figfile = 'plot_accretion_turb.png'
pl.savefig(figfile)
print('Figure saved to '+figfile)
