#!/usr/bin/env python3

"""
Given a /data/temp.txt file (or one of the rotated backups) plot the temperatures
of FPGAs in DP as a function of time, chassis, and physical slot.  These 
temperatures are plotted as color maps for (1) the mean FPGA temperature per board
and (2) the maximum FPGA temperature per board.
"""

import sys
import numpy
import pytz
from datetime import datetime

from matplotlib import pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.dates import *
from matplotlib.ticker import *

MST7MDT = pytz.timezone('US/Mountain')


# Logical to physical mapping
LPmapping = { 1:  7, 
              2:  8, 
              3:  6, 
              4:  9, 
              5:  5, 
              6: 10, 
              7:  4, 
              8: 11, 
              9:  3, 
             10: 12, 
             11:  2, 
             12: 13, 
             13:  1, 
             14: 14,}

def main(args):
    if len(args) < 1:
        print('Need a filename to plot.')
        sys.exit(1)
    
    data = []
    for filename in args:
        # Read in the data
        fh = open(filename)
        for line in fh:
            # Split the line by commas
            line = line.replace('\n', '')
            fields = line.split(',')
            if fields[-1] == '':
                fields = fields[:-1]
                
            # Skip over lines that are probably wrong
            if len(fields) > 141:
                print("WARNING: Entry at %s has too many chips" % MST7MDT.localize(datetime.fromtimestamp(float(fields[0]))))
                continue
            
            # Pad if we find less chips than expected and emit a warning
            if len(fields) < 141:
                print("WARNING: Entry at %s has only %i chips" % (MST7MDT.localize(datetime.fromtimestamp(float(fields[0]))), len(fields)-1))
                fields.extend(['0.0']*(141-len(fields)))
            
            #print(len(fields))
            # Convert all values to floats
            data.append( [float(f) for f in fields] )

    # Split out the time and interpret it
    data = numpy.array(data)
    order = numpy.argsort(data[:,0])
    data = data[order,:]
    
    dates = [MST7MDT.localize(datetime.fromtimestamp(t)) for t in data[:,0]]
    print('File spans %s to %s with %i measurements' % (dates[0], dates[-1], len(dates)))
    tRange = (data[-1,0] - data[0,0]) / 3600.0
    
    # Split into boards
    boardMean = numpy.zeros((data.shape[0], 28))
    boardMax  = numpy.zeros((data.shape[0], 28))
    for i in range(data.shape[0]):
        line = data[i,1:]
        line.shape = (28, 5)
        
        for j in range(28):
            boardMean[i,j] = line[j,:].mean()
            boardMax[i,j]  = line[j,:].max()
    
    tempMean1 = numpy.zeros((data.shape[0], 14))
    tempMax1  = numpy.zeros((data.shape[0], 14))
    tempMean2 = numpy.zeros((data.shape[0], 14))
    tempMax2  = numpy.zeros((data.shape[0], 14))
    for i in range(boardMean.shape[1]):
        if i < 14:
            # lwa16 is board #1 and it is in the bottom chassis
            p = LPmapping[i+1]
            tempMean1[:,p-1] = boardMean[:,i]
            tempMax1[:,p-1] = boardMax[:,i-1]
        else:
            # lwa34 is board #15 and it is in the top chassis
            p = LPmapping[i-13]
            tempMean2[:,p-1] = boardMean[:,i]
            tempMax2[:,p-1] = boardMax[:,i]
    
    valid = numpy.where( boardMean > 0 )
    
    # Plot all of the Chips
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    
    vmin = boardMean[valid].min()
    vmax = boardMean[valid].max()
    
    c = ax1.imshow(tempMean2, origin='lower', extent=(0.5, 14.5, 0, tRange), vmin=vmin, vmax=vmax, interpolation='nearest')
    cb = fig.colorbar(c, ax=ax1)
    cb.ax.set_ylabel('Temp. [C]')
    
    c = ax2.imshow(tempMean1, origin='lower', extent=(0.5, 14.5, 0, tRange), vmin=vmin, vmax=vmax, interpolation='nearest')
    cb = fig.colorbar(c, ax=ax2)
    cb.ax.set_ylabel('Temp. [C]')
    
    ax1.axis('auto')
    ax1.set_xticks(range(1, 15) )
    ax2.axis('auto')
    ax2.set_xticks(range(1, 15) )
    
    fig.suptitle('Mean FPGA Temperatures Per Board')
    ax2.set_xlabel('Chassis Slot')
    ax1.set_ylabel('%s' % dates[0].strftime('%m/%d/%Y %H:%M:%S %Z'))
    ax2.set_ylabel('Hours Since')
    
    plt.draw()
    
    # Plot all of the Chips
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    
    vmin = boardMax[valid].min()
    vmax = boardMax[valid].max()
    
    c = ax1.imshow(tempMax2, origin='lower', extent=(0.5, 14.5, 0, tRange), vmin=vmin, vmax=vmax, interpolation='nearest')
    cb = fig.colorbar(c, ax=ax1)
    cb.ax.set_ylabel('Temp. [C]')
    
    c = ax2.imshow(tempMax1, origin='lower', extent=(0.5, 14.5, 0, tRange), vmin=vmin, vmax=vmax, interpolation='nearest')
    cb = fig.colorbar(c, ax=ax2)
    cb.ax.set_ylabel('Temp. [C]')
    
    ax1.axis('auto')
    ax1.set_xticks(range(1, 15) )
    ax2.axis('auto')
    ax2.set_xticks(range(1, 15) )
    
    fig.suptitle('Maximum FPGA Temperatures Per Board')
    ax2.set_xlabel('Chassis Slot')
    ax1.set_ylabel('%s' % dates[0].strftime('%m/%d/%Y %H:%M:%S %Z'))
    ax2.set_ylabel('Hours Since')
    
    plt.draw()
    
    # Show
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
