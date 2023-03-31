#!/usr/bin/env python3

"""
Given a /data/temp.txt file (or one of the rotated backups) plot the temperatures
of all 140 FPGAs in DP.
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
    
    # Plot all of the Chips
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    for i in range(1,141):
        ax1.plot_date(dates, data[:,i], fmt='-', tz=MST7MDT, marker='x', linestyle=' ')

    # Label and format dates
    ax1.set_title('DP FPGA Temperatures')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Temperature [$^\circ$C]')
    fig.autofmt_xdate()

    # Show
    ax1.set_ylim(-5, 95)
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
