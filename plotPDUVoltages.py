#!/usr/bin/env python3

"""
Given a /data/rack##.txt file (or one of the rotated backups) plot the PDU input
voltage over time.  This script is designed to accept multiple files from multiple
rack if needed.
"""

# Python2 compatiability
from __future__ import print_function, division

import os
import re
import sys
import numpy
import pytz
from datetime import datetime

from matplotlib import pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.dates import *
from matplotlib.ticker import *

MST7MDT = pytz.timezone('US/Mountain')

filenameRE = re.compile(r'rack(?P<rack>\d{1,2}).txt')


def main(args):
    if len(args) < 1:
        print('Need a filename to plot.')
        sys.exit(1)
    
    data = {}
    for filename in args:
        mtch = filenameRE.match(os.path.split(filename)[1])
        rack = int(mtch.group('rack'))
        
        # Read in the data
        fh = open(filename)
        for line in fh:
            # Split the line by commas
            line = line.replace('\n', '')
            fields = line.split(',')
            if fields[-1] == '':
                fields = fields[:-1]
            
            try:
                data[rack].append( [float(f) for f in fields] )
            except KeyError:
                data[rack] = []
                data[rack].append( [float(f) for f in fields] )

    # Split out the time and interpret it
    dates = {}
    for i,k in enumerate(data.keys()):
        data[k] = numpy.array(data[k])
        order = numpy.argsort(data[k][:,0])
        data[k] = data[k][order,:]
    
        dates[k] = [MST7MDT.localize(datetime.fromtimestamp(t)) for t in data[k][:,0]]
        if i == 0:
            print('File spans %s to %s with %i measurements' % (dates[k][0], dates[k][-1], len(dates[k])))
    
    # Plot all of the inputs
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    for k in dates.keys():
        ax1.plot_date(dates[k], data[k][:,2], fmt='-', tz=MST7MDT, marker='x', linestyle=' ', label='Rack #%i' % k)

    # Label and format dates
    ax1.legend(loc=0)
    ax1.set_title('PDU Input Voltage')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Voltage [VAC]')
    fig.autofmt_xdate()

    # Show
    #if data[:,2].mean() > 160:
        #ax1.set_ylim(192, 248)
    #else:
        #ax1.set_ylim(96, 144)
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
