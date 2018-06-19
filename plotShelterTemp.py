#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plotShelterTemp.py - Script to read in the shelter.txt file and plot up the 
shelter temperature (in F) as a function of time.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
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
        print 'Need a filename to plot.'
        sys.exit(1)
    
    data = []
    inF = True
    for filename in args:
        # Read in the data
        fh = open(filename)
        for line in fh:
            line = line.replace('\n', '')
            fields = line.split()
            if len(fields) == 1:
                # What's this for?  There are two types of shelter logs:  Ones
                # from polling MCS and writing to a file and ones produced by
                # SHL MCS as part of its operation.  The first kind is space
                # separated and the second is comma separated.  The second kind
                # is also degrees C, rather than degrees F.
                inF = False
                fields = line.split(',')
            fields = [float(f) for f in fields]
            data.append( fields )

    # Split out the time and interpret it
    data = numpy.array(data)
    order = numpy.argsort(data[:,0])
    data = data[order,:]
    
    dates = [MST7MDT.localize(datetime.fromtimestamp(t)) for t in data[:,0]]
    print 'File spans %s to %s with %i measurements' % (dates[0], dates[-1], len(dates))
    
    # Convert to degree F if needed
    if not inF:
        data[:,1] = data[:,1] * 9 / 5 + 32
    
    # Plot
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot_date(dates, data[:,1], fmt='-', tz=MST7MDT, marker='x', linestyle=' ')

    # Label and format dates
    ax1.set_title('Shelter Temperature')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Temperature [$^\circ$F]')
    fig.autofmt_xdate()

    # Show
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
