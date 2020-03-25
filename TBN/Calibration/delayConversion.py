#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Combine the delay differences from a text file with the a priori knowledge of the
cable model to make a NPZ file that reflects the fully delay of the system.
"""

import sys
import numpy
from matplotlib import pyplot as plt

from lsl.common.stations import lwa1


def main(args):
    antennas = lwa1.getAntennas()
    
    try:
        filename = args[0]
    except:
        filename = 'add-delay.txt'

    fh = open(filename)

    delays = {}
    for i in xrange(520):
        delays[i] = []

    for line in fh:
        line = line.replace('\n', '')

        parts = line.split(None)
        ant = parts[0]
        delay = parts[1]
        ant = int(ant)
        delay = float(delay)
        
        if len(parts) == 3:
            chi2 = float(parts[2])
            if chi2 > 1.0:
                print "Skipping %3i due to large (%6.3f) Chi2 [status=%i]" % (ant, chi2, antennas[ant].getStatus())
                continue
        
        delays[ant].append( delay )

    fh.close()

    pols = numpy.array([i % 2 for i in xrange(520)])

    delays2 = numpy.zeros(520)
    delays3 = numpy.zeros(520) - 1000.0
    for i in xrange(520):
        part = numpy.array(delays[i])
        if part.size == 0:
            continue
        delays2[i] = numpy.median(part)
        delays3[i] = numpy.median(part)
        print i, delays2[i], part.std()

    #suspect = numpy.where( numpy.abs(delays3) > 35e-9 )[0]
    #for i in suspect:
        #print i, delays3[i]*1e9

    fig = plt.figure()
    ax = fig.gca()
    goodX = numpy.where( (delays3 != -1000) & (pols == 0) )
    goodY = numpy.where( (delays3 != -1000) & (pols == 1) )
    ax.hist(delays3[goodX]*1e9, bins=50, label='X')
    ax.hist(delays3[goodY]*1e9, bins=50, label='Y')
    ax.legend(loc=0)
    ax.set_xlabel('Residual Delay [ns]')
    ax.set_ylabel('Number of Antennas')
    plt.show()


    numpy.savez('add-delay.npz', delayDiffs=delays2)


if __name__ == "__main__":
    main(sys.argv[1:])
