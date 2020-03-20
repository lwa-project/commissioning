#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script to read in the positions of stands at LWA-1 and make a plot
of the site."""

import os
import sys
import numpy
import getopt

from lsl.common import stations
from lsl.correlator import uvutil

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def usage(exitCode=None):
    print """findMyStand.py - Plot the x, y, and z locations of stands at 
LWA-1 and mark and label particular stands.

Usage: findMyStand.py [OPTIONS] stand1 [stand2 [...]]]

Options:
-h, --help             Display this help information
-v, --verbose          Run plotStands in verbose mode
"""

    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    
# Command line flags - default values
    config['label'] = False
    config['verbose'] = False
    config['args'] = []

    # Read in and process the command line flags
    try:
        opts, arg = getopt.getopt(args, "hv", ["help", "verbose"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-v', '--verbose'):
            config['verbose'] = True
        else:
            assert False
    
    # Add in arguments
    config['args'] = [int(i) for i in arg]

    # Return configuration
    return config


def main(args):
    # Parse command line
    config = parseOptions(args)
    toMark = numpy.array(config['args'])-1
    
    # Set the LWA Station
    station = stations.lwa1
    stands = station.getStands()
    stands.sort()

    # Load in the stand position data
    data = numpy.zeros((len(stands),3))
    
    i = 0
    for stand in stands[::2]:
        if stand.x == 0 and stand.y == 0 and stand.z == 0:
            continue
        data[i,0] = stand.x
        data[i,1] = stand.y
        data[i,2] = stand.z
        i += 1
    data = data[:i,:]

    # Color-code the stands by their elevation
    color = data[:,2]

    # Plot the stands as colored circles
    fig = plt.figure(figsize=(8,8))

    ax1 = fig.gca()
    c = ax1.scatter(data[:,0], data[:,1], c=color, s=40.0, alpha=0.30)
    ax1.set_xlabel('$\Delta$X [E-W; m]')
    ax1.set_xlim([-80, 80])
    ax1.set_ylabel('$\Delta$Y [N-S; m]')
    ax1.set_ylim([-80, 80])
    ax1.set_title('%s Site:  %.3f$^\circ$N, %.3f$^\circ$W' % (station.name, station.lat*180.0/numpy.pi, -station.long*180.0/numpy.pi))
    
    # Explicitly mark those that need to be marked
    if toMark.size != 0:
        for i in xrange(toMark.size):
            ax1.plot(data[toMark[i],0], data[toMark[i],1], marker='x', linestyle='x', markersize=10, color='black')
            
            ax1.annotate('%i' % (toMark[i]+1), xy=(data[toMark[i],0], data[toMark[i],1]), xytext=(data[toMark[i],0]+1, data[toMark[i],1]+1))

    # Set the axis limits
    ax1.set_xlim([-60, 60])
    ax1.set_ylim([-60, 60])

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
