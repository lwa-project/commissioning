#!/usr/bin/env python3

"""
Example script to read in the positions of stands at LWA-1 and make a plot
of the site.
"""

import os
import sys
import numpy
import argparse

from lsl.common import stations
from lsl.correlator import uvutils
from lsl.misc import parser as aph

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def main(args):
    # Parse command line
    toMark = numpy.array(args.stand) - 1
    
    # Set the LWA Station
    station = stations.lwa1
    stands = station.stands
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
        for i in range(toMark.size):
            ax1.plot(data[toMark[i],0], data[toMark[i],1], marker='x', linestyle='', markersize=10, color='black')
            
            ax1.annotate('%i' % (toMark[i]+1), xy=(data[toMark[i],0], data[toMark[i],1]), xytext=(data[toMark[i],0]+1, data[toMark[i],1]+1))

    # Set the axis limits
    ax1.set_xlim([-60, 60])
    ax1.set_ylim([-60, 60])

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='plot the x, y, and z locations of stands at LWA-1 and mark and label particular stands',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('stand', type=aph.positive_int, nargs='+',
                        help='stand ID to plot')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='run %(prog)s in verbose mode')
    args = parser.parse_args()
    main(args)
    
