#!/usr/bin/env python3

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import numpy
import argparse

from lsl.common.stations import parse_ssmif


def main(args):
    #
    # Load in the data
    #
    site     = parse_ssmif(args.ssmif)
    dataFile = numpy.loadtxt(args.filename)
    
    #
    # Gather the station meta-data from its various sources
    #
    observer = site.get_observer()
    antennas = site.antennas
    
    #
    # Calculate the new stretch factors
    #
    output = [[None, None, None, None, None] for ant in antennas]
    for i in xrange(dataFile.shape[0]):
        ## Parse the line
        stand, ampX, addDelayX, ampY, addDelayY = dataFile[i,:]
        stand = int(stand)
        
        ## Pull out the cables
        digX, digY = None, None
        cableX, cableY = None, None
        for ant in antennas:
            if ant.stand.id == stand and ant.pol == 0:
                digX = ant.digitizer
                cableX = ant.cable
            elif ant.stand.id == stand and ant.pol == 1:
                digY = ant.digitizer
                cableY = ant.cable
                
        ## Find the current stretch factor/delay
        origX = cableX.stretch*1.0
        origY = cableY.stretch*1.0
        freq = numpy.linspace(35e6, 85e6, 101)
        baseDelayX = cableX.delay(freq, ns=True)
        baseDelayY = cableY.delay(freq, ns=True)
        
        bestX, bestY = 1e9, 1e9
        stretchX, stretchY = 0.90, 0.90
        for stretch in numpy.linspace(0.90, 1.10, 2001):
            cableX.stretch = stretch
            cableY.stretch = stretch
            
            newDelayX = cableX.delay(freq, ns=True)
            newDelayY = cableY.delay(freq, ns=True)
            
            diffX = (newDelayX - baseDelayX).mean() - addDelayX
            diffY = (newDelayY - baseDelayY).mean() - addDelayY
            
            if numpy.abs(diffX) < bestX:
                bestX = numpy.abs(diffX)
                stretchX = stretch
            if numpy.abs(diffY) < bestY:
                bestY = numpy.abs(diffY)
                stretchY = stretch
                
        output[digX-1] = [digX, stretchX, addDelayX, bestX, 9.0]
        output[digY-1] = [digY, stretchY, addDelayY, bestY, 9.0]
        
    #
    # Report
    #
    if args.output is not None:
        fh = open(args.output, 'w')
    else:
        fh = sys.stdout
        
    fh.write("##########################################\n")
    fh.write("#                                        #\n")
    fh.write("# Columns:                               #\n")
    fh.write("# 1) Digitizer                           #\n")
    fh.write("# 2) Cable length stretch factor         #\n")
    fh.write("# 3) Additional delay over previous (ns) #\n")
    fh.write("# 4) Mean delay error in new stretch     #\n")
    fh.write("# 5) The number 9                        #\n")
    fh.write("#                                        #\n")
    fh.write("##########################################\n")
    for entry in output:
        fh.write("%3i  %6.4f  %6.4f  %6.4f  %6.4f\n" % tuple(entry))
        
    if args.output is not None:
        fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="given a file containing additional delays over the current SSMIF, convert the delays to new stretch factors",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('ssmif', type=str,
                        help='station SSMIF')
    parser.add_argument('filename', type=str, 
                        help='filename to convert')
    parser.add_argument('-o', '--output', type=str,
                        help='write output to the specified filename instead of the screen')
    args = parser.parse_args()
    main(args)
    