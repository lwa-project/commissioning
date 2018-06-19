#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import getopt

from lsl.common.stations import parseSSMIF


def usage(exitCode=None):
    print """convertDelayToStretch.py - Given a file containing additional delays 
over the current SSMIF, convert the delays to new stretch factors.  These 
stretch factors can then be used the applyNewStretchFactors.py to generate an
updated SSMIF.

Valid sources for additional delay files are applySelfCalTBW2.py and 
checkConsistency.py.

Usage: convertDelayToStretch.py [OPTIONS] SSMIF delayFile

Options:
-h, --help             Display this help information
-o, --output           Write output to the specified filename 
                    (default = write to screen)
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseConfig(args):
    config = {}
    # Command line flags - default values
    config['outname'] = None

    # Read in and process the command line flags
    try:
        opts, arg = getopt.getopt(args, "ho:", ["help", "output="])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-o', '--output'):
            config['outname'] = value
        else:
            assert False
            
    # Add in arguments
    config['args'] = arg
    
    # Validate the input
    if len(config['args']) != 2:
        raise RuntimeError("Must provide both a SSMIF and delay file")
        
    # Return configuration
    return config


def main(args):
    #
    # Parse the command line
    #
    config = parseConfig(args)
    
    #
    # Load in the data
    #
    site     = parseSSMIF(config['args'][0])
    dataFile = numpy.loadtxt(config['args'][1])
    
    #
    # Gather the station meta-data from its various sources
    #
    observer = site.getObserver()
    antennas = site.getAntennas()
    
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
    if config['outname'] is not None:
        fh = open(config['outname'], 'w')
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
        
    if config['outname'] is not None:
        fh.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    