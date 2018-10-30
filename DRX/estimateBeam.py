#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read in a SSMIF file and estimate the DRX beam for a given frequency and 
topocentric pointing center.  The estimate is based off a simple delay-and-sum 
beam former so it won't be an exact match.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import sys
import time
import numpy
import getopt

from lsl.misc import beamformer
from lsl.common.stations import parseSSMIF

from matplotlib import pyplot as plt

def usage(exitCode=None):
    print """estimateBeam.py - Read in a SSMIF file and estimate the DRX beam 
for a given frequency and topocentric pointing center.  The estimate is based 
off a simple delay-and-sum beam former so it won't be an exact match.

Usage: estimateBeam.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-f, --frequency             Frequency in MHz to calculate the beam for 
                            (Default = 65 MHz)
-a, --azimuth               Azimuth east of north in degrees for the pointing center
                            (Default = 90 degrees)
-e, --elevation             Elevation above the horizon in degrees for the pointing
                            center (Default = 90 degrees)
-n, --no-plots              Do not show plots of the beam response
"""
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['freq'] = 65.0e6
    config['az'] = 90.0
    config['el'] = 90.0
    config['plots'] = True
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hf:a:e:n", ["help", "frequency=", "azimuth=", "elevation=", "no-plots"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-f', '--frequency'):
            config['freq'] = float(value)*1e6
        elif opt in ('-a', '--azimuth'):
            config['az'] = float(value)
        elif opt in ('-e', '--elevation'):
            config['el'] = float(value)
        elif opt in ('-n', '--no-plots'):
            config['plots'] = False
        else:
            assert False
    
    # Add in arguments
    config['args'] = args

    # Return configuration
    return config


def main(args):
    config = parseOptions(args)
    filename = config['args'][0]

    station = parseSSMIF(filename)
    antennas = station.getAntennas()

    digs    = numpy.array([ant.digitizer  for ant in antennas])
    ants    = numpy.array([ant.id         for ant in antennas])
    stands  = numpy.array([ant.stand.id   for ant in antennas])
    pols    = numpy.array([ant.pol        for ant in antennas])
    antStat = numpy.array([ant.status     for ant in antennas])
    feeStat = numpy.array([ant.fee.status for ant in antennas])

    badStands = numpy.where( antStat != 3 )[0]
    badFees   = numpy.where( feeStat != 3 )[0]
    bad = numpy.where( (stands > 256) | (antStat != 3) | (feeStat != 3) )[0]
    print "Number of bad stands:   %3i" % len(badStands)
    print "Number of bad FEEs:     %3i" % len(badFees)
    print "---------------------------"
    print "Total number bad inputs: %3i" % len(bad)
    print " "
    
    beams = []
    for p,name in enumerate(('NS', 'EW')):
        antennas2 = []
        for ant in antennas:
            if ant.status == 3 and ant.fee.status == 3 and ant.pol == p:
                antennas2.append(ant)
                

        print "Calculating beam for az. %.2f, el. %.2f at %.2f MHz - %s with %i antennas" % (config['az'], config['el'], config['freq']/1e6, name, len(antennas2))
        tStart = time.time()
        beam = beamformer.intBeamShape(antennas2, azimuth=config['az'], elevation=config['el'], freq=config['freq'], progress=True)
        print "-> Finished in %.3f seconds" % (time.time() - tStart)
        
        # Older version of the intBeamShape function return the square root of the beam profile
        if float(beamformer.__version__) < 0.5:
            beam = beam**2
        beams.append(beam)
        
        numpy.savez('%s_%iMHz_%iaz_%iel_%s.npz' % (station.name, config['freq']/1e6, config['az'], config['el'], name), station=station.name.lower(), beam=beam, freq=config['freq'], pol=name, az=config['az'], el=config['el'])
        
    if config['plots']:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        ax1.imshow( beams[0], interpolation='nearest' )
        ax2.imshow( beams[1], interpolation='nearest' )
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])

