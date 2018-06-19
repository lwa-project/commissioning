#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read in SSMIF file and create a STEPPED/SPEC_DELAYS_GAINS SDF that puts a 
single dipole or beam on the X pol and the outlier on the other.  The 
gain files are constructed such that all data is from X pol.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import re
import sys
import numpy
import getopt

from lsl.misc import beamformer
from lsl.common.stations import parseSSMIF
from lsl.common import sdf
from lsl.common.mcs import applyPointingCorrection
from lsl.common.dp import delaytoDPD, gaintoDPG

def usage(exitCode=None):
    print """fringeSets.py - Read in SSMIF file and create a STEPPED/SPEC_DELAYS_GAINS SDF
that puts a single dipole or beam on the X pol and the outlier on the other.

Usage: fringeSDF.py [OPTIONS] SSMIF YYYY/MM/DD HH:MM:SS.SSS

Options:
-h, --help                  Display this help information
-a, --azimuth               Beam only, azimuth east of north in degrees for the 
                            pointing center (Default = 90 degrees)
-e, --elevation             Beam only, elevation above the horizon in degrees for 
                            the pointing center (Default = 90 degrees)
-d, --dipole                Using a dipole instead of the beam (Default = use beam)
-y, --y-pol                 Generate an SDF for the Y polarization (Default = X)
-r, --reference             Reference for the fringing (Default = stand #258)
-b, --drx-beam              DP beam to run the observation on (Default = 2)
-l, --obs-length            Duration of the observation in seconds (Default = 3600.0)
-1, --frequency1            Frequency in MHz for Tuning #1 (Default = 37.9 MHz)
-2, --frequency2            Frequency in MHz for Tuning #2 (Default = 74.0 MHz)
-f, --filter                DRX filter code (Default = 7)
-s, --spec-setup            Spectrometer setup to use, i.e., "32 6144{Stokes=IV}" 
                            (Default = do not use DR spectrometer)
-o, --output                Filename to save the SDF to (Default = fringe.sdf)
"""
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['beam'] = True
    config['az'] = 90.0
    config['el'] = 90.0
    config['dipole'] = 1
    config['xPol'] = True
    config['ref'] = 258
    config['drxBeam'] = 2
    config['duration'] = 3600.0
    config['freq1'] = 37.9e6
    config['freq2'] = 74.0e6
    config['filter'] = 7
    config['spcSetup'] = [0, 0]
    config['spcMetatag'] = ""
    config['output'] = 'fringe.sdf'
    config['args'] = []
    
    # Create the metatag regular expression to deal with spectrometer mode settings
    metaRE = re.compile(r'\{.*\}')
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "ha:e:d:yr:b:l:1:2:f:s:o:", ["help", "azimuth=", "elevation=", "dipole=", "y-pol", "reference=", "drx-beam=", "obs-length=", "frequency1=", "frequency2=", "filter=", "spec-setup=", "output="])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-a', '--azimuth'):
            config['az'] = float(value)
        elif opt in ('-e', '--elevation'):
            config['el'] = float(value)
        elif opt in ('-d', '--dipole'):
            config['beam'] = False
            config['dipole'] = int(value)
        elif opt in ('-y', '--y-pol'):
            config['xPol'] = False
        elif opt in ('-r', '--reference'):
            config['ref'] = int(value)
        elif opt in ('-b', '--dp-beam'):
            config['drxBeam'] = int(value)
        elif opt in ('-l', '--obs-length'):
            config['duration'] = float(value)
        elif opt in ('-1', '--frequency1'):
            config['freq1'] = float(value)*1e6
        elif opt in ('-2', '--frequency2'):
            config['freq2'] = float(value)*1e6
        elif opt in ('-f', '--filter'):
            config['filter'] = int(value)
        elif opt in ('-s', '--spec-setup'):
            # Remove the ' marks
            value = value.replace("'", "")
            # Excise the metatags
            mtch = metaRE.search(value)
            if mtch is not None:
                metatag = mtch.group(0)
                value = metaRE.sub('', value)
            else:
                metatag = None
            
            config['spcSetup'] = [int(i) for i in value.lstrip().rstrip().split(None, 1)]
            config['spcMetatag'] = metatag
        elif opt in ('-o', '--output'):
            config['output'] = value
        else:
            assert False
    
    # Add in arguments
    config['args'] = args

    # Return configuration
    return config


def twoByteSwap(i):
    """
    gaintoDPG and delaytoDPD return values that are ready for big-
    endian packing, MCS is expecting little-endian.
    """

    return ((i & 0xFF) << 8) | ((i >> 8) & 0xFF)


def getPointingTerms(filename):
    """
    Return a three-element tuple of the pointing correction terms (theta, 
    phi, psi) stored in the SSMIF.
    """
    
    theta = 0.0
    phi = 0.0
    psi = 0.0
    
    fh = open(filename)
    for line in fh:
        line = line.replace('\n', '')
        if len(line) < 3:
            continue
        
        fields = line.split()
        if fields[0] == 'PC_AXIS_TH':
            theta = float(fields[1])
        elif fields[0] == 'PC_AXIS_PH':
            phi = float(fields[1])
        elif fields[0] == 'PC_ROT':
            psi = float(fields[1])
        else:
            pass
    fh.close()
    
    return (theta, phi, psi)


def main(args):
    config = parseOptions(args)
    filename = config['args'][0]
    
    # Observation start time
    config['args'][1] = config['args'][1].replace('-', '/')
    tStart = "%s %s" % (config['args'][1], config['args'][2])

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
    
    # Adjust the gain so that it matches the outlier better
    if config['beam']:
        bgain = 20.0 / (520 - len(bad))
    else:
        bgain = 1.0
        
    # Setup the base gain lists for the different scenarios
    baseEmptyGain = [0.0000, 0.0000, 0.0000, 0.0000]
    if config['xPol']:
        baseBeamGain    = [bgain,  0.0000, 0.0000, 0.0000]
        baseDipoleGain  = [0.0000, 1.0000, 0.0000, 0.0000]
    else:
        baseBeamGain    = [0.0000, 0.0000, bgain,  0.0000]
        baseDipoleGain  = [0.0000, 0.0000, 0.0000, 1.0000]
        
    if config['beam']:
        freq = max([config['freq1'], config['freq2']])
        
        # Load in the pointing correction
        pcTerms = getPointingTerms(filename)
        print "Applying Pointing Correction Terms: theta=%.2f, phi=%.2f, psi=%.2f" % pcTerms
        az, el = applyPointingCorrection(config['az'], config['el'], *pcTerms)
        print "-> az %.2f, el %.2f to az %.2f, el %.2f" % (config['az'], config['el'], az, el)
        print " "
        
        print "Calculating delays for az. %.2f, el. %.2f at %.2f MHz" % (az, el, freq/1e6)
        delays = beamformer.calcDelay(antennas, freq=freq, azimuth=az, elevation=el)
        delays *= 1e9
        delays = delays.max() - delays
        delays = [twoByteSwap(delaytoDPD(d)) for d in delays]
        
        print "Setting gains for %i good inputs, %i bad inputs" % (len(antennas)-len(bad), len(bad))
        print "-> Using gain setting of %.4f for the beam" % bgain
        
        gains = [[twoByteSwap(gaintoDPG(g)) for g in baseBeamGain] for i in xrange(260)] # initialize gain list
        for d in digs[bad]:
            # Digitizers start at 1, list indicies at 0
            i = d - 1
            gains[i/2] = [twoByteSwap(gaintoDPG(g)) for g in baseEmptyGain]
            
        for i in xrange(len(stands)/2):
            # Put the reference stand in there all by itself
            if stands[2*i] == config['ref']:
                gains[i] = [twoByteSwap(gaintoDPG(g)) for g in baseDipoleGain]
    else:
        print "Setting all delays to zero"
        delays = [0 for i in antennas]
        delays = [twoByteSwap(delaytoDPD(d)) for d in delays]
        
        print "Setting gains for dipoles %i and %i" % (config['dipole'], config['ref'])
        
        gains = [[twoByteSwap(gaintoDPG(g)) for g in baseEmptyGain] for i in xrange(260)] # initialize gain list
        for i in xrange(len(stands)/2):
            # Put the fringing stand in there all by itself
            if stands[2*i] == config['dipole']:
                gains[i] = [twoByteSwap(gaintoDPG(g)) for g in baseBeamGain]
            
            # Put the reference stand in there all by itself
            if stands[2*i] == config['ref']:
                gains[i] = [twoByteSwap(gaintoDPG(g)) for g in baseDipoleGain]
    
    # Resort the gains into a list of 260 2x2 matrices
    newGains = []
    for gain in gains:
        newGains.append([[gain[0], gain[1]], [gain[2], gain[3]]])
    gains = newGains
    
    # Create the SDF
    sessionComment = 'Input Pol.: %s; Output Pol.: beam -> X, reference -> Y' % ('X' if config['xPol'] else 'Y',)
    observer = sdf.Observer("fringeSDF.py Observer", 99)
    session = sdf.Session("fringeSDF.py Session", 1, comments=sessionComment)
    project = sdf.Project(observer, "fringeSDF.py Project", "FRINGSDF", [session,])
    obs = sdf.Stepped("fringeSDF.py Target", "Custom", tStart, config['filter'], RADec=False)
    stp = sdf.BeamStep(config['az'], config['el'], str(config['duration']), config['freq1'], config['freq2'], RADec=False, SpecDelays=delays, SpecGains=gains)
    obs.append(stp)
    obs.gain = 1
    project.sessions[0].observations.append(obs)
    project.sessions[0].drxBeam = config['drxBeam']
    project.sessions[0].spcSetup = config['spcSetup']
    project.sessions[0].spcMetatag = config['spcMetatag']
    
    # Write it out
    if os.path.exists(config['output']):
        raise RuntimeError("File '%s' already exists" % config['output'])
    project.render(verbose=True)
    fh = open(config['output'], 'w')
    fh.write(project.render())
    fh.close()


if __name__ == "__main__":
    main(sys.argv[1:])
