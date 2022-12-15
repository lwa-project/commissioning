#!/usr/bin/env python3

"""
Read in SSMIF file and create a STEPPED/SPEC_DELAYS_GAINS SDF that puts a 
single dipole or beam on the X pol and the outlier on the other.  The 
gain files are constructed such that all data is from X pol.
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    range = xrange
except NameError:
    pass
    
import os
import re
import sys
import numpy
import argparse

from lsl.misc import beamformer
from lsl.common.stations import parse_ssmif
from lsl.common import sdf, sdfADP
from lsl.common.mcs import apply_pointing_correction
from lsl.common.dp import delay_to_dpd, gain_to_dpg
from lsl.misc import parser as aph


# Create the metatag regular expression to deal with spectrometer mode settings
metaRE = re.compile(r'\{.*\}')


def twoByteSwap(i):
    """
    gain_to_dpg and delay_to_dpd return values that are ready for big-
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
    filename = args.filename
    
    # Observation start time
    tStart = "%s %s" % (args.date, args.time)

    station = parse_ssmif(filename)
    #Import the right version of the sdf module for the desired station.
    if station.name == 'LWASV':
        sdf_module = sdfADP
    else:
        sdf_module = sdf

    antennas = station.antennas

    digs    = numpy.array([ant.digitizer  for ant in antennas])
    ants    = numpy.array([ant.id         for ant in antennas])
    stands  = numpy.array([ant.stand.id   for ant in antennas])
    pols    = numpy.array([ant.pol        for ant in antennas])
    antStat = numpy.array([ant.status     for ant in antennas])
    feeStat = numpy.array([ant.fee.status for ant in antennas])

    badStands = numpy.where( antStat != 3 )[0]
    badFees   = numpy.where( feeStat != 3 )[0]
    bad = numpy.where( (stands > 256) | (antStat != 3) | (feeStat != 3) )[0]
    print("Number of bad stands:   %3i" % len(badStands))
    print("Number of bad FEEs:     %3i" % len(badFees))
    print("---------------------------")
    print("Total number bad inputs: %3i" % len(bad))
    print(" ")
    
    # Adjust the gain so that it matches the outlier better
    if args.dipole == 0:
        bgain = 20.0 / (len(antennas) - len(bad))
    else:
        bgain = 1.0
        
    # Setup the base gain lists for the different scenarios
    baseEmptyGain = [0.0000, 0.0000, 0.0000, 0.0000]
    if not args.y_pol:
        baseBeamGain    = [bgain,  0.0000, 0.0000, 0.0000]
        baseDipoleGain  = [0.0000, 1.0000, 0.0000, 0.0000]
    else:
        baseBeamGain    = [0.0000, 0.0000, bgain,  0.0000]
        baseDipoleGain  = [0.0000, 0.0000, 0.0000, 1.0000]
        
    if (args.dipole == 0):
        freq = max([args.frequency1, args.frequency2])
        
        # Load in the pointing correction
        pcTerms = getPointingTerms(filename)
        print("Applying Pointing Correction Terms: theta=%.2f, phi=%.2f, psi=%.2f" % pcTerms)
        az, el = apply_pointing_correction(args.azimuth, args.elevation, *pcTerms)
        print("-> az %.2f, el %.2f to az %.2f, el %.2f" % (args.azimuth, args.elevation, az, el))
        print(" ")
        
        print("Calculating delays for az. %.2f, el. %.2f at %.2f MHz" % (az, el, freq/1e6))
        delays = beamformer.calc_delay(antennas, freq=freq, azimuth=az, elevation=el)
        delays *= 1e9
        delays = delays.max() - delays
        delays = [twoByteSwap(delay_to_dpd(d)) for d in delays]
        
        print("Setting gains for %i good inputs, %i bad inputs" % (len(antennas)-len(bad), len(bad)))
        print("-> Using gain setting of %.4f for the beam" % bgain)
        
        gains = [[twoByteSwap(gain_to_dpg(g)) for g in baseBeamGain] for i in xrange(int(len(antennas)/2))] # initialize gain list 
        
        for d in digs[bad]:
            # Digitizers start at 1, list indicies at 0
            i = d - 1
            gains[i//2] = [twoByteSwap(gain_to_dpg(g)) for g in baseEmptyGain]
            
        for i in range(len(stands)//2):
            # Put the reference stand in there all by itself
            if stands[2*i] == args.reference:
                gains[i] = [twoByteSwap(gain_to_dpg(g)) for g in baseDipoleGain]
    else:
        print("Setting all delays to zero")
        delays = [0 for i in antennas]
        delays = [twoByteSwap(delay_to_dpd(d)) for d in delays]
        
        print("Setting gains for dipoles %i and %i" % (args.dipole, args.reference))
        
        gains = [[twoByteSwap(gain_to_dpg(g)) for g in baseEmptyGain] for i in range(int(len(antennas)/2))] # initialize gain list
        for i in range(len(stands)//2):
            # Put the fringing stand in there all by itself
            if stands[2*i] == args.dipole:
                gains[i] = [twoByteSwap(gain_to_dpg(g)) for g in baseBeamGain]
            
            # Put the reference stand in there all by itself
            if stands[2*i] == args.reference:
                gains[i] = [twoByteSwap(gain_to_dpg(g)) for g in baseDipoleGain]
    
    # Resort the gains into a list of 2x2 matrices
    newGains = []
    for gain in gains:
        newGains.append([[gain[0], gain[1]], [gain[2], gain[3]]])
    gains = newGains
    
    # Create the SDF
    sessionComment = 'Input Pol.: %s; Output Pol.: beam -> X, reference -> Y' % ('X' if not args.y_pol else 'Y',)
    observer = sdf_module.Observer("fringeSDF.py Observer", 99)
    session = sdf_module.Session("fringeSDF.py Session", 1, comments=sessionComment)
    project = sdf_module.Project(observer, "fringeSDF.py Project", "FRINGSDF", [session,])
    obs = sdf_module.Stepped("fringeSDF.py Target", "Custom", tStart, args.filter, is_radec=False)
    stp = sdf_module.BeamStep(args.azimuth, args.elevation, str(args.obs_length), args.frequency1, args.frequency2, is_radec=False, spec_delays=delays, spec_gains=gains)
    obs.append(stp)
    obs.gain = 1
    project.sessions[0].observations.append(obs)
    project.sessions[0].drx_beam = args.dp_beam
    ## Spectrometer setup
    if args.spec_setup is not None:
        # Remove the ' marks
        args.spec_setup = args.spec_setup.replace("'", "")
        # Excise the metatags
        mtch = metaRE.search(args.spec_setup)
        if mtch is not None:
            metatag = mtch.group(0)
            args.spec_setup = metaRE.sub('', args.spec_setup)
        else:
            metatag = None
            
        project.sessions[0].spcSetup = [int(i) for i in args.spec_setup.lstrip().rstrip().split(None, 1)]
        project.sessions[0].spcMetatag = metatag
        
    # Write it out
    if os.path.exists(args.output):
        raise RuntimeError("File '%s' already exists" % args.output)
    project.render(verbose=True)
    fh = open(args.output, 'w')
    fh.write(project.render())
    fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="read in SSMIF file and create a STEPPED/SPEC_DELAYS_GAINS SDF that puts a single dipole or beam on the X pol and the outlier on the other",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str,
                        help='SSMIF to use')
    parser.add_argument('date', type=aph.date,
                        help='UTC date for the start of the observation in YYYY/MM/DD')
    parser.add_argument('time', type=aph.time,
                        help='UTC time for the start the observation in HH:MM:SS')
    parser.add_argument('-a', '--azimuth', type=float, default=90.0,
                        help='beam only, azimuth east of north in degrees for the pointing center')
    parser.add_argument('-e', '--elevation', type=aph.positive_or_zero_float, default=90.0,
                        help='beam only, elevation above the horizon in degrees for the pointing center')
    parser.add_argument('-d', '--dipole', type=aph.positive_or_zero_int, default=0,
                        help='use the specified dipole instead of the beam; 0 = use beam')
    parser.add_argument('-y', '--y-pol', action='store_true', 
                        help='generate an SDF for the Y polarization instead of X')
    parser.add_argument('-r', '--reference', type=aph.positive_int, default=258,
                        help='reference dipole for the fringing')
    parser.add_argument('-b', '--dp-beam', type=aph.positive_int, default=2,
                        help='DP beam to run the observation on')
    parser.add_argument('-l', '--obs-length', type=aph.positive_float, default=3600.0,
                        help='duration of the observation in seconds')
    parser.add_argument('-1', '--frequency1', type=aph.positive_float, default=37.9,
                        help='frequency in MHz for Tuning #1')
    parser.add_argument('-2', '--frequency2', type=aph.positive_float, default=74.0,
                        help='frequency in MHz for Tuning #2')
    parser.add_argument('-f', '--filter', type=aph.positive_int, default=7,
                        help='DRX filter code')
    parser.add_argument('-s', '--spec-setup', type=str,
                        help='DR spectrometer setup to use, i.e., "32 6144{Stokes=IV}"') 
    parser.add_argument('-o', '--output', type=str, default='fringe.sdf',
                        help='filename to save the SDF to')
    args = parser.parse_args()
    args.frequency1 *= 1e6
    args.frequency2 *= 1e6
    main(args)
    
