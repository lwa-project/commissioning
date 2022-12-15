#!/usr/bin/env python3

"""
Read in SSMIF file and create a set of DRX gain and delay files for a given 
frequency and topogentric pointing center.
"""

# Python2 compatibility
from __future__ import print_function, division

import sys
import numpy
import argparse

import gain
import delay
from lsl.misc import beamformer
from lsl.common.stations import parse_ssmif
from lsl.misc import parser as aph


def main(args):
    filename = args.filename

    station = parse_ssmif(filename)
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
    print("Total number bad inuts: %3i" % len(bad))
    print(" ")

    dftBase = 'beams_%iMHz_%iaz_%iel_%03ibg' % (args.frequency/1e6, args.azimuth, args.elevation, args.gain*100)
    gftBase = 'beams_%iMHz_%iaz_%iel_%03ibg' % (args.frequency/1e6, args.azimuth, args.elevation, args.gain*100)

    print("Calculating delays for az. %.2f, el. %.2f at %.2f MHz" % (args.azimuth, args.elevation, args.frequency/1e6))
    delays = beamformer.calc_delay(antennas, freq=args.frequency, azimuth=args.azimuth, elevation=args.elevation)
    delays *= 1e9
    delays = delays.max() - delays
    junk = delay.list2delayfile('.', dftBase, delays)

    print("Setting gains for %i good inputs, %i bad inputs" % (len(antennas)-len(bad), len(bad)))
    bgain = args.gain
    bgain_cross = 0.0000
    gains = [[bgain, bgain_cross, bgain_cross, bgain]]*260 # initialize gain list
    for d in digs[bad]:
        # Digitizers start at 1, list indicies at 0
        i = d - 1
        gains[i/2] = [0,0,0,0]
    junk = gain.list2gainfile('.', gftBase, gains)

    print("\nDelay and gain files are:\n %s.dft\n %s.gft" % (dftBase, gftBase))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="read in a SSMIF file and create a set of DRX gain and delay files for a given frequency and topogentric pointing center",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename of SSMIF to use')
    parser.add_argument('-f', '--frequency', type=aph.positive_float, default=65.0,
                        help='frequency in MHz to calculate the gain/delays for')
    parser.add_argument('-a', '--azimuth', type=float, default=90.0,
                        help='azimuth east of north in degrees for the pointing center')
    parser.add_argument('-e', '--elevation', type=aph.positive_or_zero_float, default=90.0,
                        help='elevation above the horizon in degrees for the pointing')
    parser.add_argument('-g', '--gain', type=float, default=1.0,
                        help='DRX antenna gain')
    args = parser.parse_args()
    args.frequency *= 1e6
    main(args)
    
