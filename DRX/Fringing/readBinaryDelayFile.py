#!/usr/bin/env python3

"""
Simple script to read in a MCS binary packed DP delay file (.df) and print
out the delays in ns.
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    range = xrange
except NameError:
    pass
    
import os
import sys
import numpy
import struct
import argparse

from lsl.common.dp import dpd_to_delay
from lsl.common import stations


def main(args):

    # Read in the entire file
    fh = open(args.filename, 'rb')
    data = fh.read()
    fh.close()
    
    # Unpack the raw delays (520 unsigned short ints)
    rawDelays = struct.unpack('<520H', data)

    # Convert to delays in ns
    delays = [dpd_to_delay(d) for d in rawDelays]
    
    #Build up the station
    if args.lwasv:
        site = stations.lwasv
    else:
        site = stations.lwa1

    # Report
    ants = site.antennas[0::2]
    
    print("Std   X [ns]    Y [ns]")
    print("----------------------")
    for i in range(len(ants)):
        dx, dy = delays[2*i+0], delays[2*i+1]
        print("%3i   %7.2f  %7.2f" % (ants[i].stand.id, dx, dy))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description='Simple script to read in a MCS binary packed DP delay file (.df) and print out the delays in ns.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('filename', type=str,
            help='Binary delay file (.df)')
    parser.add_argument('-v','--lwasv', action='store_true',
            help='Station is LWA-SV')

    args = parser.parse_args()
    main(args)
