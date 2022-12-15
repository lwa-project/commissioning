#!/usr/bin/env python3

"""
Simple script to read in a MCS binary packed DP gain file (.gf) and print 
out the gains.
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

from lsl.common.dp import dpg_to_gain
from lsl.common import stations


def main(args):

    # Read in the entire file
    fh = open(args.filename, 'rb')
    data = fh.read()
    fh.close()
    
    # Unpack the raw delays (260*4 signed short ints)
    rawGains = struct.unpack('<1040h', data)

    # Convert to delays in ns
    gains = [dpg_to_gain(g) for g in rawGains]

    #Build up the station
    if args.lwasv:
        site = stations.lwasv
    else:
        site = stations.lwa1

    # Report
    ants = site.antennas[0::2]
    
    print("Std   X->x   X->y    Y->x   Y->y")
    print("----------------------------------")
    for i in range(len(ants)):
        xofx, xofy, yofx, yofy = gains[4*i:4*(i+1)]
        print("%3i   %6.4f %6.4f  %6.4f %6.4f" % (ants[i].stand.id, xofx, xofy, yofx, yofy))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Simple script to read in a MCS binary packed DP gain file (.gf) and print out the gains.',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('filename', type=str,
            help='Binary gain file (.gf)')
    parser.add_argument('-v','--lwasv', action='store_true',
            help='Station is LWA-SV')

    args = parser.parse_args()
    main(args)
