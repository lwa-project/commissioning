#!/usr/bin/env python3

"""
Get ARX information about a stand.
"""

import os
import sys
import argparse
from datetime import datetime

from lsl.common import stations
from lsl.misc import parser as aph


def main(args):
    # Get a list of stands to work on
    stands = args.stand
    
    # Set the station
    if args.metadata is not None:
        station = stations.parse_ssmif(args.metadata)
    else:
        station = stations.lwa1
        
    # Match the stands to ASP channels
    ants = []
    antennas = station.antennas
    for stand in stands:
        for antenna in antennas:
            if antenna.stand.id == stand:
                ants.append( antenna )
                
    # Report
    for ant in ants:
        c = ant.arx.asp_channel
        print("Stand %i, pol. %i" % (ant.stand.id, ant.pol))
        print("  Antenna: %i" % ant.id)
        print("  ARX Board: %i" % (c/16+1,))
        print("      SN: %s" % ant.arx.id)
        print("      Channel: %s" % ant.arx.channel)
        print("      Control: %s" % ((c+c%2)/2,))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='given a stand, display information about the ARX board and channel',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('stand', type=aph.positive_int, nargs='+',
                        help='stand ID to query')
    parser.add_argument('-m', '--metadata', type=str,
                        help='name of SSMIF file to use for mappings')
    args = parser.parse_args()
    main(args)
    
