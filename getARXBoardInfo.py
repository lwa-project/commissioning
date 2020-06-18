#!/usr/bin/env python

"""
Get ARX information about a stand.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import getopt
from datetime import datetime

from lsl.common import stations


def usage(exitCode=None):
    print("""getARXBoardInfo.py - Given a stand, display information about the
ARX board and channel.

Usage: getARXBoardInfo.py [OPTIONS] stand [stand [...]]

Options:
-h, --help                  Display this help information
-m, --metadata              Name of SSMIF file to use for mappings
""")
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['SSMIF'] = ''
    config['args'] = []

    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hm:", ["help", "metadata="])
    except getopt.GetoptError as err:
        # Print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-m', '--metadata'):
            config['SSMIF'] = value
        else:
            assert False
    
    # Add in arguments
    config['args'] = args

    # Return configuration
    return config


def main(args):
    # Parse the command line and get a list of stands to work on
    config = parseOptions(args)
    stands = [int(stand, 10) for stand in config['args']]
    
    # Set the station
    if config['SSMIF'] != '':
        station = stations.parse_ssmif(config['SSMIF'])
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
        c = ant.arx.aspChannel
        print("Stand %i, pol. %i" % (ant.stand.id, ant.pol))
        print("  Antenna: %i" % ant.id)
        print("  ARX Board: %i" % (c/16+1,))
        print("      SN: %s" % ant.arx.id)
        print("      Channel: %s" % ant.arx.channel)
        print("      Control: %s" % ((c+c%2)/2,))


if __name__ == "__main__":
    main(sys.argv[1:])
    
