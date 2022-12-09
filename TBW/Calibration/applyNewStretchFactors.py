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
from lsl.misc import parser as aph

# List of stands *not* to update
EXCLUDED_STANDS = []


def main(args):
    #
    # Load in the data
    #
    ssmifContents = open(args.ssmif, 'r').readlines()
    site     = parse_ssmif(args.ssmif)
    dataFile = numpy.loadtxt(args.filename)
    
    #
    # Gather the station meta-data from its various sources
    #
    observer = site.get_observer()
    antennas = site.antennas
    
    #
    # Match the new stretch factors to the antennas
    #
    factors = [1.0 for i in xrange(len(antennas))]
    for i in xrange(dataFile.shape[0]):
        dig, stretch, addDelay, rms, chi2 = dataFile[i,:]
        dig = int(dig)
        antenna = antennas[dig-1]
        if antenna.stand.id in args.exclude:
            continue
            
        factors[antenna.id-1] = stretch
        
    #
    # Final results
    #
    if args.output is not None:
        fh = open(args.output, 'w')
    else:
        fh = sys.stdout
        
    for line in ssmifContents:
        if line[0:8] == 'RPD_STR[':
            start = line.find('[')
            stop  = line.find(']')
            try:
                junk, toSave = line.split('#', 1)
                toSave = " # %s" % toSave
            except ValueError:
                toSave = "\n"
            
            antID = int(line[start+1:stop])
            fh.write("RPD_STR[%i]  %.4f%s" % (antID, factors[antID-1], toSave))
        else:
            fh.write(line)
            
    if args.output is not None:
        fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="given an existing SSMIF and new stretch factors, build a new SSMIF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('ssmif', type=str,
                        help='station SSMIF')
    parser.add_argument('filename', type=str, 
                        help='filename to stretch factors')
    parser.add_argument('-e', '--exclude', type=aph.csv_int_list,
                        help='comma seperated list of stands not to update')
    parser.add_argument('-o', '--output', type=str,
                        help='write output to the specified filename instead of the screen')
    args = parser.parse_args()
    main(args)
    