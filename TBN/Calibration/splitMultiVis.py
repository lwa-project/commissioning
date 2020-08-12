#!/usr/bin/env python

"""
Split multi-vis NPZ files that are generated from combined TBN observations into
single NPZ files, one for each frequency.

Usage:
./splitMultiVis.py <NPZ multi-vis. file> [<NPZ multi-vis. file> [...]]
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import numpy
import argparse

from datetime import datetime

from lsl.statistics import robust
from lsl.misc import parser as aph


def main(args):
    filenames = args.filename

    for filename in filenames:
        dataDict = numpy.load(filename)
    
        print("Working on file '%s'" % filename)

        # Load in the data
        ref_ant = dataDict['ref'].item()
        refX   = dataDict['refX'].item()
        refY   = dataDict['refY'].item()
        tInt = dataDict['tInt'].item()
    
        times = dataDict['times']
        phase = dataDict['simpleVis']
    
        central_freqs = dataDict['central_freqs']

        # Load in the SSMIF
        ssmifContents = dataDict['ssmifContents']

        # Find the unique sets of (non-zero) frequencies and report
        uFreq = numpy.unique(central_freqs)
        uFreq = uFreq[numpy.where(uFreq != 0)]
        print("  Found %i unique frequencies from %.3f to %.3f MHz" % (len(uFreq), uFreq.min()/1e6, uFreq.max()/1e6))
    
        # Report on the start time
        beginDate = datetime.utcfromtimestamp(times[0])
        print("  Start date/time of data: %s UTC" % beginDate)
        
        # Gather
        tInts = []
        for i,f in enumerate(uFreq):
            ## Select what we need and trim off the last index to deal 
            ## with frequency changes
            toKeep = numpy.where( central_freqs == f )[0]
            toKeep = toKeep[:-1]
            
            tInts.append(len(toKeep))
        cutLength = numpy.median(tInts) - 2
        if args.auto_integrations:
            print("  Setting minimum integration count to %i (%.1f s)" % (cutLength, cutLength*tInt))
            args.min_integrations = cutLength
    
        # Split
        for i,f in enumerate(uFreq):
            ## Select what we need and trim off the last index to deal 
            ## with frequency changes
            toKeep = numpy.where( central_freqs == f )[0]
            toKeep = toKeep[:-1]
        
            ## Sub-sets of `times` and `phase`
            subTimes = times[toKeep]
            subPhase = phase[toKeep,:]
        
            if len(toKeep) < args.min_integrations:
                print("  -> Skipping %.3f MHz with only %i integrations (%.1f s)" % (f, len(toKeep), len(toKeep)*tInt))
                continue
        
            ## Save the split data to its own file
            outname = os.path.splitext(filename)[0]
            outname = outname.replace('-multi', '')
            outname = "%s-%03i.npz" % (outname, i+1)
            print("  Saving visibility data for %.3f MHz to '%s'" % (f/1e6, outname))
            numpy.savez(outname, ref=ref_ant, refX=refX, refY=refY, tInt=tInt, central_freq=f, 
                        times=subTimes, simpleVis=subPhase, ssmifContents=ssmifContents)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="split multi-vis NPZ files that are generated from combined TBN observations into single NPZ files, one for each frequency",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+', 
                        help='file to spllit')
    parser.add_argument('-a', '--auto-integrations', action='store_true',
                        help='use automatic integration selection (median count, less 2 integrations)')
    parser.add_argument('-m', '--min-integrations', type=aph.positive_int, default=20,
                        help='minimum number of integrations needed to keep a split out frequency')
    args = parser.parse_args()
    main(args)
    