#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Check the time times in a DR spectrometer file for flow.

$Rev: 1658 $
$LastChangedBy: jayce $
$LastChangedDate: 2014-05-12 09:25:49 -0600 (Mon, 12 May 2014) $
"""

import os
import sys
import numpy
import gc
import argparse
from datetime import datetime

from lsl.reader import drx, drspec, errors
from lsl.misc import parser as aph


def main(args):
    skip = args.skip
    fh = open(args.filename, "rb")
        
    try:
        for i in xrange(5):
            junkFrame = drx.readFrame(fh)
        raise RuntimeError("ERROR: '%s' appears to be a raw DRX file, not a DR spectrometer file" % args.filename)
    except errors.syncError:
        fh.seek(0)
        
    # Interrogate the file to figure out what frames sizes to expect, now many 
    # frames there are, and what the transform length is
    FrameSize = drspec.getFrameSize(fh)
    nFrames = os.path.getsize(args.filename) / FrameSize
    nChunks = nFrames
    LFFT = drspec.getTransformSize(fh)

    # Read in the first frame to figure out the DP information
    junkFrame = drspec.readFrame(fh)
    fh.seek(-FrameSize, 1)
    srate = junkFrame.getSampleRate()
    t0 = junkFrame.getTime()
    tInt = junkFrame.header.nInts*LFFT/srate
    
    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(skip / tInt))
    fh.seek(offset*FrameSize, 1)
    
    # Iterate on the offsets until we reach the right point in the file.  This
    # is needed to deal with files that start with only one tuning and/or a 
    # different sample rate.  
    while True:
        ## Figure out where in the file we are and what the current tuning/sample 
        ## rate is
        junkFrame = drspec.readFrame(fh)
        srate = junkFrame.getSampleRate()
        t1 = junkFrame.getTime()
        tInt = junkFrame.header.nInts*LFFT/srate
        fh.seek(-FrameSize, 1)
        
        ## See how far off the current frame is from the target
        tDiff = t1 - (t0 + skip)
        
        ## Half that to come up with a new seek parameter
        tCorr = -tDiff / 2.0
        cOffset = int(round(tCorr / tInt))
        offset += cOffset
        
        ## If the offset is zero, we are done.  Otherwise, apply the offset
        ## and check the location in the file again/
        if cOffset is 0:
            break
        fh.seek(cOffset*FrameSize, 1)
        
    # Update the offset actually used
    skip = t1 - t0
    nChunks = (os.path.getsize(args.filename) - fh.tell()) / FrameSize
    
    # Update the file contents
    beam = junkFrame.parseID()
    centralFreq1 = junkFrame.getCentralFreq(1)
    centralFreq2 = junkFrame.getCentralFreq(2)
    srate = junkFrame.getSampleRate()
    dataProducts = junkFrame.getDataProducts()
    t0 = junkFrame.getTime()
    tInt = junkFrame.header.nInts*LFFT/srate
    beginDate = datetime.utcfromtimestamp(junkFrame.getTime())
        
    # Report
    print "Filename: %s" % args.filename
    print "Date of First Frame: %s" % beginDate
    print "Beam: %i" % beam
    print "Sample Rate: %i Hz" % srate
    print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
    print "Data Products: %s" % ','.join(dataProducts)
    print "Frames: %i (%.3f s)" % (nFrames, nFrames*tInt)
    print "---"
    print "Transform Length: %i" % LFFT
    print "Integration: %.3f s" % tInt
    
    for i in xrange(nChunks):
        frame = drspec.readFrame(fh)
        
        cTime = frame.getTime()
        if i % 1000 == 0:
            print "Frame %i: %s" % (i, datetime.utcfromtimestamp(cTime))		
            
        try:
            if cTime > oTime + 1.001*tInt:
                print 'Warning: Time tag error at frame %i; %.3f > %.3f + %.3f' % (i, cTime, oTime, tInt)
        except NameError:
            pass
        oTime = frame.getTime()
        
        cFreq1 = frame.getCentralFreq(1)
        cFreq2 = frame.getCentralFreq(2)
        try:
            if cFreq1 != oFreq1:
                print 'Warning: Tuning 1 frequncy changed at frame %i; %.3f Hz != %.3f Hz' % (i, cFreq1, oFreq1)
            if cFreq2 != oFreq2:
                print 'Warning: Tuning 2 frequncy changed at frame %i; %.3f Hz != %.3f Hz' % (i, cFreq2, oFreq2)
        except NameError:
            pass
        oFreq1 = frame.getCentralFreq(1)
        oFreq2 = frame.getCentralFreq(2)
        
        del frame
        if i % 100 == 0:
            gc.collect()

        
    # Done
    fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a DR spectrometer file and check the flow of time', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip period in seconds between chunks')
    args = parser.parse_args()
    main(args)
    