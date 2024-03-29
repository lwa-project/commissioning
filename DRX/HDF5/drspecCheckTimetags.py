#!/usr/bin/env python3

"""
Check the time times in a DR spectrometer file for flow.
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
        for i in range(5):
            junkFrame = drx.read_frame(fh)
        raise RuntimeError("ERROR: '%s' appears to be a raw DRX file, not a DR spectrometer file" % args.filename)
    except errors.SyncError:
        fh.seek(0)
        
    # Interrogate the file to figure out what frames sizes to expect, now many 
    # frames there are, and what the transform length is
    FRAME_SIZE = drspec.get_frame_size(fh)
    nFrames = os.path.getsize(args.filename) / FRAME_SIZE
    nChunks = nFrames
    LFFT = drspec.get_transform_size(fh)

    # Read in the first frame to figure out the DP information
    junkFrame = drspec.read_frame(fh)
    fh.seek(-FRAME_SIZE, 1)
    srate = junkFrame.sample_rate
    t0 = junkFrame.time
    tInt = junkFrame.header.nints*LFFT/srate
    
    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(skip / tInt))
    fh.seek(offset*FRAME_SIZE, 1)
    
    # Iterate on the offsets until we reach the right point in the file.  This
    # is needed to deal with files that start with only one tuning and/or a 
    # different sample rate.  
    while True:
        ## Figure out where in the file we are and what the current tuning/sample 
        ## rate is
        junkFrame = drspec.read_frame(fh)
        srate = junkFrame.sample_rate
        t1 = junkFrame.time
        tInt = junkFrame.header.nints*LFFT/srate
        fh.seek(-FRAME_SIZE, 1)
        
        ## See how far off the current frame is from the target
        tDiff = t1 - (t0 + skip)
        
        ## Half that to come up with a new seek parameter
        tCorr = -tDiff / 2.0
        cOffset = int(round(tCorr / tInt))
        offset += cOffset
        
        ## If the offset is zero, we are done.  Otherwise, apply the offset
        ## and check the location in the file again/
        if cOffset == 0:
            break
        fh.seek(cOffset*FRAME_SIZE, 1)
        
    # Update the offset actually used
    skip = t1 - t0
    nChunks = (os.path.getsize(args.filename) - fh.tell()) // FRAME_SIZE
    
    # Update the file contents
    beam = junkFrame.id
    central_freq1, central_freq2 = junkFrame.central_freq
    srate = junkFrame.sample_rate
    data_products = junkFrame.data_products
    t0 = junkFrame.time
    tInt = junkFrame.header.nints*LFFT/srate
    beginDate = junkFrame.time.datetime
        
    # Report
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % beginDate)
    print("Beam: %i" % beam)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (central_freq1, central_freq2))
    print("Data Products: %s" % ','.join(data_products))
    print("Frames: %i (%.3f s)" % (nFrames, nFrames*tInt))
    print("---")
    print("Transform Length: %i" % LFFT)
    print("Integration: %.3f s" % tInt)
    
    for i in range(nChunks):
        frame = drspec.read_frame(fh)
        
        cTime = frame.time
        if i % 1000 == 0:
            print("Frame %i: %s" % (i, cTime.datetime))
            
        try:
            if cTime > oTime + 1.001*tInt:
                print('Warning: Time tag error at frame %i; %.3f > %.3f + %.3f' % (i, cTime, oTime, tInt))
        except NameError:
            pass
        oTime = frame.time
        
        cFreq1, cFreq2 = frame.central_freq
        try:
            if cFreq1 != oFreq1:
                print('Warning: Tuning 1 frequncy changed at frame %i; %.3f Hz != %.3f Hz' % (i, cFreq1, oFreq1))
            if cFreq2 != oFreq2:
                print('Warning: Tuning 2 frequncy changed at frame %i; %.3f Hz != %.3f Hz' % (i, cFreq2, oFreq2))
        except NameError:
            pass
        oFreq1, oFreq2 = frame.central_freq
        
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
    
