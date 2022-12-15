#!/usr/bin/env python3

"""
Given a TBF file, check the time tags.
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    range = xrange
except NameError:
    pass
        
import os
import sys
import math
import numpy
import argparse

from lsl.common import stations
from lsl.reader import tbf
from lsl.reader import errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def main(args):
    fh = open(args.filename, "rb")
    nFrames = os.path.getsize(args.filename) / tbf.FRAME_SIZE
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbf.read_frame(fh)
    fh.seek(0)
    beginDate = junkFrame.time.datetime
    
    # Figure out how many frames there are per observation and the number of
    # channels that are in the file
    nFramesPerObs = tbf.get_frames_per_obs(fh)
    nchannels = tbf.get_channel_count(fh)
    nSamples = 7840
    
    # Figure out how many chunks we need to work with
    nChunks = nFrames / nFramesPerObs
    
    # Pre-load the channel mapper
    mapper = []
    for i in range(2*nFramesPerObs):
        cFrame = tbf.read_frame(fh)
        if cFrame.header.first_chan not in mapper:
            mapper.append( cFrame.header.first_chan )
    fh.seek(-2*nFramesPerObs*tbf.FRAME_SIZE, 1)
    mapper.sort()
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Frames per Observation: %i" % nFramesPerObs)
    print("Channel Count: %i" % nchannels)
    print("Frames: %i" % nFrames)
    print("===")
    print("Chunks: %i" % nChunks)
    
    # Master loop over all of the file chunks
    timetags = numpy.zeros((nFramesPerObs, nChunks), dtype=numpy.int64) - 1
    for i in range(nChunks):
        # Inner loop that actually reads the frames into the data array
        for j in range(nFramesPerObs):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = tbf.read_frame(fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbf.FRAME_SIZE-1))
                continue
            if not cFrame.header.is_tbf:
                continue
                
            first_chan = cFrame.header.first_chan
            
            # Figure out where to map the channel sequence to
            try:
                aStand = mapper.index(first_chan)
            except ValueError:
                mapper.append(first_chan)
                aStand = mapper.index(first_chan)
            
            if cFrame.header.frame_count % 10000 == 0:
                print("%4i -> %4i  %7i  %i" % (first_chan, aStand, cFrame.header.frame_count, cFrame.payload.timetag))
                
            # Actually load the data.  x pol goes into the even numbers, y pol into the 
            # odd numbers
            if i == 0 and j == 0:
                refCount = cFrame.header.frame_count
            count = cFrame.header.frame_count - refCount
            timetags[aStand,   count] = cFrame.payload.timetag
            
    # Check for missing frames
    missing = numpy.where( timetags < 0 )
    if len(missing) != 0:
        print("Found %i missing frames.  Missing data from:" % len(missing[0]))
        for i,f in zip(missing[0], missing[1]):
            print("  channel set %4i @ frame %5i" % (mapper[i], f+1))
            
    # Check time tags to make sure every ant/pol as the same time as each frame
    for f in range(timetags.shape[1]):
        ## For each frame count value, get the median time tag and use this for comparison.
        ## If things are really bad, we will get a lot of errors.
        frameTime = numpy.median( timetags[:,f] )

        ## Compare all of the antpols at a particular frame count, ignoring the ones that
        ## are missing.
        missing = numpy.where( (timetags[:,f] != frameTime) & (timetags[:,f]>=0) )[0]

        ## Report any errors
        for m in missing:
            print("ERROR: t.t. %i @ frame %i != frame median of %i" % (timetags[m,f], f+1, frameTime))
            print("       -> difference: %i" % (timetags[m,f]-frameTime,))

    # Check time tags to make sure the times increment correctly between frames
    for i in range(timetags.shape[0]):
        for f in range(1,timetags.shape[1]):
            ## Skip missing frames since they always fail
            if timetags[i,f] < 0 or timetags[i,f-1] < 0:
                continue

            ## Compare the current time tag with previous and report an error if there
            ## is a discrepancy between the two modulo the expected skip.
            if timetags[i,f] > (timetags[i,f-1] + nSamples):
                ## Too far into the future
                print("ERROR: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (timetags[i,f], f+1, timetags[i,f-1], f))
                print("       -> difference: %i" % (timetags[i,f]-timetags[i,f-1],))
            elif timetags[i,f] < (timetags[i,f-1] + nSamples):
                ## Not far enough into the future
                print("ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i + skip" % (timetags[i,f], f+1, timetags[i,f-1], f))
                print("       -> difference: %i" % (timetags[i,f]-timetags[i,f-1],))
            else:
                ## Everything is good if we make it here
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='given a TBF file, check the time tags', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    args = parser.parse_args()
    main(args)
    
    
