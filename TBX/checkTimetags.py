#!/usr/bin/env python3

"""
Given a TBX file, check the time tags.
"""

import os
import sys
import math
import numpy as np
import argparse

from lsl.common import stations
from lsl.reader import tbx, errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def main(args):
    fh = open(args.filename, "rb")
    tbx.FRAME_SIZE = tbx.get_frame_size(fh)
    nFrames = os.path.getsize(args.filename) // tbx.FRAME_SIZE
    
    # Read in the first frame and get the number of stands and date/time of
    # the first sample of the frame.  This is needed to get the list of stands.
    junkFrame = tbx.read_frame(fh)
    nAnt = junkFrame.nstand
    fh.seek(0)
    beginDate = junkFrame.time.datetime
    
    # Figure out how many frames there are per observation and the number of
    # channels that are in the file
    nFramesPerObs = tbx.get_frames_per_obs(fh)
    nchannels = tbx.get_channel_count(fh)
    nSamples = 8192
    
    # Figure out how many chunks we need to work with
    nChunks = nFrames // nFramesPerObs
    
    # Pre-load the channel mapper
    mapper = []
    nread = 0
    while len(mapper) < nFramesPerObs:
        cFrame = tbx.read_frame(fh)
        if cFrame.header.first_chan not in mapper:
            mapper.append( cFrame.header.first_chan )
        nread += 1
    fh.seek(-nread*tbx.FRAME_SIZE, 1)
    mapper.sort()
    
    # File summary
    print(f"Filename: {args.filename}")
    print(f"Date of First Frame: {str(beginDate)}")
    print(f"Frames per Observation: {nFramesPerObs}")
    print(f"Channel Count: {nchannels}")
    print(f"Frames: {nFrames}")
    print("===")
    print(f"Chunks: {nChunks}")
    
    # Master loop over all of the file chunks
    timetags = np.zeros((nFramesPerObs, nChunks), dtype=np.int64) - 1
    for i in range(nChunks):
        # Inner loop that actually reads the frames into the data array
        for j in range(nFramesPerObs):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = tbx.read_frame(fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                print(f"WARNING: Mark 5C sync error on frame #{fh.tell()//tbx.FRAME_SIZE-1}")
                continue
            if not cFrame.header.is_tbx:
                continue
                
            first_chan = cFrame.header.first_chan
            
            # Figure out where to map the channel sequence to
            try:
                aStand = mapper.index(first_chan)
            except ValueError:
                mapper.append(first_chan)
                aStand = mapper.index(first_chan)
            
            if cFrame.header.frame_count % 10000 == 0:
                print(f"{first_chan:4d} -> {aStand:4d}  {cFrame.header.frame_count:7d}  {cFrame.payload.timetag}")
                
            # Actually load the data.  x pol goes into the even numbers, y pol into the 
            # odd numbers
            if i == 0 and j == 0:
                refCount = cFrame.header.frame_count
            count = cFrame.header.frame_count - refCount        # pylint: disable=possibly-used-before-assignment
            timetags[aStand,   count] = cFrame.payload.timetag
            
    # Check for missing frames
    missing = np.where( timetags < 0 )
    if len(missing[0]) != 0:
        print(f"Found {len(missing[0])} missing frames.  Missing data from:")
        for i,f in zip(missing[0], missing[1]):
            print(f"  channel set {mapper[i]:4d} @ frame {f+1:5d}")
            
    # Check time tags to make sure every ant/pol as the same time as each frame
    for f in range(timetags.shape[1]):
        ## For each frame count value, get the median time tag and use this for comparison.
        ## If things are really bad, we will get a lot of errors.
        frameTime = np.median( timetags[:,f] )

        ## Compare all of the antpols at a particular frame count, ignoring the ones that
        ## are missing.
        missing = np.where( (timetags[:,f] != frameTime) & (timetags[:,f]>=0) )[0]

        ## Report any errors
        for m in missing:
            print(f"ERROR: t.t. {timetags[m,f]} @ frame {f+1} != frame median of {frameTime}")
            print(f"       -> difference: {timetags[m,f]-frameTime}")

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
                print(f"ERROR: t.t. {timetags[i,f]} @ frame {f+1} > t.t. {timetags[i,f-1]} @ frame {f} + skip")
                print(f"       -> difference: {timetags[i,f]-timetags[i,f-1]}")
            elif timetags[i,f] < (timetags[i,f-1] + nSamples):
                ## Not far enough into the future
                print(f"ERROR: t.t. {timetags[i,f]} @ frame {f+1} < t.t. {timetags[i,f-1]} @ frame {f} + skip")
                print(f"       -> difference: {timetags[i,f]-timetags[i,f-1]}")
            else:
                ## Everything is good if we make it here
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='given a TBX file, check the time tags', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    args = parser.parse_args()
    main(args)
