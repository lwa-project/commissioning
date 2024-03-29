#!/usr/bin/env python3

"""
Check the time tags in a full 520 antenna stand data set.
"""

import os
import sys
import numpy
import argparse
from functools import reduce

from lsl import astro
from lsl.common import stations
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.common.dp import fS


def main(args):
    # Set the station
    if args.lwasv:
        station = stations.lwasv
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    fh = open(args.filename, "rb", buffering=tbn.FRAME_SIZE*10000)

    # Get the first frame and find out what the firt time tag is, which the
    # first frame number is, and what the sample rate it.  From the sample 
    # rate, estimate how the time tag should advance between frames.
    junkFrame = tbn.read_frame(fh)
    sample_rate = tbn.get_sample_rate(fh)
    antpols = len(antennas)
    tagSkip = fS // sample_rate * junkFrame.payload.data.shape[0]
    fh.seek(0)

    # Store the information about the first frame.
    prevTime = junkFrame.payload.timetag
    prevDate = junkFrame.time.datetime
    prevFrame = junkFrame.header.frame_count

    # Report on the file
    print("Filename: %s" % os.path.basename(args.filename))
    print("Date of first frame: %i -> %s" % (prevTime, str(prevDate)))
    print("Sample rate: %i Hz" % sample_rate)
    print("Time tag skip per frame: %i" % tagSkip)

    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols//2+1), pols=[0, 1])
    
    j = 0
    k = 0
    while True:
        try:
            cFrame = tbn.read_frame(fh)
            k += 1
        except errors.EOFError:
            break
        except errors.SyncError:
            #print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FRAME_SIZE-1))
            continue
                
        buffer.append(cFrame)
        cFrames = buffer.get()

        if cFrames is None:
            continue
        
        valid = reduce(lambda x,y: x+int(y.valid), cFrames, 0)
        if valid != antpols:
            print("WARNING: frame count %i at %i missing %.2f%% of frames" % (cFrames[0].header.frame_count, cFrames[0].payload.timetag, float(antpols - valid)/antpols*100))
        
        timetags = numpy.zeros(len(cFrames), dtype=numpy.int64) - 1
        for cFrame in cFrames:
            stand,pol = cFrame.id
            timetags[2*(stand-1)+pol] = cFrame.payload.timetag
            
        if j == 0:
            prevTime  = numpy.median(timetags)
            prevDate  = cFrames[0].time.datetime
            prevFrame = cFrames[0].header.frame_count
            
            j += 1
            continue
        else:
            currTime = numpy.median(timetags)
            currDate  = cFrames[0].time.datetime
            currFrame = cFrames[0].header.frame_count
            
        if currFrame % 1000 == 0:
            print("At frame %i t.t. is %i -> %s" % (currFrame, currTime, currDate))

        if currTime < prevTime:
            print("ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i" % (currTime, currFrame, prevTime, prevFrame))
            print("       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate)))
        if (currTime-prevTime) > tagSkip:
            print("ERROR: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (currTime, currFrame, prevTime, prevFrame))
            print("       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate)))
        for i in range(timetags.size):
            if timetags[i] != currTime:
                print("ERROR: t.t. of dig. %i != frame set median of %i" % (i, currTime))
                print("       -> difference: %i" % (currTime-timetags[i],))
        
        prevTime  = currTime
        prevData  = currDate
        prevFrame = currFrame
        
        j += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a TBN file and check the flow of time', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-v', '--lwasv', action='store_true', 
                        help='use mapping from LWA-SV instead of LWA1')
    args = parser.parse_args()
    main(args)
    
