#!/usr/bin/env python

"""
Given a TBW file, check the time tags.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import math
import ephem
import numpy
import argparse

from lsl.common import stations
from lsl.reader import tbw, tbn
from lsl.reader import errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt

def main(args):
    # Set the station
    station = stations.lwa1
    antennas = station.antennas

    fh = open(args.filename, "rb")
    nFrames = os.path.getsize(args.filename) // tbw.FRAME_SIZE
    dataBits = tbw.get_data_bits(fh)
    # The number of ant/pols in the file is hard coded because I cannot figure out 
    # a way to get this number in a systematic fashion
    maxFrames = 30000*260
    antpols = len(antennas)
    nChunks = int(math.ceil(1.0*nFrames/maxFrames))
    if dataBits == 12:
        nSamples = 400
    else:
        nSamples = 1200

    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbw.read_frame(fh)
    fh.seek(0)
    beginDate = ephem.Date(unix_to_utcjd(junkFrame.get_time()) - DJD_OFFSET)

    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Ant/Pols: %i" % antpols)
    print("Sample Length: %i-bit" % dataBits)
    print("Frames: %i" % nFrames)
    print("Chunks: %i" % nChunks)
    print("===")

    nChunks = 1

    # Skip over any non-TBW frames at the beginning of the file
    i = 0
    junkFrame = tbw.read_frame(fh)
    while not junkFrame.header.is_tbw:
        try:
            junkFrame = tbw.read_frame(fh)
        except errors.SyncError:
            fh.seek(0)
            while True:
                try:
                    junkFrame = tbn.read_frame(fh)
                    i += 1
                except errors.SyncError:
                    break
            fh.seek(-2*tbn.FRAME_SIZE, 1)
            junkFrame = tbw.read_frame(fh)
        i += 1
    fh.seek(-tbw.FRAME_SIZE, 1)
    print("Skipped %i non-TBW frames at the beginning of the file" % i)

    # Master loop over all of the file chunks
    timetags = numpy.zeros((antpols, 30000), dtype=numpy.int64) - 1
    for i in range(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFrames - i*maxFrames
        if framesRemaining > maxFrames:
            framesWork = maxFrames
        else:
            framesWork = nFrames
        print("Working on chunk %i, %i frames remaining" % ((i+1), framesRemaining))

        # Inner loop that actually reads the frames into the data array
        for j in range(framesWork):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = tbw.read_frame(fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbw.FRAME_SIZE-1))
                continue
            if not cFrame.header.is_tbw:
                continue
            
            stand = cFrame.header.id
            # In the current configuration, stands start at 1 and go up to 10.  So, we
            # can use this little trick to populate the data array
            aStand = 2*(stand-1)
            if cFrame.header.frame_count % 10000 == 0:
                print("%3i -> %3i  %5i  %i" % (stand, aStand, cFrame.header.frame_count, cFrame.data.timetag))

            # Actually load the data.  x pol goes into the even numbers, y pol into the 
            # odd numbers
            count = cFrame.header.frame_count - 1
            timetags[aStand,   count] = cFrame.data.timetag
            timetags[aStand+1, count] = cFrame.data.timetag

    # Check for missing frames
    missing = numpy.where( timetags < 0 )
    if len(missing) != 0:
        dp1Boards = {}
        print("Found %i missing frames (%i missing time tags).  Missing data from:" % (len(missing[0])/2, len(missing[0])))
        for i,f in zip(missing[0], missing[1]):
            try:
                dp1Boards[antennas[i].board] += 1
            except KeyError:
                dp1Boards[antennas[i].board] = 1

            print("  stand %3i, pol. %1i (dig. %3i) @ frame %5i" % (antennas[i].stand.id, antennas[i].pol, antennas[i].digitizer, f+1))
        print("-> DP1 boards with missing frames:")
        for k in dp1Boards.keys():
            v = dp1Boards[k]
            print("   %2i %6i (%7.3f%%)" % (k, v, 100.0*v/(30000*10)))

    # Check time tags to make sure every ant/pol as the same time as each frame
    for f in xrange(timetags.shape[1]):
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
    for i in xrange(timetags.shape[0]):
        for f in xrange(1,timetags.shape[1]):
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
        description='read in a TBW file and check the flow of time', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    args = parser.parse_args()
    main(args)
    