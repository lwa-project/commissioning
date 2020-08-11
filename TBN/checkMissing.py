#!/usr/bin/env python

"""
Given a TBN file, check for missing frames (or frames considered missing by the
ring buffer) and plot what the missing packet rate and what packets might be 
missing.  Rather than do this for a whole file, it is done for some small portion
of the file that is controlled by the -s/--skip and -a/--average flags.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import math
import numpy
import argparse

from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.correlator import fx as fxc
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.misc import parser as aph

import matplotlib.pyplot as plt


def plotMissing(ax1, ax2, missingPackets, missingList, antpols):
    d = ax1.imshow(missingPackets, vmin=0, vmax=1, cmap=plt.cm.gray)
    ax1.set_xlabel('Missing Packets Set')
    ax1.set_ylabel('Digitizer-1')
    #fig.colorbar(d, ax=ax1)
                    
    ax2.plot(numpy.array(missingList)/float(antpols))
    ax2.set_xlabel('Frame Count Relative to File Start')
    ax2.set_ylabel('Fraction Missing Packets')


def main(args):
    # Set the station
    if args.lwasv:
        station = stations.lwasv
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    fh = open(args.filename, "rb")
    nFramesFile = os.path.getsize(args.filename) // tbn.FRAME_SIZE
    srate = tbn.get_sample_rate(fh)
    #antpols = tbn.get_frames_per_obs(fh)
    antpols = len(antennas)

    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(args.skip * srate / 512 * antpols)
    offset = int(1.0 * offset / antpols) * antpols
    args.skip = 1.0 * offset / antpols * 512 / srate
    fh.seek(offset*tbn.FRAME_SIZE)

    # Number of frames to integrate over
    nFrames = int(args.average * srate / 512 * antpols)
    args.average = 1.0 * nFrames / antpols * 512 / srate

    # Number of remaining chunks
    nChunks = int(math.ceil(1.0*(nFrames)/(200*520)))

    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbn.read_frame(fh)
    fh.seek(0)
    beginDate = junkFrame.time.datetime

    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Ant/Pols: %i" % antpols)
    print("Sample Rate: %i Hz" % srate)
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate))
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.skip, offset))
    print("Integration: %.3f s (%i frames; %i frames per stand/pol)" % (args.average, nFrames, nFrames / antpols))
    print("Chunks: %i" % nChunks)

    # Sanity check
    if offset > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
    if nFrames > (nFramesFile - offset):
        raise RuntimeError("Requested integration time+offset is greater than file length")

    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols//2+1), pols=[0, 1])

    # Master loop over all of the file chunks
    masterCount = [0 for a in xrange(len(antennas))]
    
    # Missing packet control variables
    missingPackets = numpy.ones((antpols, 2048), dtype=numpy.int8)
    pc = 0
    missing = 0
    missingList = []
    
    # Figure
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    
    k = 0
    for i in xrange(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time ((200*520)),
        # only deal with that chunk
        framesRemaining = nFrames - k
        if framesRemaining > (200*520):
            framesWork = (200*520)
        else:
            framesWork = framesRemaining
        
        count = [0 for a in xrange(len(antennas))]
        
        j = 0
        fillsWork = framesWork // antpols
        # Inner loop that actually reads the frames into the data array
        while j < fillsWork:
            try:
                cFrame = tbn.read_frame(fh)
                k = k + 1
            except errors.EOFError:
                break
            except errors.SyncError:
                print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FRAME_SIZE-1))
                continue
            
            #print(cFrame.header.frame_count, cFrame.payload.timetag, cFrame.id)
            
            buffer.append(cFrame)
            cFrames = buffer.get()

            if cFrames is None:
                continue
            
            valid = sum(lambda x,y: x+int(y.valid), cFrames, 0)
            print("Frame #%5i:  %.4f seconds with %i valid ant/pols%s" % (cFrames[0].header.frame_count, cFrames[0].time, valid, '!' if valid != antpols else ''))
            if valid != antpols:
                bad = []
                for cFrame in cFrames:
                    if not cFrame.valid:
                        bad.append(cFrame.id)
                        
                        missingPackets[2*(bad[-1][0]-1)+bad[-1][1], pc] = 0
                bad.sort()
                
                pc += 1
                if pc == missingPackets.shape[1]:
                    plotMissing(ax1, ax2, missingPackets, missingList, antpols)
                    plt.show()
                    sys.exit(0)

                missing += (antpols-valid)
                total = (buffer.full + buffer.partial)*antpols
                #print(j, valid, antpols-valid, cFrames[0].header.frame_count, 1.0*missing / total* 100, bad[0], bad[-1], buffer.dropped)
                #print(buffer.status())
                
                missingList.append( antpols - valid )
            else:
                total = (buffer.full + buffer.partial)*antpols
                missingList.append(0)
                
            times = numpy.array([f.payload.timetag for f in cFrames], dtype=numpy.int64)
            #print(cFrames[0].header.frame_count, times.min(), times.max(), times.max()-times.min(), "%6.3f%%" % (1.0*missing/total*100,))
            for cFrame in cFrames:
                stand,pol = cFrame.header.id
                
                # In the current configuration, stands start at 1 and go up to 260.  So, we
                # can use this little trick to populate the data array
                aStand = 2*(stand-1)+pol
                
                # Update the counters so that we can average properly later on
                count[aStand] = count[aStand] + 1
                masterCount[aStand] = masterCount[aStand] + 1
            
            j += 1
    
    # Empty the remaining portion of the buffer and integrate what's left
    for cFrames in buffer.flush():
        valid = sum(lambda x,y: x+int(y.valid), cFrames, 0)
        print("Frame #%5i:  %.4f seconds with %i valid ant/pols" % (cFrames[0].header.frame_count, cFrames[0].time, valid))
        if valid != antpols:
            bad = []
            for cFrame in cFrames:
                if not cFrame.valid:
                    bad.append(cFrame.id)
                    
                    missingPackets[2*(bad[-1][0]-1)+bad[-1][1], pc] = 0
            bad.sort()
            
            pc += 1
            if pc == missingPackets.shape[1]:
                plotMissing(ax1, ax2, missingPackets, missingList, antpols)
                plt.show()
                sys.exit(0)

            missing += (antpols-valid)
            total = (buffer.full + buffer.partial)*antpols
            #print(j, valid, antpols-valid, cFrames[0].header.frame_count, 1.0*missing / total* 100, bad[0], bad[-1], buffer.dropped)
            #print(buffer.status())
            
            missingList.append( antpols - valid )
        else:
            total = (buffer.full + buffer.partial)*antpols
            missingList.append(0)
        
        # Inner loop that actually reads the frames into the data array
        for cFrame in cFrames:
            stand,pol = cFrame.header.id
            # In the current configuration, stands start at 1 and go up to 10.  So, we
            # can use this little trick to populate the data array
            aStand = 2*(stand-1)+pol
                            
            # Update the counters so that we can average properly later on
            count[aStand] = count[aStand] + 1
            masterCount[aStand] = masterCount[aStand] + 1

        j += 1

    plotMissing(ax1, ax2, missingPackets, missingList, antpols)
    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a TBN file and find out what is missing or appears to be missing', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-v', '--lwasv', action='store_true', 
                        help='use mapping from LWA-SV instead of LWA1')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-a', '--average', type=aph.positive_float, default=10.0, 
                        help='number of seconds of data to examine for frame loss')
    args = parser.parse_args()
    main(args)
    
