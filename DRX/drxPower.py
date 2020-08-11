#!/usr/bin/env python

"""
Given a DRX file, plot the instantaneous power as a function of time.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import math
import time
import numpy
import argparse

from lsl import astro
import lsl.reader.drx as drx
import lsl.reader.errors as errors
from lsl.misc import parser as aph

import matplotlib.pyplot as plt


def main(args):
    fh = open(args.filename, "rb")
    nFramesFile = os.path.getsize(args.filename) // drx.FRAME_SIZE
    
    while True:
        junkFrame = drx.read_frame(fh)
        try:
            srate = junkFrame.sample_rate
            t0i, t0f = junkFrame.time
            break
        except ZeroDivisionError:
            pass
    fh.seek(-drx.FRAME_SIZE, 1)
    
    beams = drx.get_beam_count(fh)
    tunepols = drx.get_frames_per_obs(fh)
    tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
    beampols = tunepol

    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(args.skip * srate / 4096 * beampols))
    offset = int(1.0 * offset / beampols) * beampols
    fh.seek(offset*drx.FRAME_SIZE, 1)
    
    # Iterate on the offsets until we reach the right point in the file.  This
    # is needed to deal with files that start with only one tuning and/or a 
    # different sample rate.  
    while True:
        ## Figure out where in the file we are and what the current tuning/sample 
        ## rate is
        junkFrame = drx.read_frame(fh)
        srate = junkFrame.sample_rate
        t1i, t1f = junkFrame.time
        tunepols = drx.get_frames_per_obs(fh)
        tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
        beampols = tunepol
        fh.seek(-drx.FRAME_SIZE, 1)
        
        ## See how far off the current frame is from the target
        tDiff = t1i - (t0i + args.skip) + (t1f - t0f)
        
        ## Half that to come up with a new seek parameter
        tCorr = -tDiff / 2.0
        cOffset = int(tCorr * srate / 4096 * beampols)
        cOffset = int(1.0 * cOffset / beampols) * beampols
        offset += cOffset
        
        ## If the offset is zero, we are done.  Otherwise, apply the offset
        ## and check the location in the file again/
        if cOffset is 0:
            break
        fh.seek(cOffset*drx.FRAME_SIZE, 1)
    
    # Update the offset actually used
    args.skip = t1i - t0i + t1f - t0f
    offset = int(round(args.skip * srate / 4096 * beampols))
    offset = int(1.0 * offset / beampols) * beampols

    # Make sure that the file chunk size contains is an intger multiple
    # of the beampols.
    maxFrames = int(round(args.average*srate/4096))*beampols
    if maxFrames < beampols:
        maxFrames = beampols
    args.average = 1.0*maxFrames/beampols*4096/srate

    # Number of remaining chunks
    nChunks = int(round(args.duration / args.average))
    nFrames = maxFrames * nChunks

    # Store the information about the first frame.
    prevTime = junkFrame.payload.timetag
    prevDate = junkFrame.time.datetime

    # File summary
    print("Filename: %s" % args.filename)
    print("Beams: %i" % beams)
    print("Tune/Pols: %i %i %i %i" % tunepols)
    print("Sample Rate: %i Hz" % srate)
    print("Date of first frame: %i -> %s" % (prevTime, str(prevDate)))
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate))
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.skip, offset))
    print("Integration: %.4f s (%i frames; %i frames per beam/tune/pol)" % (args.average, maxFrames, maxFrames // beampols))
    print("Duration: %.4f s (%i frames; %i frames per beam/tune/pol)" % (args.average*nChunks, nFrames, nFrames // beampols))
    print(" ")

    # Sanity check
    if nFrames > (nFramesFile - offset):
        raise RuntimeError("Requested integration time+offset is greater than file length")

    # Align the file handle so that the first frame read in the
    # main analysis loop is from tuning 1, polarization 0
    junkFrame = drx.read_frame(fh)
    b,t,p = junkFrame.id
    while 2*(t-1)+p != 0:
        junkFrame = drx.read_frame(fh)
        b,t,p = junkFrame.id
    fh.seek(-drx.FRAME_SIZE, 1)

    # Master loop over all of the file chuncks
    masterTimes = numpy.zeros((nChunks, 4), dtype=numpy.float64)
    masterData  = numpy.zeros((nChunks, 4), dtype=numpy.float32)
    masterData2 = numpy.zeros((nChunks, 4), dtype=numpy.float32)
    masterProcessed = [0, 0, 0, 0]
    masterRemoved = [0, 0, 0, 0]
    for i in range(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFrames - i*maxFrames
        if framesRemaining > maxFrames:
            framesWork = maxFrames
        else:
            framesWork = framesRemaining
        #print("Working on chunk %i, %i frames remaining" % (i, framesRemaining))
        
        count = {0:0, 1:0, 2:0, 3:0}
        data  = numpy.zeros((4,framesWork*4096//beampols), dtype=numpy.float32)
        data2 = numpy.zeros((4,framesWork*4096//beampols), dtype=numpy.float32)
        
        ## Inner loop that actually reads the frames into the data array
        #print("Working on %.2f ms of data" % ((framesWork*4096/beampols/srate)*1000.0))
        t0 = time.time()
        
        for j in xrange(framesWork):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = drx.read_frame(fh, verbose=False)
            except errors.EOFError:
                break
            except errors.SyncError:
                #print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/drx.FRAME_SIZE-1))
                continue
                
            beam,tune,pol = cFrame.id
            aStand = 2*(tune-1) + pol
            
            if j < 4:
                masterTimes[i,aStand] = cFrame.time

            try:
                framePower = numpy.abs(cFrame.payload.data)**2
                
                # Calculate the clipping
                mask = numpy.where( framePower <= args.trim_level, 1, 0 )
                
                data[ aStand,  count[aStand]*4096:(count[aStand]+1)*4096] = framePower
                data2[aStand,  count[aStand]*4096:(count[aStand]+1)*4096] = framePower*mask
                
                masterProcessed[aStand] += mask.size
                masterRemoved[aStand] += (mask.size - mask.sum())
            
                # Update the counters so that we can average properly later on
                count[aStand] += 1
            except ValueError:
                pass

        # Save the data
        masterData[i,:]  = data.sum(axis=1)
        masterData2[i,:] = data2.sum(axis=1)
        
    # Really save the data to a NPZ file
    if args.write_npz:
        outfile = os.path.split(args.filename)[1]
        outfile = os.path.splitext(outfile)[0]
        outfile = '%s-power.npz' % outfile
        numpy.savez(outfile, beam=beam, avgPowerFull=masterData, avgPowerTrim=masterData2, times=masterTimes, trimLevel=args.trim_level)

    # Report on the clipping
    print("Summary:")
    for i in xrange(4):
        print("  Tuning %i, Pol. %s:" % (i/2+1, 'X' if i%2 else 'Y'))
        print("    Processed: %i samples" % masterProcessed[i])
        print("    Clipped:   %i samples" % masterRemoved[i])
        print("      -> %.1f%% blanked" % (100.0*masterRemoved[i]/masterProcessed[i],))

    # The plots:  This is setup for the current configuration of 20 beampols
    fig = plt.figure()
    figsX = int(round(math.sqrt(4)))
    figsY = 4 // figsX

    for i in xrange(masterData.shape[1]):
        ax = fig.add_subplot(figsX,figsY,i+1)
        ax.plot(numpy.arange(0, masterData.shape[0] )*args.average, masterData[:,i],  label='Full')
        ax.plot(numpy.arange(0, masterData2.shape[0])*args.average, masterData2[:,i], label='Trimmed')
        ax.set_ylim([0, masterData.max()])
        
        ax.set_title('Beam %i, Tune. %i, Pol. %i' % (beam, i//2+1, i%2))
        ax.set_xlabel('Time [seconds]')
        ax.set_ylabel('Output Power Level')

        ax.legend(loc=0)

    # Part 2, polarization stuff
    fig2 = plt.figure()
    ax = fig2.add_subplot(3, 2, 1)
    ax.plot(numpy.arange(0, masterData.shape[0])*args.average, numpy.sqrt(masterData[:,0]**2 + masterData[:,1]**2))
    ax.set_title('$\\sqrt{X1^2 + Y1^2}$')
    ax.set_xlabel('Time [seconds]')

    ax = fig2.add_subplot(3, 2, 2)
    ax.plot(numpy.arange(0, masterData.shape[0])*args.average, masterData[:,1] / masterData[:,0])
    ax.set_title('$Y1 / X1$')
    ax.set_xlabel('Time [seconds]')

    ax = fig2.add_subplot(3, 2, 3)
    ax.plot(numpy.arange(0, masterData.shape[0])*args.average, numpy.sqrt(masterData[:,2]**2 + masterData[:,3]**2))
    ax.set_title('$\\sqrt{X2^2 + Y2^2}$')
    ax.set_xlabel('Time [seconds]')

    ax = fig2.add_subplot(3, 2, 4)
    ax.plot(numpy.arange(0, masterData.shape[0])*args.average, masterData[:,3] / masterData[:,2])
    ax.set_title('$Y2 / X2$')
    ax.set_xlabel('Time [seconds]')

    ax = fig2.add_subplot(3, 2, 5)
    ax.plot(numpy.arange(0, masterData.shape[0])*args.average, numpy.sqrt(masterData[:,2]**2 + masterData[:,3]**2) / numpy.sqrt(masterData[:,0]**2 + masterData[:,1]**2))
    ax.set_title('$\\sqrt{X2^2 + Y2^2} / \\sqrt{X1^2 + Y1^2}$')
    ax.set_xlabel('Time [seconds]')

    ax = fig2.add_subplot(3, 2, 6)
    ax.plot(numpy.arange(0, masterData.shape[0])*args.average, (masterData[:,3]/masterData[:,2]) / (masterData[:,1]/masterData[:,0]))
    ax.set_title('$(Y2 / X2) / (Y1 / X1)$')
    ax.set_xlabel('Time [seconds]')

    plt.show()

    # Save image if requested
    if args.output is not None:
        fig.savefig(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="read in DRX files and create a collection of timeseries power (I*I + Q*Q) plots",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to analyze')
    parser.add_argument('-t', '--trim-level', type=int, default=49,
                        help="trim level for power analysis with clipping")
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip period at the beginning in seconds')
    parser.add_argument('-a', '--average', type=aph.positive_float, default=0.0002,
                        help='number of seconds of data to average together for power')
    parser.add_argument('-d', '--duration', type=aph.positive_float, default=10.0,
                        help='number of seconds to analyze')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false',
                        help='run %(prog)s in silent mode')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for averaged power image')
    parser.add_argument('-w', '--write-npz', action='store_true',
                        help='rite a NPZ file of the averaged power')
    args = parser.parse_args()
    main(args)
    