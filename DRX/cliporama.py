#!/usr/bin/env python3

"""
Given a DRX file, look for clip-o-rama (single samples with high instanenous power).
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

from scipy.special import erf

import lsl.reader.drx as drx
import lsl.reader.errors as errors
import lsl.statistics.robust as robust
from lsl.misc import parser as aph

import matplotlib.pyplot as plt


def main(args):
    fh = open(args.filename, "rb")
    nFramesFile = os.path.getsize(args.filename) // drx.FRAME_SIZE
    
    while True:
        junkFrame = drx.read_frame(fh)
        try:
            srate = junkFrame.sample_rate
            break
        except ZeroDivisionError:
            pass
    fh.seek(-drx.FRAME_SIZE, 1)
    
    beams = drx.get_beam_count(fh)
    tunepols = drx.get_frames_per_obs(fh)
    tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
    beampols = tunepol

    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(args.offset * srate / 4096 * beampols))
    offset = int(1.0 * offset / beampols) * beampols
    args.offset = 1.0 * offset / beampols * 4096 / srate
    fh.seek(offset*drx.FRAME_SIZE)

    # Make sure that the file chunk size contains is an intger multiple
    # of the beampols.
    maxFrames = int((19144*4)/beampols)*beampols

    # Setup the statistics data set
    if args.stats:
        if args.plot_range < 0.1:
            args.plot_range = 0.5
        
    # Number of frames to integrate over
    nFrames = int(args.plot_range * srate / 4096 * beampols)
    nFrames = int(1.0 * nFrames / beampols) * beampols
    args.plot_range = 1.0 * nFrames / beampols * 4096 / srate

    # Number of remaining chunks
    nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))

    # File summary
    print("Filename: %s" % args.filename)
    print("Beams: %i" % beams)
    print("Tune/Pols: %i %i %i %i" % tunepols)
    print("Sample Rate: %i Hz" % srate)
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate))
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.offset, offset))
    print("Plot time: %.3f s (%i frames; %i frames per beam/tune/pol)" % (args.plot_range, nFrames, nFrames // beampols))
    print("Chunks: %i" % nChunks)

    # Sanity check
    if offset > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
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
    standMapper = []
    for i in xrange(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFrames - i*maxFrames
        if framesRemaining > maxFrames:
            framesWork = maxFrames
        else:
            framesWork = framesRemaining
        print("Working on chunk %i, %i frames remaining" % (i, framesRemaining))
        
        count = {0:0, 1:0, 2:0, 3:0}
        data = numpy.zeros((beampols,framesWork*4096//beampols), dtype=numpy.csingle)
        
        # Inner loop that actually reads the frames into the data array
        print("Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0))
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
            
            data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = numpy.abs(cFrame.payload.data)**2
            
            # Update the counters so that we can average properly later on
            count[aStand] += 1
        
        # Statistics
        print("Running robust statistics")
        means = [robust.mean(data[i,:]) for i in xrange(data.shape[0])]
        stds  = [robust.std(data[i,:])  for i in xrange(data.shape[0])]
        
        if args.stats:
            ## Report statistics
            print("Mean: %s" % ' '.join(["%.3f" % m for m in means]))
            print("StdD: %s" % ' '.join(["%.3f" % s for s in stds ]))
            print("Levels:")
            
            ## Count'em up
            j = 0
            counts = [1,]*data.shape[0]
            while (means[i]+j*stds[i] <= 98) and max(counts) != 0:
                counts =[len(numpy.where( numpy.abs(data[i,:] - means[i]) >= j*stds[i] )[0]) for i in xrange(data.shape[0])]
                print(" %2isigma (%5.1f%%): %s" % (j, 100.0*(1-erf(j/numpy.sqrt(2))), ' '.join(["%7i (%5.1f%%)" % (c, 100.0*c/data.shape[1]) for c in counts])))
                j += 1
            
            ## Why j-2?  Well, j is 1 more than the last iteration.  So, that last iteration 
            ## is j-1,  which is always filled with 0s by construction.  So, the last crazy
            ## bin is j-2.
            jP = j - 2
            if jP > 20:
                counts = [len(numpy.where( numpy.abs(data[i,:] - means[i]) >= jP*stds[i] )[0]) for i in xrange(data.shape[0])]
                for i in xrange(data.shape[0]):
                    if counts[i] > 0:
                        break
                
                if counts[i] == 1:
                    print(" -> Clip-o-rama likely occuring with %i %i-sigma detection on tuning %i, pol %i" % (counts[i], jP, i//2+1, i%2))
                else:
                    print(" -> Clip-o-rama likely occuring with %i %i-sigma detections on tuning %i, pol %i" % (counts[i], jP, i//2+1, i%2))
        
        else:
            outfile = os.path.splitext(args.filename)[0]
            outfile = '%s.txt' % outfile
            fh = open(outfile, 'w')
            
            # Plot possible clip-o-rama and flag it
            print("Computing power derivatives w.r.t. time")
            deriv = numpy.zeros_like(data)
            for i in xrange(data.shape[0]):
                deriv[i,:] = numpy.roll(data[i,:], -1) - data[i,:]
            
            # The plots:  This is setup for the current configuration of 20 beampols
            print("Plotting")
            fig = plt.figure()
            figsX = int(round(math.sqrt(beampols)))
            figsY = beampols // figsX
            
            for i in xrange(data.shape[0]):
                ax = fig.add_subplot(figsX,figsY,i+1)
                ax.plot(args.offset + numpy.arange(0,data.shape[1])/srate, data[i,:])
                
                ## Mark areas of crazy derivatives
                bad = numpy.where( deriv[i,:] > 20*stds[i]*numpy.sqrt(2) )[0]
                for j in bad:
                    fh.write("Clip-o-rama on tuning %i, pol. %i at %.6f seconds\n" % (i//2+1, i%2, args.offset + j/srate))
                    print("Clip-o-rama on tuning %i, pol. %i at %.6f seconds" % (i//2+1, i%2, args.offset + j/srate))
                    ax.vlines(args.offset + j/srate, -10, 100, linestyle='--', color='red', linewidth=2.0)
                
                ## Mark areas of crazy power levels
                bad = numpy.where( data[i,:] == 98 )[0]
                for j in bad:
                    fh.write("Saturation on tuning %i, pol. %i at %.6f seconds\n" % (i//2+1, i%2, args.offset + j/srate))
                    print("Saturation on tuning %i, pol. %i at %.6f seconds" % (i//2+1, i%2, args.offset + j/srate))
                    ax.vlines(args.offset + j/srate, -10, 100, linestyle='-.', color='red')
                
                ax.set_ylim([-10, 100])
                
                ax.set_title('Beam %i, Tune. %i, Pol. %i' % (beam, i//2+1,i%2))
                ax.set_xlabel('Time [seconds]')
                ax.set_ylabel('I$^2$ + Q$^2$')
                
            fh.close()
            if args.do_plot:
                plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="read in a DRX file create a collection of timeseries (I/Q) plots",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-o', '--offset', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip period in seconds between chunks')
    parser.add_argument('-p', '--plot-range', type=aph.positive_float, default=0.5, 
                        help='number of seconds of data to show in the I/Q plots')
    parser.add_argument('-s', '--stats', action='store_true',
                        help='show power statistics for the first 0.5 seconds after outset')
    parser.add_argument('-n', '--no-plot', dest='do_plot', action='store_false',
                        help='do not plot, only identify clip-o-rama events')
    args = parser.parse_args()
    main(args)
    