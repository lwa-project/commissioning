#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a DRX file, look for clip-o-rama (single samples with high instanenous power).

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import math
import time
import numpy
import getopt

from scipy.special import erf

import lsl.reader.drx as drx
import lsl.reader.errors as errors
import lsl.statistics.robust as robust

import matplotlib.pyplot as plt


def usage(exitCode=None):
    print """cliporama.py - Read in DRX files and create a collection of 
timeseries (I/Q) plots.

Usage: cliporama.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-o, --offset                Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-p, --plot-range            Number of seconds of data to show in the I/Q plots
                            (default = 0.0001)
-s, --stats                 Power statistics for the first 0.5 seconds after offset
-n, --no-plot               Do not plot, only identify clip-o-rama events
"""

    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['offset'] = 0.0
    config['average'] = 0.5
    config['maxFrames'] = 19144*4
    config['stats'] = False
    config['doPlot'] = True
    config['args'] = []

    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hqo:sp:n", ["help", "quiet", "offset=", "stats", "plot-range=", "no-plot"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-q', '--quiet'):
            config['verbose'] = False
        elif opt in ('-o', '--offset'):
            config['offset'] = float(value)
        elif opt in ('-p', '--plot-range'):
            config['average'] = float(value)
        elif opt in ('-s', '--stats'):
            config['stats'] = True
        elif opt in ('-n', '--no-plot'):
            config['doPlot'] = False
        else:
            assert False
    
    # Add in arguments
    config['args'] = args

    # Return configuration
    return config


def main(args):
    # Parse command line options
    config = parseOptions(args)
    
    fh = open(config['args'][0], "rb")
    nFramesFile = os.path.getsize(config['args'][0]) / drx.FrameSize
    
    while True:
        junkFrame = drx.readFrame(fh)
        try:
            srate = junkFrame.getSampleRate()
            break
        except ZeroDivisionError:
            pass
    fh.seek(-drx.FrameSize, 1)
    
    beams = drx.getBeamCount(fh)
    tunepols = drx.getFramesPerObs(fh)
    tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
    beampols = tunepol

    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(config['offset'] * srate / 4096 * beampols))
    offset = int(1.0 * offset / beampols) * beampols
    config['offset'] = 1.0 * offset / beampols * 4096 / srate
    fh.seek(offset*drx.FrameSize)

    # Make sure that the file chunk size contains is an intger multiple
    # of the beampols.
    maxFrames = int(config['maxFrames']/beampols)*beampols

    # Setup the statistics data set
    if config['stats']:
        if config['average'] < 0.1:
            config['average'] = 0.5
        
    # Number of frames to integrate over
    nFrames = int(config['average'] * srate / 4096 * beampols)
    nFrames = int(1.0 * nFrames / beampols) * beampols
    config['average'] = 1.0 * nFrames / beampols * 4096 / srate

    # Number of remaining chunks
    nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))

    # File summary
    print "Filename: %s" % config['args'][0]
    print "Beams: %i" % beams
    print "Tune/Pols: %i %i %i %i" % tunepols
    print "Sample Rate: %i Hz" % srate
    print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
    print "---"
    print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
    print "Plot time: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFrames, nFrames / beampols)
    print "Chunks: %i" % nChunks

    # Sanity check
    if offset > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
    if nFrames > (nFramesFile - offset):
        raise RuntimeError("Requested integration time+offset is greater than file length")

    # Align the file handle so that the first frame read in the
    # main analysis loop is from tuning 1, polarization 0
    junkFrame = drx.readFrame(fh)
    b,t,p = junkFrame.parseID()
    while 2*(t-1)+p != 0:
        junkFrame = drx.readFrame(fh)
        b,t,p = junkFrame.parseID()
    fh.seek(-drx.FrameSize, 1)

    # Master loop over all of the file chuncks
    standMapper = []
    for i in range(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFrames - i*maxFrames
        if framesRemaining > maxFrames:
            framesWork = maxFrames
        else:
            framesWork = framesRemaining
        print "Working on chunk %i, %i frames remaining" % (i, framesRemaining)
        
        count = {0:0, 1:0, 2:0, 3:0}
        data = numpy.zeros((beampols,framesWork*4096/beampols), dtype=numpy.csingle)
        
        # Inner loop that actually reads the frames into the data array
        print "Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0)
        t0 = time.time()
        
        for j in xrange(framesWork):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = drx.readFrame(fh, Verbose=False)
            except errors.eofError:
                break
            except errors.syncError:
                #print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/drx.FrameSize-1)
                continue
                
            beam,tune,pol = cFrame.parseID()
            aStand = 2*(tune-1) + pol
            
            data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = numpy.abs(cFrame.data.iq)**2
            
            # Update the counters so that we can average properly later on
            count[aStand] += 1
        
        # Statistics
        print "Running robust statistics"
        means = [robust.mean(data[i,:]) for i in xrange(data.shape[0])]
        stds  = [robust.std(data[i,:])  for i in xrange(data.shape[0])]
        
        if config['stats']:
            ## Report statistics
            print "Mean: %s" % ' '.join(["%.3f" % m for m in means])
            print "StdD: %s" % ' '.join(["%.3f" % s for s in stds ])
            print "Levels:"
            
            ## Count'em up
            j = 0
            counts = [1,]*data.shape[0]
            while (means[i]+j*stds[i] <= 98) and max(counts) != 0:
                counts =[len(numpy.where( numpy.abs(data[i,:] - means[i]) >= j*stds[i] )[0]) for i in xrange(data.shape[0])]
                print " %2isigma (%5.1f%%): %s" % (j, 100.0*(1-erf(j/numpy.sqrt(2))), ' '.join(["%7i (%5.1f%%)" % (c, 100.0*c/data.shape[1]) for c in counts]))
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
                    print " -> Clip-o-rama likely occuring with %i %i-sigma detection on tuning %i, pol %i" % (counts[i], jP, i/2+1, i%2)
                else:
                    print " -> Clip-o-rama likely occuring with %i %i-sigma detections on tuning %i, pol %i" % (counts[i], jP, i/2+1, i%2)
        
        else:
            outfile = os.path.splitext(config['args'][0])[0]
            outfile = '%s.txt' % outfile
            fh = open(outfile, 'w')
            
            # Plot possible clip-o-rama and flag it
            print "Computing power derivatives w.r.t. time"
            deriv = numpy.zeros_like(data)
            for i in xrange(data.shape[0]):
                deriv[i,:] = numpy.roll(data[i,:], -1) - data[i,:]
            
            # The plots:  This is setup for the current configuration of 20 beampols
            print "Plotting"
            fig = plt.figure()
            figsX = int(round(math.sqrt(beampols)))
            figsY = beampols / figsX
            
            for i in xrange(data.shape[0]):
                ax = fig.add_subplot(figsX,figsY,i+1)
                ax.plot(config['offset'] + numpy.arange(0,data.shape[1])/srate, data[i,:])
                
                ## Mark areas of crazy derivatives
                bad = numpy.where( deriv[i,:] > 20*stds[i]*numpy.sqrt(2) )[0]
                for j in bad:
                    fh.write("Clip-o-rama on tuning %i, pol. %i at %.6f seconds\n" % (i/2+1, i%2, config['offset'] + j/srate))
                    print "Clip-o-rama on tuning %i, pol. %i at %.6f seconds" % (i/2+1, i%2, config['offset'] + j/srate)
                    ax.vlines(config['offset'] + j/srate, -10, 100, linestyle='--', color='red', linewidth=2.0)
                
                ## Mark areas of crazy power levels
                bad = numpy.where( data[i,:] == 98 )[0]
                for j in bad:
                    fh.write("Saturation on tuning %i, pol. %i at %.6f seconds\n" % (i/2+1, i%2, config['offset'] + j/srate))
                    print "Saturation on tuning %i, pol. %i at %.6f seconds" % (i/2+1, i%2, config['offset'] + j/srate)
                    ax.vlines(config['offset'] + j/srate, -10, 100, linestyle='-.', color='red')
                
                ax.set_ylim([-10, 100])
                
                ax.set_title('Beam %i, Tune. %i, Pol. %i' % (beam, i/2+1,i%2))
                ax.set_xlabel('Time [seconds]')
                ax.set_ylabel('I$^2$ + Q$^2$')
                
            fh.close()
            if config['doPlot']:
                plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
