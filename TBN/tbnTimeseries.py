#!/usr/bin/env python3

"""
Given a TBN file, plot the time series I and Q data as a function of time.
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

from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.correlator import fx as fxc
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.misc import parser as aph

import matplotlib.pyplot as plt


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
    antpols = len(antennas)

    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(args.skip * srate / 512 * antpols))
    offset = int(1.0 * offset / antpols) * antpols
    args.skip = 1.0 * offset / antpols * 512 / srate
    fh.seek(offset*tbn.FRAME_SIZE)

    # Make sure that the file chunk size contains is an intger multiple
    # of the beampols.
    maxFrames = int((2*260*1000)/antpols)*antpols

    # Number of frames to integrate over
    toClip = False
    oldAverage = args.plot_range
    if args.plot_range < maxFrames:		
        toClip = True
        if args.plot_range < 512/srate:
            args.plot_range = 512/srate
    nFrames = int(round(args.plot_range * srate / 512 * antpols))
    print(args.plot_range, nFrames)
    nFrames = int(1.0 * nFrames / antpols) * antpols
    args.plot_range = 1.0 * nFrames / antpols * 512 / srate

    # Number of remaining chunks
    nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))

    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbn.read_frame(fh)
    fh.seek(-tbn.FRAME_SIZE, 1)
    central_freq = junkFrame.central_freq
    beginDate = junkFrame.time.datetime

    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Ant/Pols: %i" % antpols)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz" % central_freq)
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate))
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.skip, offset))
    print("Plot time: %.3f s (%i frames; %i frames per ant)" % (args.plot_range, nFrames, nFrames // antpols))
    print("Chunks: %i" % nChunks)

    # Sanity check
    if offset > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
    if nFrames > (nFramesFile - offset):
        raise RuntimeError("Requested integration time+offset is greater than file length")

    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols//2+1), pols=[0, 1])

    # Master loop over all of the file chunks
    k = 0

    for i in range(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFrames - k
        if framesRemaining > maxFrames:
            framesWork = maxFrames
            data = numpy.zeros((antpols, framesWork*512//antpols), dtype=numpy.csingle)
        else:
            framesWork = framesRemaining + antpols*buffer.nsegments
            data = numpy.zeros((antpols, framesWork//antpols*512), dtype=numpy.csingle)
            framesWork = framesRemaining
            print("Padding from %i to %i frames" % (framesRemaining, framesWork))
        print("Working on chunk %i, %i frames remaining" % (i, framesRemaining))
        
        count = [0 for a in range(len(antennas))]
        
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
                #print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FRAME_SIZE-1))
                continue
                    
            buffer.append(cFrame)
            cFrames = buffer.get()

            if cFrames is None:
                continue
            
            valid = sum(lambda x,y: x+int(y.valid), cFrames, 0)
            if valid != antpols:
                print("WARNING: frame count %i at %i missing %.2f%% of frames" % (cFrames[0].header.frame_count, cFrames[0].payload.timetag, float(antpols - valid)/antpols*100))
                
            for cFrame in cFrames:
                stand,pol = cFrame.header.id
                
                # In the current configuration, stands start at 1 and go up to 260.  So, we
                # can use this little trick to populate the data array
                aStand = 2*(stand-1)+pol
                
                data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
                # Update the counters so that we can average properly later on
                count[aStand] = count[aStand] + 1
            
            j += 1
    
    # Empty the remaining portion of the buffer and integrate what's left
    for cFrames in buffer.flush():
        # Inner loop that actually reads the frames into the data array
        valid = sum(lambda x,y: x+int(y.valid), cFrames, 0)
        if valid != antpols:
            print("WARNING: frame count %i at %i missing %.2f%% of frames" % (cFrames[0].header.frame_count, cFrames[0].payload.timetag, float(antpols - valid)/antpols*100))
        
        for cFrame in cFrames:
            stand,pol = cFrame.header.id
            # In the current configuration, stands start at 1 and go up to 10.  So, we
            # can use this little trick to populate the data array
            aStand = 2*(stand-1)+pol
            
            data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
            # Update the counters so that we can average properly later on
            count[aStand] = count[aStand] + 1
    
    samples = int(round(oldAverage * srate))
    if toClip:
        print("Plotting only the first %i samples (%.3f ms) of data" % (samples, oldAverage*1000.0))

    # Deal with the `keep` options
    if args.keep == 'all':
        antpolsDisp = int(numpy.ceil(antpols/20))
        js = [i for i in range(antpols)]
    else:
        antpolsDisp = int(numpy.ceil(len(args.keep)*2/20))
        if antpolsDisp < 1:
            antpolsDisp = 1
        
        js = []
        for k in args.keep:
            for i,ant in enumerate(antennas):
                if ant.stand.id == k:
                    js.append(i)

    for i in range(antpolsDisp):
        # Normal plotting
        fig = plt.figure()
        figsY = 4
        figsX = 5
        fig.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.20, hspace=0.50)
        for k in range(i*20, i*20+20):
            try:
                j = js[k]
                currTS = data[j,:]
            except IndexError:
                break

            ax = fig.add_subplot(figsX, figsY, (k%20)+1)
            if toClip:
                if args.instantaneous_power:
                    ax.plot(numpy.arange(0,samples)/srate, numpy.abs(currTS[0:samples])**2)
                else:
                    ax.plot(numpy.arange(0,samples)/srate, currTS[0:samples].real, label='Real')
                    ax.plot(numpy.arange(0,samples)/srate, currTS[0:samples].imag, label='Imag')
            else:
                if args.instantaneous_power:
                    ax.plot(numpy.arange(0,data.shape[1])/srate, numpy.abs(currTS)**2)
                else:
                    ax.plot(numpy.arange(0,data.shape[1])/srate, currTS.real, label='Real')
                    ax.plot(numpy.arange(0,data.shape[1])/srate, currTS.imag, label='Imag')
            ax.set_title('Stand: %i (%i); Dig: %i [%i]' % (antennas[j].stand.id, antennas[j].pol, antennas[j].digitizer, antennas[j].combined_status))
            ax.set_xlabel('Time [seconds]')

            if args.instantaneous_power:
                ax.set_ylabel('I$^2$ + Q$^2$')
            else:
                ax.legend(loc=0)
                ax.set_ylim([-127, 127])
                ax.set_ylabel('Output Level')
            
        # Save spectra image if requested
        if args.output is not None:
            base, ext = os.path.splitext(args.output)
            outFigure = "%s-%02i%s" % (base, i+1, ext)
            fig.savefig(outFigure)
            
        plt.draw()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in TBN files and create a collection of timeseries (I/Q) plots', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to plot')
    parser.add_argument('-v', '--lwasv', action='store_true', 
                        help='use mapping from LWA-SV instead of LWA1')
    parser.add_argument('-s', '--skip', type=aph.positive_float, default=0.0, 
                        help='skip period in seconds between chunks')
    parser.add_argument('-p', '--plot-range', type=aph.positive_float, default=0.01, 
                        help='number of seconds of data to show in the I/Q plots')
    parser.add_argument('-k', '--keep', type=aph.csv_int_list, default='all', 
                        help='only display the following comma-seperated list of stands')
    parser.add_argument('-i', '--instantaneous-power', action='store_true', 
                        help='plot I*I + Q*Q instead of the raw samples')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false', 
                        help='run %(prog)s in silent mode')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for timeseries image')
    args = parser.parse_args()
    main(args)
    
