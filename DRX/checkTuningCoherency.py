#!/usr/bin/env python3

"""
Given a DRX file with both tunings set to the same parameters, check for coherency.
"""

import os
import sys
import math
import time
import numpy
import argparse

import lsl.reader.drx as drx
import lsl.reader.errors as errors
from lsl.misc import parser as aph

import matplotlib.pyplot as plt


def crossCorrelate(sig, ref):
    """
    Cross-correlate two signals to get the lag between the two
    in samples.  Returns a two-element tuple of the lag values 
    in samples and the strength of the correlation.
    """
    
    cc = numpy.fft.fft(sig)*numpy.fft.fft(ref).conj()
    cc = numpy.abs(numpy.fft.fftshift(numpy.fft.ifft(cc)))
    lag = numpy.arange(-len(cc)//2,len(cc)//2)

    sigI = sig.real
    refI = ref.real
    ccI = numpy.fft.fft(sigI)*numpy.fft.fft(refI).conj()
    ccI = numpy.abs(numpy.fft.fftshift(numpy.fft.ifft(ccI)))

    sigQ = sig.imag
    refQ = ref.imag
    ccQ = numpy.fft.fft(sigQ)*numpy.fft.fft(refQ).conj()
    ccQ = numpy.abs(numpy.fft.fftshift(numpy.fft.ifft(ccQ)))
    
    return (lag, cc), (lag, ccI), (lag, ccQ)


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
    
    print(junkFrame.header.time_offset)
    beams = drx.get_beam_count(fh)
    tunepols = drx.get_frames_per_obs(fh)
    tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
    beampols = tunepol

    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(args.skip * srate / 4096 * beampols))
    offset = int(1.0 * offset / beampols) * beampols
    args.skip = 1.0 * offset / beampols * 4096 / srate
    fh.seek(offset*drx.FRAME_SIZE)

    # Make sure that the file chunk size contains is an intger multiple
    # of the beampols.
    maxFrames = int(19144/beampols)*beampols

    # Number of frames to integrate over
    toClip = False
    oldAverage = args.plot_range
    if args.plot_range < 4096/srate:		
        toClip = True
        args.plot_range = 4096/srate
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
    print("Offset: %.3f s (%i frames)" % (args.skip, offset))
    print("Plot time: %.3f s (%i frames; %i frames per beam/tune/pol)" % (args.plot_range, nFrames, nFrames // beampols))
    print("Chunks: %i" % nChunks)

    # Sanity check
    if offset > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
    if nFrames > (nFramesFile - offset):
        raise RuntimeError("Requested integration time+offset is greater than file length")

    junkFrame = drx.read_frame(fh)
    b,t,p = junkFrame.id
    while 2*(t-1)+p != 0:
        junkFrame = drx.read_frame(fh)
        b,t,p = junkFrame.id
        print(b,t,p)
    print(fh.tell())
    fh.seek(-drx.FRAME_SIZE, 1)

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
        print("Working on chunk %i, %i frames remaining" % (i, framesRemaining))
        
        count = {}
        data = numpy.zeros((beampols,framesWork*4096//beampols), dtype=numpy.csingle)
        
        # Inner loop that actually reads the frames into the data array
        print("Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0))
        t0 = time.time()
        
        for j in range(framesWork):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = drx.read_frame(fh, verbose=False)
            except errors.EOFError:
                break
            except errors.SyncError:
                #print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/drx.FRAME_SIZE-1))
                continue
                
            beam,tune,pol = cFrame.id
            aStand = 4*(beam-1) + 2*(tune-1) + pol
            #print(aStand, beam, tune, pol)
            if aStand not in standMapper:
                standMapper.append(aStand)
                oStand = 1*aStand
                aStand = standMapper.index(aStand)
                print("Mapping beam %i, tune. %1i, pol. %1i (%2i) to array index %3i" % (beam, tune, pol, oStand, aStand))
            else:
                aStand = standMapper.index(aStand)

            if aStand not in count.keys():
                count[aStand] = 0
            #if cFrame.header.frame_count % 10000 == 0 and args.verbose:
            #	print("%2i,%1i,%1i -> %2i  %5i  %i" % (beam, tune, pol, aStand, cFrame.header.frame_count, cFrame.payload.timetag))

            #print(data.shape, count[aStand]*4096, (count[aStand]+1)*4096, cFrame.payload.data.shape)
            data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.payload.data
            # Update the counters so that we can average properly later on
            count[aStand] += 1
            
        # The plots:  This is setup for the current configuration of 20 beampols
        fig = plt.figure()
        figsX = int(round(math.sqrt(beampols)))
        figsY = beampols // figsX

        t1X = 1
        t1Y = 1

        offset = 0
        samples = 65536
        for sec in range(data.shape[1]//samples):
            if toClip:
                print("Plotting only the first %i samples (%.3f ms) of data" % (samples, oldAverage*1000.0))

            sortedMapper = sorted(standMapper)
            for k, aStand in enumerate(sortedMapper):
                i = standMapper.index(aStand)

                if standMapper[i] % 2 == 0:
                    ref = data[0,:]
                    t1R = t1X
                else:
                    ref = data[1,:]
                    t1R = t1Y
                
                (lag, cc), junkI, junkQ = crossCorrelate(data[i,sec*samples:(sec+1)*samples], 
                                    ref[offset+sec*samples:offset+(sec+1)*samples])
                best = numpy.where( cc == cc.max() )[0][0]
                if args.verbose:
                    print('tune %i pol. %s' % (standMapper[i]%4//2+1, standMapper[i]%2))
                    print(' -> best peak of %.0f at a lag of %i samples' % (cc.max(), lag[best]))
                    print(' -> NCM with tuning 1 of %.3f' % (cc.max()/t1R))
                
                # Plot
                ax = fig.add_subplot(figsX,figsY,k+1)
                ax.plot(lag,  cc,  label='Same', color='blue')

                # Center on the peak
                best = numpy.where( cc == cc.max() )[0][0]
                ax.set_xlim([lag[best-50], lag[best+50]])
            
                ax.set_title('Beam %i, Tune. %i, Pol. %i' % (standMapper[i]//4+1, standMapper[i]%4//2+1, standMapper[i]%2))
                ax.set_xlabel('Lag [samples]')
                ax.set_ylabel('Analysis Sets')

                # Save the tuning 1 values for the peak of the CC function
                if standMapper[i]%4/2+1 == 1:
                    if standMapper[i]%2 == 0:
                        t1X = cc.max()
                    else:
                        t1Y = cc.max()

        plt.show()

        # Save image if requested
        if args.output is not None:
            fig.savefig(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="read in DRX files and check for coherency",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip period at the beginning in seconds')
    parser.add_argument('-p', '--plot-range', type=aph.positive_float, default=2.0, 
                        help='number of seconds of data to show in the I/Q plots')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false', 
                        help='run %(prog)s in silent mode')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for timeseries image')
    args = parser.parse_args()
    main(args)
    
