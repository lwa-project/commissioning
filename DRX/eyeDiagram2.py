#!/usr/bin/env python3

"""
Given a DRX file, look for glitches in a DRX or TBN data by fitting a sine wave 
to the data.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import numpy
import argparse

from lsl.reader import errors, tbn, drx
from lsl.misc import parser as aph

from matplotlib import pyplot as plt


def quantizedSine(x, scale=1.0):
    return numpy.round(scale*numpy.round(2**11*numpy.sin(x)))


def qsDiff(B, x, y):
    scale = B[0]
    freq = B[1]
    phase = B[2]
    
    xPrime = 2*numpy.pi*freq*x + phase
    #yFit = quantizedSine(xPrime, scale=scale)
    yFit = numpy.sin(xPrime)*numpy.abs(scale)
    
    #return yFit - y
    return ((yFit - y)**2).sum()


def main(args):
    # Filename
    filename = args.filename
    
    # The observation and sinusoid frequencies as well as the plotting time
    cFreq = args.freq
    iFreq = args.input_freq
    
    # Reader
    rdr = drx
    
    #
    # Step 1:  Read in the file and process it
    #
    
    # File size
    sizeB = os.path.getsize(filename)
    
    # Open the file
    fh = open(filename, 'rb')

    # Get the sampel rate and number of stands for each pol
    sample_rate = rdr.get_sample_rate(fh)
    nFramesFile = sizeB / rdr.FRAME_SIZE
    nFrames = numpy.array(rdr.get_frames_per_obs(fh))
    nCaptures = nFramesFile / nFrames.sum()

    # Number of remaining chunks
    maxFrames = (1936*10)
    nChunks = int(numpy.ceil(1.0*(nCaptures*nFrames.sum())/maxFrames))

    print("Filename:    %s" % filename)
    print("Size:        %.1f MB" % (float(sizeB)/1024/1024))
    print("Captures:    %i" % nCaptures)
    print("Sample Rate: %.2f kHz" % (sample_rate/1000.0))
    print("===")
    print("Chunks: %i" % nChunks)

    frame = rdr.read_frame(fh)
    fh.seek(0)
    
    standMapper = []
    for c in range(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFramesFile - c*maxFrames
        if framesRemaining > maxFrames:
            framesWork = maxFrames
        else:
            framesWork = framesRemaining
        print("Working on chunk %i, %i frames remaining" % (c+1, framesRemaining))
        
        count = {}
        data = numpy.zeros((nFrames.sum(),framesWork*frame.payload.data.size/nFrames.sum()), dtype=numpy.csingle)

        # Inner loop that actually reads the frames into the data array
        print("Working on %.1f ms of data" % ((framesWork*frame.payload.data.size/nFrames.sum()/sample_rate)*1000.0))

        # Go...
        dtime = numpy.zeros((nFrames.sum(), framesWork*frame.payload.data.size/nFrames.sum()), dtype=numpy.float64)
        data = numpy.zeros((nFrames.sum(), framesWork*frame.payload.data.size/nFrames.sum()), dtype=numpy.complex64)
        
        count = {}
        masterCount = 0
        tStart = 0
        data *= 0
        for f in xrange(framesWork):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = rdr.read_frame(fh, verbose=False)
            except errors.EOFError:
                break
            except errors.SyncError:
                #print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/rdr.FRAME_SIZE-1))
                continue
            
            if f == 0:
                tStart = cFrame.time

            try:
                stand,pol = cFrame.header.id
                aStand = 2*(stand-1)+pol
            except:
                beam,tune,pol = cFrame.id
                aStand = 4*(beam-1) + 2*(tune-1) + pol
                
            if aStand not in standMapper:
                standMapper.append(aStand)
                oStand = 1*aStand
                aStand = standMapper.index(aStand)
                try:
                    print("Mapping beam %i, tune. %1i, pol. %1i (%2i) to array index %3i" % (beam, tune, pol, oStand, aStand))
                except:
                    print("Mapping stand %i, pol. %1i (%2i) to array index %3i" % (stand, pol, oStand, aStand))
            else:
                aStand = standMapper.index(aStand)
            
            if aStand not in count.keys():
                count[aStand] = 0

            dtime[aStand, count[aStand]*cFrame.payload.data.size:(count[aStand]+1)*cFrame.payload.data.size] =  4096*count[aStand]/sample_rate + 1.0 / sample_rate * numpy.arange(0.0, cFrame.payload.data.size, dtype=numpy.float64)
            data[aStand, count[aStand]*cFrame.payload.data.size:(count[aStand]+1)*cFrame.payload.data.size] = cFrame.payload.data
            
            count[aStand] = count[aStand] + 1
            masterCount = masterCount + 1
                
        #
        # Step 2: Create the eye diagram
        #
        period = 1.0 / numpy.abs(iFreq - cFreq)
        dtime = dtime - dtime.min()
        
        endPt = data.shape[1]/8
        print(endPt / sample_rate / period)
    
        from scipy.optimize import fmin, leastsq
        
        fig = plt.figure()
        for i in xrange(data.shape[0]):
            if i == 0 or i == 3:
                continue
            
            #p0 = [6.0, iFreq-cFreq, numpy.pi/4]
            #p1 = fmin(qsDiff, p0, args=(dtime[i,:], data[i,:].real))
            #print(p1)
            
            freq = iFreq - cFreq
            def errFunc(p, x, y):
                scale = p[0]
                freq = p[1]
                phase = p[2]
                
                xPrime = 2*numpy.pi*freq*x + phase
                yFit = numpy.abs(scale)*numpy.sin(xPrime)
                
                return (yFit-y)
                
            p0 = [6, freq, 0.0]
            p1, success = leastsq(errFunc, p0, args=(dtime[i,:], data[i,:].real), maxfev=100000)
            #print(p0, p1, success)
            
            xPrime = 2*numpy.pi*p1[1]*dtime[i,:] + p1[2]
            #yFit = quantizedSine(xPrime, scale=p1[0])
            yFit = numpy.sin(xPrime)*numpy.abs(p1[0])
            print(i, standMapper[i], p1, ((data[i,:].real - yFit)**2).sum())
            
            ax1 = fig.add_subplot(2, 2, i+1)
            ax1.plot(dtime[i,:], data[i,:].real, color='blue')
            ax1.plot(dtime[i,:], yFit[:], color='green')
            ax1.plot(dtime[i,:], (yFit-data[i,:].real), color='red')
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="look for glitches in a DRX or TBN data by fitting a sine wave to the data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to analyze')
    parser.add_argument('-f', '--freq', type=aph.positive_or_zero_float, default=38.0,
                        help='frequency of the observations in MHz')
    parser.add_argument('-i', '--input-freq', type=aph.positive_or_zero_float, default=38.25,
                        help='frequency of the input sinusoid in MHz')
    parser.add_argument('-k', '--keep', type=aph.csv_int_list, default='1,2,3,4',
                        help='data array indiece (stands/beams) to keep')
    args = parser.parse_args()
    args.freq *= 1e6
    args.input_freq *= 1e6
    main(args)
    