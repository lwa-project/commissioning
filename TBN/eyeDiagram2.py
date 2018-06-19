#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Look for glitches in a DRX or TBN data by fitting a sine wave to the data.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy

import getopt

from lsl.reader import errors, tbn, drx

from matplotlib import pyplot as plt


def quantizedSine(x, scale=1.0):
    return numpy.round(scale*numpy.round(2**11*numpy.sin(x)))


def qsDiff(B, x, y):
    scale = B[0]
    freq = B[1]
    phase = B[2]
    
    yFit = quantizedSine(2*numpy.pi*freq + phase, scale=scale)
    
    #return yFit - y
    return ((yFit - y)**2).sum()


def usage(exitCode=None):
    print """eyeDiagram2.py - Look for glitches in a DRX or TBN data by fitting a 
sine wave to the data.

Usage:  eyeDiagram2.py [OPTIONS] data_file

Options:
-h, --help                  Display this help information
-d, --drx                   DRX data is being supplied (default = TBN)
-f, --freq                  Set the frequency of the observations in MHz (default = 
                            38.00 MHz)
-i, --input-freq            Set the frequency of the input sinusoid in MHz (default = 
                            38.25 MHz)
-k, --keep                  Data array indiece (stands/beams) to keep (default = 1,2,3,4)
"""

    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['reader'] = tbn
    config['maxFrames'] = 1936
    config['freq'] = 38.00e6
    config['ifreq'] = 38.25e6
    config['keep'] = [1,2,3,4]
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hdf:i:k:", ["help", "drx", "freq=", "input-freq=", "keep="])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-d', '--drx'):
            config['reader'] = drx
        elif opt in ('-f', '--freq'):
            config['freq'] = float(value)*1e6
        elif opt in ('-i', '--input-freq'):
            config['ifreq'] = float(value)*1e6
        elif opt in ('-k', '--keep'):
            config['keep'] = value.split(',')
        else:
            assert False
            
    # Add in arguments
    config['args'] = args

    # Return configuration
    return config


def main(args):
    # Parse command line options
    config = parseOptions(args)
    
    # Filename
    filename = config['args'][0]
    
    # The observation and sinusoid frequencies as well as the plotting time
    cFreq = config['freq']
    iFreq = config['ifreq']
    
    # Reader
    rdr = config['reader']
    
    #
    # Step 1:  Read in the file and process it
    #
    
    # File size
    sizeB = os.path.getsize(filename)
    
    # Open the file
    fh = open(filename, 'rb')

    # Get the sampel rate and number of stands for each pol
    sampleRate = rdr.getSampleRate(fh)
    nFramesFile = sizeB / rdr.FrameSize
    nFrames = numpy.array(rdr.getFramesPerObs(fh))
    nCaptures = nFramesFile / nFrames.sum()

    # Number of remaining chunks
    maxFrames = config['maxFrames']
    nChunks = int(numpy.ceil(1.0*(nCaptures*nFrames.sum())/maxFrames))

    print "Filename:    %s" % filename
    print "Size:        %.1f MB" % (float(sizeB)/1024/1024)
    print "Captures:    %i" % nCaptures
    print "Sample Rate: %.2f kHz" % (sampleRate/1000.0)
    print "==="
    print "Chunks: %i" % nChunks

    frame = rdr.readFrame(fh)
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
        print "Working on chunk %i, %i frames remaining" % (c+1, framesRemaining)
        
        count = {}
        data = numpy.zeros((nFrames.sum(),framesWork*frame.data.iq.size/nFrames.sum()), dtype=numpy.csingle)

        # Inner loop that actually reads the frames into the data array
        print "Working on %.1f ms of data" % ((framesWork*frame.data.iq.size/nFrames.sum()/sampleRate)*1000.0)

        # Go...
        dtime = numpy.zeros((nFrames.sum(), framesWork*frame.data.iq.size/nFrames.sum()), dtype=numpy.float64)
        data = numpy.zeros((nFrames.sum(), framesWork*frame.data.iq.size/nFrames.sum()), dtype=numpy.complex64)
        
        count = {}
        masterCount = 0
        tStart = 0
        data *= 0
        for f in xrange(framesWork):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = rdr.readFrame(fh, Verbose=False)
            except errors.eofError:
                break
            except errors.syncError:
                #print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/rdr.FrameSize-1)
                continue
            
            if f == 0:
                tStart = cFrame.getTime()

            try:
                stand,pol = cFrame.header.parseID()
                aStand = 2*(stand-1)+pol
            except:
                beam,tune,pol = cFrame.parseID()
                aStand = 4*(beam-1) + 2*(tune-1) + pol
                
            if aStand not in standMapper:
                standMapper.append(aStand)
                oStand = 1*aStand
                aStand = standMapper.index(aStand)
                try:
                    print "Mapping beam %i, tune. %1i, pol. %1i (%2i) to array index %3i" % (beam, tune, pol, oStand, aStand)
                except:
                    print "Mapping stand %i, pol. %1i (%2i) to array index %3i" % (stand, pol, oStand, aStand)
            else:
                aStand = standMapper.index(aStand)
            
            if aStand not in count.keys():
                count[aStand] = 0

            dtime[aStand, count[aStand]*cFrame.data.iq.size:(count[aStand]+1)*cFrame.data.iq.size] = cFrame.getTime() + 1.0 / sampleRate * numpy.arange(0.0, cFrame.data.iq.size, dtype=numpy.float64)
            data[aStand, count[aStand]*cFrame.data.iq.size:(count[aStand]+1)*cFrame.data.iq.size] = cFrame.data.iq
            
            count[aStand] = count[aStand] + 1
            masterCount = masterCount + 1
                
        #
        # Step 2: Create the eye diagram
        #
        period = 1.0 / numpy.abs(iFreq - cFreq)
        dtime = dtime - dtime.min()
        
        endPt = data.shape[1]/8
        print endPt / sampleRate / period
    
        from scipy.optimize import fmin, leastsq
        p0 = [data.max(), iFreq-cFreq, 0.0]
        for i in xrange(data.shape[0]):
            freq = numpy.fft.fftfreq(4096, d=1/sampleRate)
            psd = numpy.abs(numpy.fft.fft(data[i,0:4096]))**2
            print freq[numpy.where( psd == psd.max() )[0]]
            #import pylab
            #pylab.plot(freq, numpy.log10(psd)*10)
            #pylab.show()
            
            p0 = [data[i,:].std()*1.4/2048, freq[numpy.where( psd == psd.max() )][0], 0.0]
            
            p1 = fmin(qsDiff, p0, args=(dtime[i,:], data[i,:].real))
            print p1
            
            xPrime = 2*numpy.pi*p1[1]*dtime[i,:] + p1[2]
            yFit = quantizedSine(xPrime, scale=p1[0])
            import pylab
            #pylab.plot(dtime[i,0:100], (yFit-data[i,:].real)[0:100])
            pylab.plot(dtime[i,0:1000], data[i,0:1000].real, color='blue')
            #pylab.plot(dtime[i,0:1000], yFit[0:1000], color='green')
            pylab.show()
            print data[i,:].real - yFit
            

if __name__ == "__main__":
    main(sys.argv[1:])
    
