#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple script for performing time series cross-correlation of TBN data for 
all stands relative to the outlier (#258).

This script differs from simpleFringe.py in that it is designed to deal with
multiple TBN captures (frequencies) that are stored in a single file.

Usage:
./simpleFringeDemux.py <TBN data file>
"""

import os
import sys
import ephem
import numpy
import getopt

from lsl.common.stations import lwa1
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common.paths import data as dataPath

from matplotlib import pyplot as plt

from collections import deque

import fringe
from multiStation import parseSSMIF


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",]


def usage(exitCode=None):
    print """simpleFringeDemux.py - Simple script for performing time series cross-correlation 
of TBN data for all stands relative to the outlier

Usage: simpleFringeDemux.py [OPTIONS] file

Options:
-h, --help            Display this help information
-m, --metadata        Name of SSMIF file to use for mappings
-a, --average         Integration time in seconds (default = 10)
-r, --reference	      Stand to use as a reference (default = 258)
-c, --clip            Clip level in sqrt(I*I + Q*Q) to use to exclude
                    samples in the time domain (default = 120)
"""

    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['SSMIF'] = None
    config['tInt'] = 10.0
    config['refStand'] = 258
    config['clipLevel'] = 120.0
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hm:a:r:c:", ["help", "metadata=", "average=", "reference=", "clip="])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-m', '--metadata'):
            config['SSMIF'] = value
        elif opt in ('-a', '--average'):
            config['tInt'] = float(value)
        elif opt in ('-r', '--reference'):
            config['refStand'] = int(value)
        elif opt in ('-c', '--clip'):
            config['clipLevel'] = float(value)
        else:
            assert False
            
    # Add in arguments
    config['args'] = args
    
    # Return configuration
    return config


def unitRead(fh, count=520, found={}):
    for i in xrange(count):
        frame = tbn.readFrame(fh)
        
        try:
            found[frame.data.timeTag].append(frame)
        except KeyError:
            found[frame.data.timeTag] = []
            found[frame.data.timeTag].append(frame)
        
    return found


def main(args):
    config = parseOptions(args)

    # The task at hand
    filename = config['args'][0]
    
    # The station
    if config['SSMIF'] is not None:
        site = parseSSMIF(config['SSMIF'])
        ssmifContents = open(config['SSMIF']).readlines()
    else:
        site = lwa1
        ssmifContents = open(os.path.join(dataPath, 'lwa1-ssmif.txt')).readlines()
    observer = site.getObserver()
    antennas = site.getAntennas()
    
    # The file's parameters
    fh = open(filename, 'rb')
    
    nFramesFile = os.path.getsize(filename) / tbn.FrameSize
    srate = tbn.getSampleRate(fh)
    antpols = len(antennas)
    fh.seek(0)
    if srate < 1000:
        fh.seek(len(antennas)*4*tbn.FrameSize)
        srate = tbn.getSampleRate(fh)
        antpols = len(antennas)
        fh.seek(len(antennas)*4*tbn.FrameSize)
        
    # Reference antenna
    ref = config['refStand']
    foundRef = False
    for i,a in enumerate(antennas):
        if a.stand.id == ref and a.pol == 0:
            refX = i
            foundRef = True
        elif a.stand.id == ref and a.pol == 1:
            refY = i
        else:
            pass
    if not foundRef:
        raise RuntimeError("Cannot file Stand #%i" % ref)
    
    # Integration time (seconds and frames)
    tInt = config['tInt']
    nFrames = int(round(tInt*srate/512*antpols))
    tInt = nFrames / antpols * 512 / srate
    
    # Total run length
    nChunks = int(1.0 * nFramesFile / antpols * 512 / srate / tInt)
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbn.readFrame(fh)
    fh.seek(-tbn.FrameSize, 1)
    startFC = junkFrame.header.frameCount
    try:
        centralFreq = junkFrame.getCentralFreq()
    except AttributeError:
        from lsl.common.dp import fS
        centralFreq = fS * junkFrame.header.secondsCount / 2**32
    beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
    
    observer.date = beginDate
    srcs = [ephem.Sun(),]
    for line in _srcs:
        srcs.append( ephem.readdb(line) )
    
    for i in xrange(len(srcs)):
        srcs[i].compute(observer)
        
        if srcs[i].alt > 0:
            print "source %s: alt %.1f degrees, az %.1f degrees" % (srcs[i].name, srcs[i].alt*180/numpy.pi, srcs[i].az*180/numpy.pi)

    # File summary
    print "Filename: %s" % filename
    print "Date of First Frame: %s" % str(beginDate)
    print "Ant/Pols: %i" % antpols
    print "Sample Rate: %i Hz" % srate
    print "Tuning Frequency: %.3f Hz" % centralFreq
    print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate)
    print "---"
    print "Integration: %.3f s (%i frames; %i frames per stand/pol)" % (tInt, nFrames, nFrames / antpols)
    print "Chunks: %i" % nChunks
    
    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1])
    
    # Create the phase average and times
    LFFT = 512
    times = numpy.zeros(nChunks, dtype=numpy.float64)
    simpleVis = numpy.zeros((nChunks, antpols), dtype=numpy.complex64)
    centralFreqs = numpy.zeros(nChunks, dtype=numpy.float64)
    
    # Go!
    k = 0
    for i in xrange(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFramesFile - k
        if framesRemaining > nFrames:
            framesWork = nFrames
            data = numpy.zeros((antpols, framesWork/antpols*512), dtype=numpy.complex64)
        else:
            framesWork = framesRemaining + antpols*buffer.nSegments
            data = numpy.zeros((antpols, framesWork/antpols*512), dtype=numpy.complex64)
        print "Working on chunk %i, %i frames remaining" % (i+1, framesRemaining)
        
        count = [0 for a in xrange(len(antennas))]
        
        j = 0
        fillsWork = framesWork / antpols
        # Inner loop that actually reads the frames into the data array
        done = False
        while j < fillsWork:
            cFrames = deque()
            for l in xrange(len(antennas)):
                try:
                    cFrames.append( tbn.readFrame(fh) )
                    k = k + 1
                except errors.eofError:
                    ## Exit at the EOF
                    done = True
                    break
                except errors.syncError:
                    #print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FrameSize-1)
                    ## Exit at the first sync error
                    done = True
                    break
                    
            buffer.append(cFrames)
            cFrames = buffer.get()

            if cFrames is None:
                continue
                
            for cFrame in cFrames:
                stand,pol = cFrame.header.parseID()
                
                # In the current configuration, stands start at 1 and go up to 260.  So, we
                # can use this little trick to populate the data array
                aStand = 2*(stand-1)+pol
                
                # Save the time
                if j == 0 and aStand == 0:
                    times[i] = cFrame.getTime()
                    try:
                        centralFreqs[i] = cFrame.getCentralFreq()
                    except AttributeError:
                        centralFreqs[i] = fS * cFrame.header.secondsCount / 2**32
                    if i > 0:
                        if centralFreqs[i] != centralFreqs[i-1]:
                            print "Frequency change from %.3f to %.3f MHz at chunk %i" % (centralFreqs[i-1]/1e6, centralFreqs[i]/1e6, i+1)
                
                data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
                
                # Update the counters so that we can average properly later on
                count[aStand] = count[aStand] + 1
            
            j += 1
            
            if done:
                break

        if done:
            break
        
        # Time-domain blanking and cross-correlation with the outlier
        simpleVis[i,:] = fringe.Simple(data, refX, refY, config['clipLevel'])
    
    fh.close()
    
    # Save the data
    outname = os.path.split(filename)[1]
    outname = os.path.splitext(outname)[0]
    outname = "%s-ref%03i-multi-vis.npz" % (outname, config['refStand'])
    numpy.savez(outname, ref=ref, refX=refX, refY=refY, tInt=tInt, centralFreqs=centralFreqs, times=times, 
            simpleVis=simpleVis, ssmifContents=ssmifContents)


if __name__ == "__main__":
    main(sys.argv[1:])
    
