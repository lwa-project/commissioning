#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple script (with hard coded integration times) for performing time series
cross-correlation of TBN data for all stands relative to the outlier (#258).

Usage:
./simpleFringe.py <TBN data file>
"""

import os
import sys
import ephem
import numpy
import getopt

from lsl.common.stations import parse_ssmif, lwa1
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common.paths import data as dataPath

from matplotlib import pyplot as plt


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",]


def usage(exitCode=None):
    print """simpleFringe.py - Simple script for performing time series cross-correlation 
of TBN data for all stands relative to the outlier

Usage: simpleFringe.py [OPTIONS] file

Options:
-h, --help            Display this help information
-m, --metadata        Name of SSMIF file to use for mappings
-a, --average         Integration time in seconds (default = 10)
-r, --reference	      Stand to use as a reference (default = 258)
-c, --clip            Clip level in sqrt(I*I + Q*Q) to use to exclude
                    samples in the time domain (default = 0 = no 
                    excision)
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
    config['clip_level'] = 0
    
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
            config['clip_level'] = float(value)
        else:
            assert False
    
    # Add in arguments
    config['args'] = args

    # Return configuration
    return config


def main(args):
    config = parseOptions(args)

    # The task at hand
    filename = config['args'][0]
    
    # The station
    if config['SSMIF'] is not None:
        site = parse_ssmif(config['SSMIF'])
        ssmifContents = open(config['SSMIF']).readlines()
    else:
        site = lwa1
        ssmifContents = open(os.path.join(dataPath, 'lwa1-ssmif.txt')).readlines()
    observer = site.get_observer()
    antennas = site.antennas
    
    # The file's parameters
    fh = open(filename, 'rb')
    nFramesFile = os.path.getsize(filename) / tbn.FRAME_SIZE
    srate = tbn.get_sample_rate(fh)
    antpols = len(antennas)
    
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
    nChunks = 1
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbn.read_frame(fh)
    fh.seek(-tbn.FRAME_SIZE, 1)
    startFC = junkFrame.header.frame_count
    try:
        central_freq = junkFrame.central_freq
    except AttributeError:
        from lsl.common.dp import fS
        central_freq = fS * junkFrame.header.second_count / 2**32
    beginDate = ephem.Date(unix_to_utcjd(junkFrame.get_time()) - DJD_OFFSET)
    
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
    print "Tuning Frequency: %.3f Hz" % central_freq
    print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate)
    print "---"
    print "Integration: %.3f s (%i frames; %i frames per stand/pol)" % (tInt, nFrames, nFrames / antpols)
    print "Chunks: %i" % nChunks
    
    junkFrame = tbn.read_frame(fh)
    while junkFrame.header.frame_count < startFC+3:
        junkFrame = tbn.read_frame(fh)
    fh.seek(-tbn.FRAME_SIZE, 1)
    
    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1])
    
    # Create the phase average and times
    LFFT = 512
    times = numpy.zeros(nChunks, dtype=numpy.float64)
    fullVis   = numpy.zeros((nChunks, antpols, LFFT), dtype=numpy.complex64)
    simpleVis = numpy.zeros((nChunks, antpols), dtype=numpy.complex64)
    
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
            framesWork = framesRemaining + antpols*buffer.nsegments
            data = numpy.zeros((antpols, framesWork/antpols*512), dtype=numpy.complex64)
        print "Working on chunk %i, %i frames remaining" % (i+1, framesRemaining)
        
        count = [0 for a in xrange(len(antennas))]
        
        j = 0
        fillsWork = framesWork / antpols
        # Inner loop that actually reads the frames into the data array
        while j < fillsWork:
            try:
                cFrame = tbn.read_frame(fh)
                k = k + 1
            except errors.EOFError:
                break
            except errors.SyncError:
                #print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FRAME_SIZE-1)
                continue
                    
            buffer.append(cFrame)
            cFrames = buffer.get()

            if cFrames is None:
                continue
                
            for cFrame in cFrames:
                stand,pol = cFrame.header.id
                
                # In the current configuration, stands start at 1 and go up to 260.  So, we
                # can use this little trick to populate the data array
                aStand = 2*(stand-1)+pol
                
                # Save the time
                if j == 0 and aStand == 0:
                    times[i] = cFrame.get_time()
                
                data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
                
                # Update the counters so that we can average properly later on
                count[aStand] = count[aStand] + 1
            
            j += 1
            
        # Mask
        if config['clip_level'] > 0:
            bad = numpy.where( numpy.abs(data) >= config['clip_level'] )
            data[bad] *= 0.0
        
        # Simple correlation
        for l in xrange(520):
            if l % 2 == 0:
                simpleVis[i,l] = (data[l,:]*data[refX,:].conj()).mean()
            else:
                simpleVis[i,l] = (data[l,:]*data[refY,:].conj()).mean()
    
    # Save the data
    outname = os.path.split(filename)[1]
    outname = os.path.splitext(outname)[0]
    outname = "%s-ref%03i-vis.npz" % (outname, config['refStand'])
    numpy.savez(outname, ref=ref, refX=refX, refY=refY, tInt=tInt, central_freq=central_freq, times=times, 
            fullVis=fullVis, simpleVis=simpleVis, ssmifContents=ssmifContents)


if __name__ == "__main__":
    main(sys.argv[1:])
    
