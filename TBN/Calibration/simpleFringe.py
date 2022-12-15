#!/usr/bin/env python3

"""
Simple script (with hard coded integration times) for performing time series
cross-correlation of TBN data for all stands relative to the outlier (#258).

Usage:
./simpleFringe.py <TBN data file>
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    range = xrange
except NameError:
    pass
    
import os
import sys
import ephem
import numpy
import argparse

from lsl.common.stations import parse_ssmif, lwa1
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common.paths import DATA as dataPath
from lsl.misc import parser as aph

from matplotlib import pyplot as plt


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",]


def main(args):
    # The task at hand
    filename = args.filename
    
    # The station
    if args.metadata is not None:
        site = parse_ssmif(args.metadata)
        ssmifContents = open(args.metadata).readlines()
    else:
        site = lwa1
        ssmifContents = open(os.path.join(dataPath, 'lwa1-ssmif.txt')).readlines()
    observer = site.get_observer()
    antennas = site.antennas
    
    # The file's parameters
    fh = open(filename, 'rb')
    nFramesFile = os.path.getsize(filename) // tbn.FRAME_SIZE
    srate = tbn.get_sample_rate(fh)
    antpols = len(antennas)
    
    # Reference antenna
    ref = args.reference
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
    tInt = args.average
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
    beginDate = junkFrame.time.datetime
    
    observer.date = beginDate
    srcs = [ephem.Sun(),]
    for line in _srcs:
        srcs.append( ephem.readdb(line) )
    
    for i in range(len(srcs)):
        srcs[i].compute(observer)
        
        if srcs[i].alt > 0:
            print("source %s: alt %.1f degrees, az %.1f degrees" % (srcs[i].name, srcs[i].alt*180/numpy.pi, srcs[i].az*180/numpy.pi))

    # File summary
    print("Filename: %s" % filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Ant/Pols: %i" % antpols)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz" % central_freq)
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate))
    print("---")
    print("Integration: %.3f s (%i frames; %i frames per stand/pol)" % (tInt, nFrames, nFrames // antpols))
    print("Chunks: %i" % nChunks)
    
    junkFrame = tbn.read_frame(fh)
    while junkFrame.header.frame_count < startFC+3:
        junkFrame = tbn.read_frame(fh)
    fh.seek(-tbn.FRAME_SIZE, 1)
    
    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols//2+1), pols=[0, 1])
    
    # Create the phase average and times
    LFFT = 512
    times = numpy.zeros(nChunks, dtype=numpy.float64)
    fullVis   = numpy.zeros((nChunks, antpols, LFFT), dtype=numpy.complex64)
    simpleVis = numpy.zeros((nChunks, antpols), dtype=numpy.complex64)
    
    # Go!
    k = 0
    for i in range(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFramesFile - k
        if framesRemaining > nFrames:
            framesWork = nFrames
            data = numpy.zeros((antpols, framesWork//antpols*512), dtype=numpy.complex64)
        else:
            framesWork = framesRemaining + antpols*buffer.nsegments
            data = numpy.zeros((antpols, framesWork//antpols*512), dtype=numpy.complex64)
        print("Working on chunk %i, %i frames remaining" % (i+1, framesRemaining))
        
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
                
            for cFrame in cFrames:
                stand,pol = cFrame.header.id
                
                # In the current configuration, stands start at 1 and go up to 260.  So, we
                # can use this little trick to populate the data array
                aStand = 2*(stand-1)+pol
                
                # Save the time
                if j == 0 and aStand == 0:
                    times[i] = cFrame.time
                
                data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
                
                # Update the counters so that we can average properly later on
                count[aStand] = count[aStand] + 1
            
            j += 1
            
        # Mask
        if args.clip > 0:
            bad = numpy.where( numpy.abs(data) >= args.clip )
            data[bad] *= 0.0
        
        # Simple correlation
        for l in range(520):
            if l % 2 == 0:
                simpleVis[i,l] = (data[l,:]*data[refX,:].conj()).mean()
            else:
                simpleVis[i,l] = (data[l,:]*data[refY,:].conj()).mean()
    
    # Save the data
    outname = os.path.split(filename)[1]
    outname = os.path.splitext(outname)[0]
    outname = "%s-ref%03i-vis.npz" % (outname, args.reference)
    numpy.savez(outname, ref=ref, refX=refX, refY=refY, tInt=tInt, central_freq=central_freq, times=times, 
            fullVis=fullVis, simpleVis=simpleVis, ssmifContents=ssmifContents)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="simple script for performing time series cross-correlation of TBN data for all stands relative to the outlier",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str,
                        help='file to fringe')
    parser.add_argument('-m', '--metadata', type=str,
                        help='name of SSMIF file to use for mappings')
    parser.add_argument('-a', '--average', type=aph.positive_float, default=10.0,
                        help='integration time in seconds')
    parser.add_argument('-r', '--reference', type=aph.positive_int, default=258,
                        help='stand to use as a reference')
    parser.add_argument('-c', '--clip', type=aph.positive_or_zero_float, default=0.0,
                        help='clip level in sqrt(I*I + Q*Q) to use to exclude samples in the time domain; 0 = no excision')
    args = parser.parse_args()
    main(args)
    
