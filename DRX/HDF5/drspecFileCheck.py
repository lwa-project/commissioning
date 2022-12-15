#!/usr/bin/env python3

"""
Run through a DR spectrometer file and determine if it is bad or not.
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    range = xrange
except NameError:
    pass
    
import os
import sys
import numpy
import argparse
from datetime import datetime

from lsl import astro
from lsl.reader import drx, drspec, errors
from lsl.misc import parser as aph


def main(args):
    fh = open(args.filename, "rb")
    
    try:
        for i in range(5):
            junkFrame = drx.read_frame(fh)
        raise RuntimeError("ERROR: '%s' appears to be a raw DRX file, not a DR spectrometer file" % args.filename)
    except errors.SyncError:
        fh.seek(0)
        
    # Interrogate the file to figure out what frames sizes to expect, now many 
    # frames there are, and what the transform length is
    FRAME_SIZE = drspec.get_frame_size(fh)
    nFrames = os.path.getsize(args.filename) / FRAME_SIZE
    nChunks = nFrames
    LFFT = drspec.get_transform_size(fh)
    
    # Read in the first frame to figure out the DP information
    junkFrame = drspec.read_frame(fh)
    fh.seek(-FRAME_SIZE, 1)
    srate = junkFrame.sample_rate
    t0 = junkFrame.time
    tInt = junkFrame.header.nints*LFFT/srate
    
    # Update the file contents
    beam = junkFrame.id
    central_freq1, central_freq2 = junkFrame.central_freq
    srate = junkFrame.sample_rate
    data_products = junkFrame.data_products
    t0 = junkFrame.time
    tInt = junkFrame.header.nints*LFFT/srate
    beginDate = junkFrame.time.datetime
        
    # Report
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % beginDate)
    print("Beam: %i" % beam)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (central_freq1, central_freq2))
    print("Data Products: %s" % ','.join(data_products))
    print("Frames: %i (%.3f s)" % (nFrames, nFrames*tInt))
    print("---")
    print("Transform Length: %i" % LFFT)
    print("Integration: %.3f s" % tInt)
    
    # Convert chunk length to total frame count
    chunkLength = int(args.length / tInt)
    
    # Convert chunk skip to total frame count
    chunkSkip = int(args.skip / tInt)
    
    # Output arrays
    clipFraction = []
    meanPower = []
    
    # Go!
    i = 1
    done = False
    print("   |%sClipping%s |%sPower %s |" % (" "*(8*len(data_products)-4), " "*(8*len(data_products)-4), " "*(6*len(data_products)-3), " "*(6*len(data_products)-3)))
    out = "   |      1X      1Y      2X      2Y |"
    for t in (1, 2):
        for dp in data_products:
            out += "%6s" % ("%i%s" % (t, dp))
    out += " |"
    print(out)
    print("-"*len(out))
    
    while True:
        count = {0:0, 1:0, 2:0, 3:0}
        sats = numpy.empty((4,chunkLength), dtype=numpy.float32)
        data = numpy.empty((2*len(data_products),chunkLength*LFFT), dtype=numpy.float32)
        for j in range(chunkLength):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = frame = drspec.read_frame(fh)
            except errors.EOFError:
                done = True
                break
            except errors.SyncError:
                continue
                
            for t in (1,2):
                for p,dp in enumerate(data_products):
                    l = len(data_products)*(t-1) + p
                    data[l,j*LFFT:(j+1)*LFFT] = getattr(cFrame.payload, '%s%i' % (dp,t-1))
            sats[:,j] = numpy.array(cFrame.payload.saturations) / (tInt*srate)
                    
        if done:
            break
            
        else:
            clipFraction.append( sats.mean(axis=1) )
            meanPower.append( data.mean(axis=1) )
            
            clip = clipFraction[-1]
            power = meanPower[-1]
            
            out = "%2i | %6.2f%% %6.2f%% %6.2f%% %6.2f%% |" % (i, clip[0]*100.0, clip[1]*100.0, clip[2]*100.0, clip[3]*100.0)
            for t in (1, 2):
                for p in range(len(data_products)):
                    out += " %5.2f" % (power[len(data_products)*(t-1)+p],)
            out += " |"
            print(out)
        
            i += 1
            fh.seek(FRAME_SIZE*chunkSkip, 1)
            
    clipFraction = numpy.array(clipFraction)
    meanPower = numpy.array(meanPower)
    
    clip = clipFraction.mean(axis=0)
    power = meanPower.mean(axis=0)
    
    print("-"*len(out))
    out = "%2s | %6.2f%% %6.2f%% %6.2f%% %6.2f%% |" % ('M', clip[0]*100.0, clip[1]*100.0, clip[2]*100.0, clip[3]*100.0)
    for t in (1, 2):
        for p in range(len(data_products)):
            out += " %5.2f" % (power[len(data_products)*(t-1)+p],)
    out += " |"
    print(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run through a DR spectrometer file and determine if it is bad or not.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-l', '--length', type=aph.positive_float, default=1.0, 
                        help='length of time in seconds to analyze')
    parser.add_argument('-s', '--skip', type=aph.positive_float, default=900.0, 
                        help='skip period in seconds between chunks')
    parser.add_argument('-t', '--trim-level', type=aph.positive_float, default=49, 
                        help='trim level for power analysis with clipping')
    args = parser.parse_args()
    main(args)
    
