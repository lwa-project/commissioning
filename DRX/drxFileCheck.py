#!/usr/bin/env python

"""
Run through a DRX file and determine if it is bad or not.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import ephem
import numpy
import argparse

from lsl import astro
from lsl.reader import drx, errors
from lsl.misc import parser as aph


def main(args):
    filename = args.filename
    
    fh = open(filename, "rb")
    nFramesFile = os.path.getsize(filename) / drx.FRAME_SIZE
    while True:
        try:
            junkFrame = drx.read_frame(fh)
            try:
                srate = junkFrame.sample_rate
                break
            except ZeroDivisionError:
                pass
        except errors.SyncError:
            fh.seek(-drx.FRAME_SIZE+1, 1)
            
    fh.seek(-drx.FRAME_SIZE, 1)
    
    beam, tune, pol = junkFrame.id
    tunepols = max(drx.get_frames_per_obs(fh))
    
    # Date & Central Frequnecy
    beginDate = ephem.Date(astro.unix_to_utcjd(sum(junkFrame.time)) - astro.DJD_OFFSET)
    central_freq1 = 0.0
    central_freq2 = 0.0
    for i in xrange(32):
        junkFrame = drx.read_frame(fh)
        b,t,p = junkFrame.id
        if p == 0 and t == 1:
            central_freq1 = junkFrame.central_freq
        elif p == 0 and t == 2:
            central_freq2 = junkFrame.central_freq
        else:
            pass
    fh.seek(-32*drx.FRAME_SIZE, 1)
    
    # Report on the file
    print("Filename: %s" % filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Beam: %i" % beam)
    print("Tune/Pols: %i" % tunepols)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (central_freq1, central_freq2))
    print(" ")
    
    # Convert chunk length to total frame count
    chunkLength = int(args.length * srate / 4096 * tunepols)
    chunkLength = int(1.0 * chunkLength / tunepols) * tunepols
    
    # Convert chunk skip to total frame count
    chunkSkip = int(args.skip * srate / 4096 * tunepols)
    chunkSkip = int(1.0 * chunkSkip / tunepols) * tunepols
    
    # Output arrays
    clipFraction = []
    meanPower = []
    
    # Go!
    i = 1
    done = False
    print("   |           Clipping              |          Power          |")
    print("   |      1X      1Y      2X      2Y |    1X    1Y    2X    2Y |")
    print("---+---------------------------------+-------------------------+")
    
    while True:
        count = {0:0, 1:0, 2:0, 3:0}
        data = numpy.empty((4,chunkLength*4096/tunepols), dtype=numpy.csingle)
        for j in xrange(chunkLength):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = drx.read_frame(fh, verbose=False)
            except errors.EOFError:
                done = True
                break
            except errors.SyncError:
                continue
            
            beam,tune,pol = cFrame.id
            aStand = 2*(tune-1) + pol
            
            try:
                data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.payload.data
                
                # Update the counters so that we can average properly later on
                count[aStand] += 1
            except ValueError:
                pass
        
        if done:
            break
            
        else:
            data = numpy.abs(data)**2
            data = data.astype(numpy.int32)
            
            clipFraction.append( numpy.zeros(4) )
            meanPower.append( data.mean(axis=1) )
            for j in xrange(4):
                bad = numpy.nonzero(data[j,:] > args.trim_level)[0]
                clipFraction[-1][j] = 1.0*len(bad) / data.shape[1]
            
            clip = clipFraction[-1]
            power = meanPower[-1]
            print("%2i | %6.2f%% %6.2f%% %6.2f%% %6.2f%% | %5.2f %5.2f %5.2f %5.2f |" % (i, clip[0]*100.0, clip[1]*100.0, clip[2]*100.0, clip[3]*100.0, power[0], power[1], power[2], power[3]))
        
            i += 1
            fh.seek(drx.FRAME_SIZE*chunkSkip, 1)
            
    clipFraction = numpy.array(clipFraction)
    meanPower = numpy.array(meanPower)
    
    clip = clipFraction.mean(axis=0)
    power = meanPower.mean(axis=0)
    
    print("---+---------------------------------+-------------------------+")
    print("%2s | %6.2f%% %6.2f%% %6.2f%% %6.2f%% | %5.2f %5.2f %5.2f %5.2f |" % ('M', clip[0]*100.0, clip[1]*100.0, clip[2]*100.0, clip[3]*100.0, power[0], power[1], power[2], power[3]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run through a DRX file and determine if it is bad or not.', 
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
    
