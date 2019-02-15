#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run through a TBN file and determine if it is bad or not.

$Rev: 965 $
$LastChangedBy: jdowell $
$LastChangedDate: 2012-07-25 17:33:24 -0600 (Wed, 25 Jul 2012) $
"""

import os
import sys
import ephem
import numpy
import argparse

from lsl import astro
from lsl.common import stations
from lsl.reader import tbn, errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.misc import parser as aph


def main(args):
    filename = args.filename
    
    # Set the station
    if args.lwasv:
        station = stations.lwasv
    else:
        station = stations.lwa1
    antennas = station.getAntennas()

    fh = open(filename, "rb")
    nFramesFile = os.path.getsize(filename) / tbn.FrameSize
    srate = tbn.getSampleRate(fh)
    antpols = len(antennas)

    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbn.readFrame(fh)
    fh.seek(-tbn.FrameSize, 1)
    centralFreq = junkFrame.getCentralFreq()
    beginDate = ephem.Date(astro.unix_to_utcjd(junkFrame.getTime()) - astro.DJD_OFFSET)

    # File summary
    print "Filename: %s" % filename
    print "Date of First Frame: %s" % str(beginDate)
    print "Ant/Pols: %i" % antpols
    print "Sample Rate: %i Hz" % srate
    print "Tuning Frequency: %.3f Hz" % centralFreq
    print " "

    # Convert chunk length to total frame count
    chunkLength = int(args.length * srate / 512 * antpols)
    chunkLength = int(1.0 * chunkLength / antpols) * antpols
    
    # Convert chunk skip to total frame count
    chunkSkip = int(args.skip * srate / 512 * antpols)
    chunkSkip = int(1.0 * chunkSkip / antpols) * antpols
    
    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1])

    # Output arrays
    clipFraction = []
    meanPower = []
    
    # Find stands #10
    toUse = []
    for i in xrange(antpols):
        ant = antennas[i]
        if ant.stand.id == 10:
            toUse.append(i)

    # Go!
    i = 1
    done = False
    print "   |     Clipping    |        Power      |"
    print "   |   10X     10Y   |    10X      10Y   |"
    print "---+-----------------+-------------------+"
    
    while True:
        count = [0 for j in xrange(antpols)]
        data = numpy.zeros((antpols, chunkLength*512/antpols), dtype=numpy.csingle)
        for j in xrange(chunkLength):
            try:
                cFrame = tbn.readFrame(fh)
            except errors.eofError:
                done = True
                break
            except errors.syncError:
                continue
                    
            buffer.append(cFrame)
            cFrames = buffer.get()

            if cFrames is None:
                continue
            
            for cFrame in cFrames:
                stand,pol = cFrame.header.parseID()
                
                # In the current configuration, stands start at 1 and go up to 260.  So, we
                # can use this little trick to populate the data array
                aStand = 2*(stand-1)+pol
                
                try:
                    data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
                    # Update the counters so that we can average properly later on
                    count[aStand] = count[aStand] + 1
                except ValueError:
                    pass
        # Empty the remaining portion of the buffer
        for cFrames in buffer.flush():
            # Inner loop that actually reads the frames into the data array
            for cFrame in cFrames:
                stand,pol = cFrame.header.parseID()
                # In the current configuration, stands start at 1 and go up to 10.  So, we
                # can use this little trick to populate the data array
                aStand = 2*(stand-1)+pol
            
                try:
                    data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
                    # Update the counters so that we can average properly later on
                    count[aStand] = count[aStand] + 1
                except ValueError:
                    pass

        if done:
            break
            
        else:
            data = numpy.abs(data)**2
            data = data.astype(numpy.int32)
            
            clipFraction.append( numpy.zeros(antpols) )
            meanPower.append( data.mean(axis=1) )
            for j in xrange(antpols):
                bad = numpy.nonzero(data[j,:] > args.trim_level)[0]
                clipFraction[-1][j] = 1.0*len(bad) / data.shape[1]
            
            clip = clipFraction[-1]
            power = meanPower[-1]
            print "%2i | %6.2f%% %6.2f%% | %8.2f %8.2f |" % (i, clip[toUse[0]]*100.0, clip[toUse[1]]*100.0, power[toUse[0]], power[toUse[1]])
        
            i += 1
            fh.seek(tbn.FrameSize*chunkSkip, 1)

    clipFraction = numpy.array(clipFraction)
    meanPower = numpy.array(meanPower)
    
    clip = clipFraction.mean(axis=0)
    power = meanPower.mean(axis=0)
    
    print "---+-----------------+-------------------+"
    print "%2s | %6.2f%% %6.2f%% | %8.2f %8.2f |" % ('M', clip[toUse[0]]*100.0, clip[toUse[1]]*100.0, power[toUse[0]], power[toUse[1]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='run through a TBN file and determine if it is bad or not', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-v', '--lwasv', action='store_true', 
                        help='use mapping from LWA-SV instead of LWA1')
    parser.add_argument('-l', '--length', type=aph.positive_float, default=1.0, 
                        help='length of time in seconds to analyze')
    parser.add_argument('-s', '--skip', type=aph.positive_float, default=900.0, 
                        help='skip period in seconds between chunks')
    parser.add_argument('-t', '--trim-level', type=aph.positive_float, default=16129.0, 
                        help='trim level for power analysis with clipping')
    args = parser.parse_args()
    main(args)
    