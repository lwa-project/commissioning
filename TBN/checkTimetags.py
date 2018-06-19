#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Check the time tags in a full 520 antenna stand data set.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import ephem
import numpy
import getopt

from lsl import astro
from lsl.common import stations
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.common.dp import fS


def usage(exitCode=None):
    print """checkTimetags.py - Read in a TBN file and check the flow of time.

Usage: checkTimetags.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-v, --lwasv                 Use mapping from LWA-SV instead of LWA1
"""

    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['site'] = 'lwa1'
    config['args'] = []

    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hv", ["help", "lwasv"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-v', '--lwasv'):
            config['site'] = 'lwasv'
        else:
            assert False
    
    # Add in arguments
    config['args'] = args

    # Return configuration
    return config


def main(args):
    # Parse command line options
    config = parseOptions(args)
    
    # Set the station
    if config['site'] == 'lwa1':
        station = stations.lwa1
    elif config['site'] == 'lwasv':
        station = stations.lwasv
    else:
        raise RuntimeError("Unknown site name: %s" % config['site'])
    antennas = station.getAntennas()
    
    fh = open(config['args'][0], "rb", buffering=tbn.FrameSize*10000)

    # Get the first frame and find out what the firt time tag is, which the
    # first frame number is, and what the sample rate it.  From the sample 
    # rate, estimate how the time tag should advance between frames.
    junkFrame = tbn.readFrame(fh)
    sampleRate = tbn.getSampleRate(fh)
    antpols = len(antennas)
    tagSkip = fS / sampleRate * junkFrame.data.iq.shape[0]
    fh.seek(0)

    # Store the information about the first frame and convert the timetag to 
    # an ephem.Date object.
    prevTime = junkFrame.data.timeTag
    prevDate = ephem.Date(astro.unix_to_utcjd(junkFrame.getTime()) - astro.DJD_OFFSET)
    prevFrame = junkFrame.header.frameCount

    # Report on the file
    print "Filename: %s" % os.path.basename(config['args'][0])
    print "Date of first frame: %i -> %s" % (prevTime, str(prevDate))
    print "Sample rate: %i Hz" % sampleRate
    print "Time tag skip per frame: %i" % tagSkip

    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1])
    
    j = 0
    k = 0
    while True:
        try:
            cFrame = tbn.readFrame(fh)
            k += 1
        except errors.eofError:
            break
        except errors.syncError:
            #print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FrameSize-1)
            continue
                
        buffer.append(cFrame)
        cFrames = buffer.get()

        if cFrames is None:
            continue
        
        valid = reduce(lambda x,y: x+int(y.valid), cFrames, 0)
        if valid != antpols:
            print "WARNING: frame count %i at %i missing %.2f%% of frames" % (cFrames[0].header.frameCount, cFrames[0].data.timeTag, float(antpols - valid)/antpols*100)
        
        timeTags = numpy.zeros(len(cFrames), dtype=numpy.int64) - 1
        for cFrame in cFrames:
            stand,pol = cFrame.parseID()
            timeTags[2*(stand-1)+pol] = cFrame.data.timeTag
            
        if j == 0:
            prevTime  = numpy.median(timeTags)
            prevDate  = ephem.Date(astro.unix_to_utcjd(cFrames[0].getTime()) - astro.DJD_OFFSET)
            prevFrame = cFrames[0].header.frameCount
            
            j += 1
            continue
        else:
            currTime = numpy.median(timeTags)
            currDate  = ephem.Date(astro.unix_to_utcjd(cFrames[0].getTime()) - astro.DJD_OFFSET)
            currFrame = cFrames[0].header.frameCount
            
        if currFrame % 1000 == 0:
            print "At frame %i t.t. is %i -> %s" % (currFrame, currTime, currDate)

        if currTime < prevTime:
            print "ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i" % (currTime, currFrame, prevTime, prevFrame)
            print "       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate))
        if (currTime-prevTime) > tagSkip:
            print "ERROR: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (currTime, currFrame, prevTime, prevFrame)
            print "       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate))
        for i in xrange(timeTags.size):
            if timeTags[i] != currTime:
                print "ERROR: t.t. of dig. %i != frame set median of %i" % (i, currTime)
                print "       -> difference: %i" % (currTime-timeTags[i],)
        
        prevTime  = currTime
        prevData  = currDate
        prevFrame = currFrame
        
        j += 1


if __name__ == "__main__":
    main(sys.argv[1:])
