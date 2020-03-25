#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a TBW file, look for the weird RFI bursts that we have been seeing.  The
bursts are likely 'microsparking' from the powerline.
"""

import os
import sys
import math
import numpy
import ephem
import getopt

from lsl.common import stations
from lsl.common.dp import fS
from lsl.reader import tbw, tbn
from lsl.reader import errors
from lsl.correlator import fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common.progress import ProgressBar

import matplotlib.pyplot as plt


def usage(exitCode=None):
    print """burstMovie.py - Read in TBW files and create a NPZ file of raw time 
data centered around saturation events.

Usage: burstMovie.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-m, --metadata              Name of SSMIF file to use for mappings
-q, --quiet                 Run burstMovie.py in silent mode
-t, --threshold             Minimum digitizer value to consider a burst 
                            (default = 2000)
-n, --no-movie              Do not create movie frames (burst-*.png)
"""

    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['SSMIF'] = ''
    config['maxFrames'] = 30000*260
    config['verbose'] = True
    config['threshold'] = 2000
    config['movie'] = True
    config['args'] = []

    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hm:qt:n", ["help", "metadata=", "quiet", "threshold=", "no-movie"])
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
        elif opt in ('-q', '--quiet'):
            config['verbose'] = False
        elif opt in ('-t', '--threshold'):
            config['threshold'] = int(value)
        elif opt in ('-n', '--no-movie'):
            config['movie'] = False
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
    if config['SSMIF'] != '':
        station = stations.parseSSMIF(config['SSMIF'])
    else:
        station = stations.lwa1
    antennas = station.getAntennas()

    # Make sure that the file chunk size contains is an integer multiple
    # of the FFT length so that no data gets dropped
    maxFrames = int(config['maxFrames'])
    # It seems like that would be a good idea, however...  TBW data comes one
    # capture at a time so doing something like this actually truncates data 
    # from the last set of stands for the first integration.  So, we really 
    # should stick with
    maxFrames = config['maxFrames']

    fh = open(config['args'][0], "rb")
    nFrames = os.path.getsize(config['args'][0]) / tbw.FrameSize
    dataBits = tbw.getDataBits(fh)
    # The number of ant/pols in the file is hard coded because I cannot figure out 
    # a way to get this number in a systematic fashion
    antpols = len(antennas)
    nChunks = int(math.ceil(1.0*nFrames/maxFrames))
    if dataBits == 12:
        nSamples = 400
    else:
        nSamples = 1200

    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbw.readFrame(fh)
    fh.seek(0)
    beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)

    # File summary
    print "Filename: %s" % config['args'][0]
    print "Date of First Frame: %s" % str(beginDate)
    print "Ant/Pols: %i" % antpols
    print "Sample Length: %i-bit" % dataBits
    print "Frames: %i" % nFrames
    print "Chunks: %i" % nChunks
    print "==="

    nChunks = 1

    # Skip over any non-TBW frames at the beginning of the file
    i = 0
    junkFrame = tbw.readFrame(fh)
    while not junkFrame.header.isTBW():
        try:
            junkFrame = tbw.readFrame(fh)
        except errors.syncError:
            fh.seek(0)
            while True:
                try:
                    junkFrame = tbn.readFrame(fh)
                    i += 1
                except errors.syncError:
                    break
            fh.seek(-2*tbn.FrameSize, 1)
            junkFrame = tbw.readFrame(fh)
        i += 1
    fh.seek(-tbw.FrameSize, 1)
    print "Skipped %i non-TBW frames at the beginning of the file" % i

    # Master loop over all of the file chunks
    for i in range(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFrames - i*maxFrames
        if framesRemaining > maxFrames:
            framesWork = maxFrames
        else:
            framesWork = framesRemaining
        print "Working on chunk %i, %i frames remaining" % ((i+1), framesRemaining)

        #
        # NOTE
        # Major change here from tbwSpectra.py/stationMaster.py.  We are only keeping
        # the first 30,000 frames of the TBW file since we don't really need to read
        # in all of it to find bursts
        #
        data = numpy.zeros((antpols, 30000*nSamples), dtype=numpy.int16)
        
        # Inner loop that actually reads the frames into the data array
        for j in range(framesWork):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = tbw.readFrame(fh)
            except errors.eofError:
                break
            except errors.syncError:
                print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbw.FrameSize-1)
                continue
            if not cFrame.header.isTBW():
                    continue
                    
            # Skip frames over 30,000 on all stands
            if cFrame.header.frameCount > 30000:
                continue
            
            stand = cFrame.header.parseID()
            # In the current configuration, stands start at 1 and go up to 10.  So, we
            # can use this little trick to populate the data array
            aStand = 2*(stand-1)
            if cFrame.header.frameCount % 5000 == 0 and config['verbose']:
                print "%3i -> %3i  %6.3f  %5i  %i" % (stand, aStand, cFrame.getTime(), cFrame.header.frameCount, cFrame.data.timeTag)

            # Actually load the data.  x pol goes into the even numbers, y pol into the 
            # odd numbers
            count = cFrame.header.frameCount - 1
            data[aStand,   count*nSamples:(count+1)*nSamples] = cFrame.data.xy[0,:]
            data[aStand+1, count*nSamples:(count+1)*nSamples] = cFrame.data.xy[1,:]
            
        # Compute the power
        data = numpy.abs(data)

        # We need to various time series data to be aligned so we need to do a 
        # correction for the cable delays.  Using the various Antenna instances, 
        # we create a array of delays (in seconds) and do everything relative to 
        # the longest delay. 
        #
        # After this, we can align the data from all of the antennas.
        delays = numpy.array([a.cable.delay(30e6) for a in antennas])
        delays = numpy.round(delays * fS).astype(numpy.int32)
        delays = delays.max() - delays
        
        alignedData = numpy.zeros((data.shape[0], data.shape[1]-delays.max()), dtype=data.dtype)
        for s in xrange(data.shape[0]):
            alignedData[s,:] = data[s,delays[s]:(delays[s]+alignedData.shape[1])]
        del(data)

        # Using the good antennas (Antenna.getStatus() == 33), we need to find the 
        # start of the RFI pulses by looking for "large" data values.  To make sure 
        # that we aren't getting stuck on a first partial burst, skip the first one
        # and use the second set of "large" data values found.
        #
        # I was using "large" as saturation/clipping (>= 2047), but the new lower 
        # value for the ARX gain makes me want to lower to something more like 
        #
        inOne = False
        first = 0
        status = numpy.array([ant.getStatus() for ant in antennas])
        good = numpy.where( status == 33 )[0]
        while first < alignedData.shape[1]:
            mv = alignedData[good,first].max()
            if mv >= config['threshold']:
                if not inOne:
                    first += 5000
                    inOne = True
                else:
                    break
            else:
                first += 1
        print "Second burst at %i" % first
        
        # Keep only what would be interesting (200 samples before and 2,800 samples
        # afterward) around the burst.  This corresponds to a time range from 1 
        # microsecond before the start of the pulse to 14 microseconds later.  Save
        # the aligned data snippet to a NPZ file.
        alignedData = alignedData[:,first-200:first+2800]
        standPos = numpy.array([[ant.stand.x, ant.stand.y, ant.stand.z] for ant in antennas])
        junk, basename = os.path.split(config['args'][0])
        shortname, ext = os.path.splitext(basename)
        numpy.savez('%s-burst.npz' % shortname, data=alignedData, ssmif=config['SSMIF'])

        # Make the movie (if needed)
        if config['movie']:
            if config['verbose']:
                print "Creating movie frames"
                pb = ProgressBar(max=alignedData.shape[1]/2)
            else:
                pb = None
                
            fig = plt.figure(figsize=(12,6))
            for i in xrange(0,alignedData.shape[1],2):
                fig.clf()
                axX = fig.add_subplot(1, 2, 1)
                axY = fig.add_subplot(1, 2, 2)
                
                colorsX = 1.0*(alignedData[0::2,i] + alignedData[0::2,i+1]) / 2
                colorsY = 1.0*(alignedData[1::2,i] + alignedData[1::2,i+1]) / 2

                axX.scatter(standPos[0::2,0], standPos[0::2,1], c=colorsX, s=40.0, vmin=0, vmax=config['threshold'])
                axY.scatter(standPos[1::2,0], standPos[1::2,1], c=colorsY, s=40.0, vmin=0, vmax=config['threshold'])
                
                ## Add the fence as a dashed line
                axX.plot([-59.827, 59.771, 60.148, -59.700, -59.827], 
                        [59.752, 59.864, -59.618, -59.948, 59.752], linestyle='--', color='k')
                axY.plot([-59.827, 59.771, 60.148, -59.700, -59.827], 
                        [59.752, 59.864, -59.618, -59.948, 59.752], linestyle='--', color='k')

                ## Add the shelter
                axX.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
                        [45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')
                axY.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
                        [45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')
                        
                axX.set_xlim([-65, 65])
                axX.set_ylim([-65, 65])
                axX.set_title('X pol. $\Delta$t = %.1f ns' % (1.0*i/fS*1e9,))
                axX.set_xlabel('$\Delta$ X [m]')
                axX.set_ylabel('$\Delta$ Y [m]')
                
                axY.set_xlim([-65, 65])
                axY.set_ylim([-65, 65])
                axY.set_title('Y pol.')
                axY.set_xlabel('$\Delta$ X [m]')
                #axY.set_ylabel('$\Delta$ Y [m]')

                fig.savefig('burst-%05i.png' % (i/2,))
                
                if pb is not None:
                    pb.inc(amount=1)
                    if pb.amount != 0 and pb.amount % 10 == 0:
                        sys.stdout.write(pb.show()+'\r')
                        sys.stdout.flush()
            del(fig)
            
            if pb is not None:
                sys.stdout.write(pb.show()+'\r')
                sys.stdout.write('\n')
                sys.stdout.flush()
            
            if config['verbose']:
                print "Creating movie"
            os.system("mencoder mf://burst-*.png -mf fps=20:type=png -ovc lavc -lavcopts vcodec=mpeg4:aspect=2/1 -o %s-burst.avi -ffourcc DX50 -vf scale=600:1200,expand=600:1200" % shortname)


if __name__ == "__main__":
    main(sys.argv[1:])

