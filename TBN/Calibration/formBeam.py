#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple TBN beam forming script based on Steve's "Fun with TBN" memo.

Usage:
./formBeam.py <cln_file> <azimuth> <elevation> <TBN_file>

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import ephem
import numpy

from lsl.common.constants import c as vLight
from lsl.common.stations import lwa1
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.astro import unix_to_utcjd, DJD_OFFSET

from matplotlib import pyplot as plt

from multiprocessing import Pool


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",]
        
        
def formBeam(data, bln):
    """
    Apply the beam forming coefficients to a section of TBN data, sum across
    the various inputs, and intergrate.
    """
    
    temp = 1.0*data
    for i in xrange(bln.size):
        temp[i,:] *= bln[i]
    
    temp = temp.sum(axis=0)
    
    return (numpy.abs(temp)**2).mean()


def getGeoDelay(antenna, az, el, freq, Degrees=False):
    """
    Get the geometrical delay (relative to the center of the array)
    for the specified antenna for a source at azimuth az, elevation el.
    """

    if Degrees:
        az = az*numpy.pi/180.0
        el = el*numpy.pi/180.0
    
    source = numpy.array([numpy.cos(el)*numpy.sin(az), 
                    numpy.cos(el)*numpy.cos(az), 
                    numpy.sin(el)])
    
    cableDelay = antenna.cable.delay(freq)
    xyz = numpy.array([antenna.stand.x, antenna.stand.y, antenna.stand.z])
    return numpy.dot(source, xyz) / vLight - 0*cableDelay


def main(args):
    # The task at hand
    clnfile = args[0]
    az = float(args[1])
    el = float(args[2])
    filename = args[3]
    
    # The station
    observer = lwa1.getObserver()
    antennas = lwa1.getAntennas()
    
    # The file's parameters
    fh = open(filename, 'rb')
    nFramesFile = os.path.getsize(filename) / tbn.FrameSize
    srate = tbn.getSampleRate(fh)
    antpols = len(antennas)
    
    # Reference antenna
    ref = 258
    for a in antennas:
        if a.stand.id == ref and a.pol == 0:
            refX = a.digitizer
        elif a.stand.id == ref and a.pol == 1:
            refY = a.digitizer
        else:
            pass
    
    # Integration time (seconds and frames)
    tInt = 5.0
    nFrames = int(round(tInt*srate/512*antpols))
    tInt = nFrames / antpols * 512 / srate
    
    # Total run length
    #nChunks = int(round(1.0*nFramesFile / nFrames))
    nChunks = 240
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbn.readFrame(fh)
    fh.seek(-tbn.FrameSize, 1)
    startFC = junkFrame.header.frameCount
    centralFreq = junkFrame.getCentralFreq()
    beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
    
    observer.date = beginDate
    srcs = []
    for line in _srcs:
        srcs.append( ephem.readdb(line) )
        srcs[-1].compute(observer)
        
        if srcs[-1].alt > 0:
            print "source %s: alt %.1f degrees, az %.1f degrees" % (srcs[-1].name, srcs[-1].alt*180/numpy.pi, srcs[-1].az*180/numpy.pi)

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
    
    junkFrame = tbn.readFrame(fh)
    while junkFrame.header.frameCount < startFC+3:
        junkFrame = tbn.readFrame(fh)
    fh.seek(-tbn.FrameSize, 1)
    
    # Get the beamformer coefficients - three sets:
    #  (1) at the requested az, el
    #  (2) at az, el - 15 degrees
    #  (3) at the transit location of Cyg A
    dataDict = numpy.load(clnfile)
    cln = dataDict['cln']
    aln1 = []
    aln2 = []
    aln3 = []
    for i in xrange(cln.shape[1]):
        gd = getGeoDelay(antennas[i], az, el, centralFreq, Degrees=True)
        aln1.append( numpy.exp(2j*numpy.pi*centralFreq*gd) )
        
        gd = getGeoDelay(antennas[i], az, el-15, centralFreq, Degrees=True)
        aln2.append( numpy.exp(2j*numpy.pi*centralFreq*gd) )
        
        gd = getGeoDelay(antennas[i], 0.5, 83.3, centralFreq, Degrees=True)
        aln3.append( numpy.exp(2j*numpy.pi*centralFreq*gd) )
        
    aln1 = numpy.array(aln1)
    aln2 = numpy.array(aln2)
    aln3 = numpy.array(aln3)
    
    bln1 = (cln*aln1).conj() / numpy.abs(cln*aln1)
    bln2 = (cln*aln2).conj() / numpy.abs(cln*aln2)
    bln3 = (cln*aln3).conj() / numpy.abs(cln*aln3)
    for i in xrange(cln.shape[1]):
        if antennas[i].getStatus() != 33 or antennas[i].stand.id == ref:
            bln1[:,i] = 0.0
            bln2[:,i] = 0.0
            bln3[:,i] = 0.0
    
    # Create the FrameBuffer instance
    buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1], ReorderFrames=False)
    
    # Create the beam
    times = numpy.zeros(nChunks, dtype=numpy.float64)
    beam1 = numpy.zeros((nChunks, 2), dtype=numpy.float64)
    beam2 = numpy.zeros((nChunks, 2), dtype=numpy.float64)
    beam3 = numpy.zeros((nChunks, 2), dtype=numpy.float64)
    
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
        
        count = [0 for a in xrange(antpols)]
        
        j = 0
        fillsWork = framesWork / antpols
        # Inner loop that actually reads the frames into the data array
        while j < fillsWork:
            try:
                cFrame = tbn.readFrame(fh)
                k = k + 1
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
                continue
            
            for cFrame in cFrames:
                stand,pol = cFrame.header.parseID()
                
                # In the current configuration, stands start at 1 and go up to 260.  So, we
                # can use this little trick to populate the data array
                aStand = 2*(stand-1)+pol
                
                # Save the time
                if j == 0 and aStand == 0:
                    times[i] = cFrame.getTime()
                
                data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
                
                # Update the counters so that we can average properly later on
                count[aStand] = count[aStand] + 1
                
            j += 1
            
        # Mask
        bad = numpy.where( numpy.abs(data) >= 90 )
        data[bad] = 0.0
        
        # Beam forming
        taskPool = Pool(processes=6)
        taskList = []
        taskList.append( (i,1,0,taskPool.apply_async(formBeam, args=(data[0::2,:], bln1[0,0::2]))) )
        taskList.append( (i,1,1,taskPool.apply_async(formBeam, args=(data[1::2,:], bln1[0,1::2]))) )
        taskList.append( (i,2,0,taskPool.apply_async(formBeam, args=(data[0::2,:], bln2[0,0::2]))) )
        taskList.append( (i,2,1,taskPool.apply_async(formBeam, args=(data[1::2,:], bln2[0,1::2]))) )
        taskList.append( (i,3,0,taskPool.apply_async(formBeam, args=(data[0::2,:], bln3[0,0::2]))) )
        taskList.append( (i,3,1,taskPool.apply_async(formBeam, args=(data[1::2,:], bln3[0,1::2]))) )
        taskPool.close()
        taskPool.join()
        
        for i,b,p,task in taskList:
            if b == 1:
                beam1[i,p] = task.get()
            elif b == 2:
                beam2[i,p] = task.get()
            else:
                beam3[i,p] = task.get()
        
        print '1', beam1[i,0], '2', beam2[i,0], '3', beam3[i,0], '1/2', beam1[i,0]/beam2[i,0], '3/2', beam3[i,0]/beam2[i,0]
        del data
    
    # Plot the data
    print 'CygA      :', beam1[:,0]
    print 'Pointing 2:', beam2[:,0]
    print 'Pointing 1:', beam3[:,0]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    ax1.plot(times-times[0], beam1[:,0])
    ax1.plot(times-times[0], beam2[:,0])
    ax1.plot(times-times[0], beam3[:,0])
    ax2.plot(times-times[0], beam1[:,1])
    ax2.plot(times-times[0], beam2[:,1])
    ax2.plot(times-times[0], beam3[:,1])
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
    
