#!/usr/bin/env python

"""
Given a TBF file, plot the time averaged spectra for each digitizer input.  Save 
the data for later review with smGUI as an NPZ file.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import math
import numpy
import ephem
import argparse

from lsl.common import stations
from lsl.reader import tbf
from lsl.reader import errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common.progress import ProgressBar
from lsl.common.paths import DATA as dataPath

import matplotlib.pyplot as plt


def main(args):
    # Set the station
    if args.metadata is not None:
        station = stations.parse_ssmif(args.metadata)
        ssmifContents = open(args.metadata).readlines()
    else:
        station = stations.lwasv
        ssmifContents = open(os.path.join(dataPath, 'lwa1-ssmif.txt')).readlines()
    antennas = station.antennas
    antpols = len(antennas)
    
    fh = open(args.filename, "rb")
    nFrames = os.path.getsize(args.filename) / tbf.FRAME_SIZE
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbf.read_frame(fh)
    fh.seek(0)
    beginTime = junkFrame.get_time()
    beginDate = ephem.Date(unix_to_utcjd(junkFrame.get_time()) - DJD_OFFSET)
    
    # Figure out how many frames there are per observation and the number of
    # channels that are in the file
    nFramesPerObs = tbf.get_frames_per_obs(fh)
    nchannels = tbf.get_channel_count(fh)
    
    # Figure out how many chunks we need to work with
    nChunks = nFrames / nFramesPerObs
    
    # Pre-load the channel mapper
    mapper = []
    freq = []
    for i in xrange(2*nFramesPerObs):
        cFrame = tbf.read_frame(fh)
        mapper.append( cFrame.header.first_chan )
        freq.extend( list(cFrame.header.channel_freqs) )
    fh.seek(-2*nFramesPerObs*tbf.FRAME_SIZE, 1)
    mapper.sort()
    freq.sort()
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Frames per Observation: %i" % nFramesPerObs)
    print("Channel Count: %i" % nchannels)
    print("Frames: %i" % nFrames)
    print("===")
    print("Chunks: %i" % nChunks)
    
    outfile = os.path.split(args.filename)[1]
    outfile = os.path.splitext(outfile)[0]
    outfile = "%s.npz" % outfile	
    if (not os.path.exists(outfile)) or args.force:
        # Master loop over all of the file chunks
        masterSpectra = numpy.zeros((nChunks, 512, nchannels), numpy.float32)
        
        for i in range(nChunks):
            # Inner loop that actually reads the frames into the data array
            for j in range(nFramesPerObs):
                # Read in the next frame and anticipate any problems that could occur
                try:
                    cFrame = tbf.read_frame(fh)
                except errors.EOFError:
                    break
                except errors.SyncError:
                    print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbf.FRAME_SIZE-1))
                    continue
                if not cFrame.header.is_tbf:
                    continue
                    
                aStand = mapper.index(cFrame.header.first_chan)
                
                # In the current configuration, stands start at 1 and go up to 10.  So, we
                # can use this little trick to populate the data array
                if cFrame.header.frame_count % 10000 == 0 and args.verbose:
                    print("%4i -> %3i  %6.3f  %5i  %i" % (cFrame.header.first_chan, aStand, cFrame.get_time(), cFrame.header.frame_count, cFrame.data.timetag))
                    
                # Actually load the data.  x pol goes into the even numbers, y pol into the 
                # odd numbers
                if i == 0 and j == 0:
                    refCount = cFrame.header.frame_count
                count = cFrame.header.frame_count - refCount
                masterSpectra[count,0::2,aStand*12:(aStand+1)*12] = numpy.abs( numpy.rollaxis(cFrame.payload.data[:,:,0], 1) )**2
                masterSpectra[count,1::2,aStand*12:(aStand+1)*12] = numpy.abs( numpy.rollaxis(cFrame.payload.data[:,:,1], 1) )**2

                
            # Compute the 1 ms average power and the data range within each 1 ms window
            subSize = 1960
            nsegments = masterSpectra.shape[1] / subSize
            
            print("Computing average power and data range in %i-sample intervals, ADC histogram" % subSize)
            pb = ProgressBar(max=masterSpectra.shape[0])
            avgPower = numpy.zeros((antpols, nsegments), dtype=numpy.float32)
            dataRange = numpy.zeros((antpols, nsegments, 3), dtype=numpy.int16)
            adcHistogram = numpy.zeros((antpols, 4096), dtype=numpy.int32)
            histBins = range(-2048, 2049)
            
        # Apply the cable loss corrections, if requested
        if True:
            for s in xrange(masterSpectra.shape[1]):
                currGain = antennas[s].cable.gain(freq)
                for c in xrange(masterSpectra.shape[0]):
                    masterSpectra[c,s,:] /= currGain
                    
        # Now that we have read through all of the chunks, perform the final averaging by
        # dividing by all of the chunks
        spec = masterSpectra.mean(axis=0)
        
        # Estimate the dipole resonance frequencies
        print("Computing dipole resonance frequencies")
        pb = ProgressBar(max=spec.shape[0])
        resFreq = numpy.zeros(spec.shape[0])
        toCompare = numpy.where( (freq>31e6) & (freq<70e6) )[0]
        for i in xrange(spec.shape[0]):
            bestOrder = 0
            bestRMS = 1e34
            for j in xrange(3, 12):
                coeff = numpy.polyfit(freq[toCompare]/1e6, numpy.log10(spec[i,toCompare])*10, j)
                fit = numpy.polyval(coeff, freq[toCompare]/1e6)
                rms = ((fit - numpy.log10(spec[i,toCompare])*10)**2).sum()
                if rms < bestRMS:
                    bestOrder = j
                    bestRMS = rms
                    
            coeff = numpy.polyfit(freq[toCompare]/1e6, numpy.log10(spec[i,toCompare])*10, bestOrder)
            fit = numpy.polyval(coeff, freq[toCompare]/1e6)
            try:
                resFreq[i] = freq[toCompare[numpy.where( fit == fit.max() )[0][0]]] / 1e6
            except:
                pass
                
            pb.inc(amount=1)
            if pb.amount != 0 and pb.amount % 10 == 0:
                sys.stdout.write(pb.show()+'\r')
                sys.stdout.flush()
        sys.stdout.write(pb.show()+'\r')
        sys.stdout.write('\n')
        sys.stdout.flush()
        
        numpy.savez(outfile, date=str(beginDate), freq=freq, masterSpectra=masterSpectra, resFreq=resFreq, 
                    avgPower=avgPower, dataRange=dataRange, adcHistogram=adcHistogram, ssmifContents=ssmifContents)
    else:
        dataDict = numpy.load(outfile)
        freq = dataDict['freq']
        masterSpectra = dataDict['masterSpectra']
        resFreq = dataDict['resFreq']
        
        # Now that we have read through all of the chunks, perform the final averaging by
        # dividing by all of the chunks
        spec = masterSpectra.mean(axis=0)
    
    # Create a good template spectra
    specTemplate = numpy.median(spec, axis=0)
    specDiff = numpy.zeros(spec.shape[0])
    toCompare = numpy.where( (freq>32e6) & (freq<50e6) )[0]
    print(len(toCompare))
    for i in xrange(spec.shape[0]):
        specDiff[i] = (spec[i,toCompare] / specTemplate[toCompare]).mean()
    specDiff = numpy.where( specDiff < 2, specDiff, 2)
    
    # Get the station
    standPos = numpy.array([[ant.stand.x, ant.stand.y, ant.stand.z] for ant in antennas if ant.pol == 0])
    
    # Plots
    if args.verbose:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.scatter(standPos[:,0], standPos[:,1], c=specDiff[0::2], s=40.0, alpha=0.50)
        
        ## Set the limits to just zoom in on the main stations
        ax1.set_xlim([-65, 65])
        ax1.set_ylim([-65, 65])		
        
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(freq/1e6, numpy.log10(specTemplate)*10, alpha=0.50)
        
        print("RBW: %.1f Hz" % (freq[1]-freq[0]))
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a TBF file and create a collection of time-averaged spectra', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to convert')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of SSMIF file to use for mappings')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false', 
                        help='run %(prog)s in silent mode')
    args = parser.parse_args()
    main(args)
    
