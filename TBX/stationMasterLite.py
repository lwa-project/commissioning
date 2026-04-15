#!/usr/bin/env python3

"""
Given a TBX file, plot the time averaged spectra for each digitizer input.  Save 
the data for later review with smGUI as an NPZ file.
"""

import os
import sys
import math
import numpy as np
import argparse

from lsl.common import stations
from lsl.reader import ldp, errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common.progress import ProgressBar
from lsl.common.data_access import DataAccess

import matplotlib.pyplot as plt


def main(args):
    # Set the station
    if args.metadata is not None:
        station = stations.parse_ssmif(args.metadata)
        ssmifContents = open(args.metadata).readlines()
    else:
        station = stations.lwa1
        with DataAccess.open('lwa1-ssmif.txt', 'r') as fh:
            ssmifContents = fh.readlines()
    antennas = station.antennas
    antpols = len(antennas)
    
    # Open the file and get its basic info
    idf = ldp.TBXFile(args.filename)
    nFrames = idf.get_info('nframe')
    nAnt = idf.get_info('nantenna')
    freq = idf.get_info('freq1')
    nchannels = len(freq)
    nFramesPerObs = nchannels // idf.get_info('frame_channel_count')
    beginDate = idf.get_info('start_time')
    nChunks = 1
    
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
        masterSpectra = np.zeros((nChunks, nAnt, nchannels), np.float32)
        
        tInt, tStart, data = idf.read(0.0)
        masterSpectra = np.abs(data)**2
        masterSpectra = masterSpectra.transpose(2,0,1)
        
        # Compute the 1 ms average power and the data range within each 1 ms window
        subSize = 1960
        nsegments = masterSpectra.shape[1] // subSize
        
        print("Computing average power and data range in %i-sample intervals, ADC histogram" % subSize)
        pb = ProgressBar(max=masterSpectra.shape[0])
        avgPower = np.zeros((antpols, nsegments), dtype=np.float32)
        dataRange = np.zeros((antpols, nsegments, 3), dtype=np.int16)
        adcHistogram = np.zeros((antpols, 4096), dtype=np.int32)
        histBins = range(-2048, 2049)
        
        # Apply the cable loss corrections, if requested
        if True:
            for s in range(masterSpectra.shape[1]):
                currGain = antennas[s].cable.gain(freq)
                for c in range(masterSpectra.shape[0]):
                    masterSpectra[c,s,:] /= currGain
                    
        # Now that we have read through all of the chunks, perform the final averaging by
        # dividing by all of the chunks
        spec = masterSpectra.mean(axis=0)
        
        # Estimate the dipole resonance frequencies
        print("Computing dipole resonance frequencies")
        pb = ProgressBar(max=spec.shape[0])
        resFreq = np.zeros(spec.shape[0])
        toCompare = np.where( (freq>31e6) & (freq<70e6) )[0]
        for i in range(spec.shape[0]):
            bestOrder = 0
            bestRMS = 1e34
            for j in range(3, 12):
                coeff = np.polyfit(freq[toCompare]/1e6, np.log10(spec[i,toCompare])*10, j)
                fit = np.polyval(coeff, freq[toCompare]/1e6)
                rms = ((fit - np.log10(spec[i,toCompare])*10)**2).sum()
                if rms < bestRMS:
                    bestOrder = j
                    bestRMS = rms
                    
            coeff = np.polyfit(freq[toCompare]/1e6, np.log10(spec[i,toCompare])*10, bestOrder)
            fit = np.polyval(coeff, freq[toCompare]/1e6)
            try:
                resFreq[i] = freq[toCompare[np.where( fit == fit.max() )[0][0]]] / 1e6
            except:
                pass
                
            pb.inc(amount=1)
            if pb.amount != 0 and pb.amount % 10 == 0:
                sys.stdout.write(pb.show()+'\r')
                sys.stdout.flush()
        sys.stdout.write(pb.show()+'\r')
        sys.stdout.write('\n')
        sys.stdout.flush()
        
        np.savez(outfile, date=str(beginDate), freq=freq, masterSpectra=masterSpectra, resFreq=resFreq, 
                    avgPower=avgPower, dataRange=dataRange, adcHistogram=adcHistogram, ssmifContents=ssmifContents)
    else:
        dataDict = np.load(outfile)
        freq = dataDict['freq']
        masterSpectra = dataDict['masterSpectra']
        resFreq = dataDict['resFreq']
        
        # Now that we have read through all of the chunks, perform the final averaging by
        # dividing by all of the chunks
        spec = masterSpectra.mean(axis=0)
    
    # Create a good template spectra
    specTemplate = np.median(spec, axis=0)
    specDiff = np.zeros(spec.shape[0])
    toCompare = np.where( (freq>32e6) & (freq<50e6) )[0]
    print(len(toCompare))
    for i in range(spec.shape[0]):
        specDiff[i] = (spec[i,toCompare] / specTemplate[toCompare]).mean()
    specDiff = np.where( specDiff < 2, specDiff, 2)
    
    # Get the station
    standPos = np.array([[ant.stand.x, ant.stand.y, ant.stand.z] for ant in antennas if ant.pol == 0])
    
    # Plots
    if args.verbose:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.scatter(standPos[:,0], standPos[:,1], c=specDiff[0::2], s=40.0, alpha=0.50)
        
        ## Set the limits to just zoom in on the main stations
        ax1.set_xlim([-65, 65])
        ax1.set_ylim([-65, 65])
        
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(freq/1e6, np.log10(specTemplate)*10, alpha=0.50)
        
        print("RBW: %.1f Hz" % (freq[1]-freq[0]))
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a TBX file and create a collection of time-averaged spectra', 
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
