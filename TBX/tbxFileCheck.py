#!/usr/bin/env python3

"""
Given a TBX file, check the time tags.
"""

import os
import sys
import math
import numpy as np
import argparse

from lsl.common.ndp import fS, fC
from lsl.common import stations
from lsl.reader import ldp, errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.misc import parser as aph


def main(args):
    # Set the station
    if args.lwasv:
        station = stations.lwasv
    elif args.lwana:
        station = stations.lwana
    else:
        station = stations.lwa1
    antennas = station.antennas
    antpols = len(antennas)
    
    # Open the file and get its basic info
    idf = ldp.TBXFile(args.filename)
    nFrames = idf.get_info('nframe')
    nAnt = idf.get_info('nantenna')
    freq = idf.get_info('freq1')
    nchannels = len(freq)
    nFramesPerObs = nchannels // idf.get_info('frame_channel_count')
    nsamp = nFrames // nFramesPerObs
    beginDate = idf.get_info('start_time')
    
    # File summary
    print(f"Filename: {args.filename}")
    print(f"Date of First Frame: {str(beginDate)}")
    print(f"Ant/Pols: {antpols}")
    print(f"Frequency Range: {freq[0]/1e6:.3f} to {freq[-1]/1e6:.3f} MHz")
    print(" ")
    
    chunkLength = int(round(args.length * fC)) / fC
    chunkSkip = int(round(args.skip * fC)) / fC
    
    # Output arrays
    clipFraction = []
    meanPower = []
    
    # Find what stands to report on
    nrpt = 0
    toUse = []
    toName = []
    polMap = {0: 'X', 1: 'Y'}
    for i,a in enumerate(antennas):
        if a.stand.id in (10,):
            nrpt += 1
            toUse.append(i)
            toName.append(f"{a.stand.id}{polMap[a.pol]}")
            
    i = 1
    print("   |     Clipping    |        Power      |")
    print("   |   10X     10Y   |    10X      10Y   |")
    print("---+-----------------+-------------------+")
    while True:
        try:
            tint, tstart, data = idf.read(chunkLength)
            
            clipFraction.append(np.zeros(nrpt))
            meanPower.append(np.zeros(nrpt))
            for j,k in enumerate(toUse):
                pwr = np.abs(data[k,:,:])**2
                pwr = pwr.astype(np.int32)
                
                bad = np.nonzero(pwr > args.trim_level)[0]
                clipFraction[-1][j] = len(bad)/pwr.shape[0]/pwr.shape[1]
                meanPower[-1][j] = pwr.mean()
                
            clip = clipFraction[-1]
            power = meanPower[-1]
            print("%2i | %6.2f%% %6.2f%% | %8.2f %8.2f |" % (i, clip[0]*100.0, clip[1]*100.0, power[0], power[1]))
            
            i += 1
            idf.offset(chunkSkip)
        except Exception as e:
            break
            
    clipFraction = np.array(clipFraction)
    meanPower = np.array(meanPower)
    
    clip = clipFraction.mean(axis=0)
    power = meanPower.mean(axis=0)
    
    print("---+-----------------+-------------------+")
    print("%2s | %6.2f%% %6.2f%% | %8.2f %8.2f |" % ('M', clip[0]*100.0, clip[1]*100.0, power[0], power[1]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='run through a TBX file and determine if it is bad or not', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-v', '--lwasv', action='store_true', 
                        help='use mapping from LWA-SV instead of LWA1')
    parser.add_argument('-n', '--lwana', action='store_true', 
                        help='use mapping from LWA-NA instead of LWA1')
    parser.add_argument('-l', '--length', type=aph.positive_float, default=1.0, 
                        help='length of time in seconds to analyze')
    parser.add_argument('-s', '--skip', type=aph.positive_float, default=900.0, 
                        help='skip period in seconds between chunks')
    parser.add_argument('-t', '--trim-level', type=aph.positive_float, default=49.0, 
                        help='trim level for power analysis with clipping')
    args = parser.parse_args()
    main(args)
