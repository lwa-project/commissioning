#!/usr/bin/env python3

import os
import sys
import time
import ephem
import numpy as np
import argparse
from datetime import datetime, timedelta, tzinfo

from lsl import astro
from lsl.common.stations import parse_ssmif
from lsl.reader.ldp import TBXFile
from lsl.misc import beamformer
from lsl.misc import parser as aph
from lsl.common.progress import ProgressBarPlus

from matplotlib import pyplot as plt


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",
         "3C123,f|J,04:37:04.38,+29:40:13.8,1",
         "3C295,f|J,14:11:20.47,+52:12:09.5,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1",
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1"]


def main(args):
    # Break out the files we need
    ssmif = args.ssmif
    filenames = args.filename
    
    # Setup the LWA station information
    station = parse_ssmif(ssmif)
    antennas = station.antennas
    if args.cal_file is not None:
        ## Apply additional calibration as cable clock offsets
        caldata = np.loadtxt(args.cal_file)
        calstd = caldata[:,0]
        caldx = caldata[:,2]*1e-9
        caldy = caldata[:,4]*1e-9
        for a in antennas:
            for s,x,y in zip(calstd, caldx, caldy):
                if a.stand.id == s:
                    if a.pol == 0:
                        a.cable.clock_offset = x
                    else:
                        a.cable.clock_offset = y
                        
    # Sort out the integration time
    int_time = args.avg_time
    if int_time == 0:
        int_time = 5.0
        
    # Get an observer reader for calculations
    obs = station.get_observer()
    
    # Load in the sources and pull out the correct one
    sources = {}
    for line in _srcs:
        bdy = ephem.readdb(line)
        sources[bdy.name] = bdy
    try:
        args.source = sources[args.source]
    except KeyError:
        raise RuntimeError(f"Unknown target source '{args.source}'")
        
    # Find the target azimuth/altitude to use
    idf = TBXFile(filenames[0])
    tStart = idf.get_info('start_time').datetime
    idf.close()
    
    obs.date = tStart.strftime("%Y/%m/%d %H:%M:%S")
    tTransit = obs.next_transit(args.source)
    obs.date = tTransit
    args.source.compute(obs)
    targetAz = args.source.az*180/np.pi
    targetAlt = args.source.alt*180/np.pi
    
    # Preliminary report
    print(f"Working on {len(filenames)} TBX files using SSMIF '{os.path.basename(ssmif)}'")
    print(f"  Source: '{args.source.name}'")
    print(f"    Transit time: {tTransit}")
    print(f"    Transit azimuth: {targetAz:.2f} degrees")
    print(f"    Transit altitude: {targetAlt:.2f} degrees")
    print(" ")
    
    # Loop over input files
    unx, lst, pwrX, pwrY = [], [], [], []
    for filename in filenames:
        ## Get the file reader
        idf = TBXFile(filename)
        
        ## Pull out some metadata and update the observer
        jd = astro.unix_to_utcjd(idf.get_info('start_time'))
        obs.date = ephem.Date(jd - astro.DJD_OFFSET)
        args.source.compute(obs)
        sample_rate = idf.get_info('sample_rate')
        transitOffset = (obs.date-tTransit)*86400.0
        if transitOffset > 12*3600:
            transitOffset -= 86400.0
            
        # Get the frequencies in the file as well
        freqs = idf.get_info('freq1')
        
        ## Metadata report
        print(f"Filename: {os.path.basename(filename)}")
        print(f"  Data type:  {type(idf)}")
        print(f"  Station: {station.name}")
        print(f"  Date observed: {str(obs.date)}")
        print(f"  Frequency Range: {freqs[0]/1e6:.1f} to {freqs[-1]/1e6:.1f} MHz")
        print(f"  MJD: {jd-astro.MJD_OFFSET:.5f}")
        print(f"  LST: {str(obs.sidereal_time())}")
        print("    %.1f s %s transit" % (abs(transitOffset), 'before' if transitOffset < 0 else 'after'))
        print(f"  {args.source.name}: az {args.source.az}, el {args.source.alt}")
        print(" ")
        
        ## Load in the data
        readT, t, data = idf.read(int_time)
        
        ## Downselect to 74 +/- 2 MHz
        
        valid_freq = np.where((freqs > (args.freq-2e6)) & (freqs < (args.freq+2e6)))[0]
        if not len(valid_freq):
            raise RuntimeError(f"Cannot find data around {args.freq/1e6:.1f} MHz")
        freqs = freqs[valid_freq]
        data = data[:,valid_freq,:]
        
        ## Come up with the pattern
        pnts = []
        ### Scale for whether or not it is a mini-station
        pm_range = args.swing_range
        ### First, declination
        for offset in np.linspace(-pm_range, pm_range, args.nstep):
            pnts.append( (args.source._ra, ephem.degrees(args.source._dec+offset*np.pi/180)) )
        ### Now, RA
        for offset in np.linspace(-pm_range, pm_range, args.nstep):
            offset = offset / np.cos(args.source.dec)
            pnts.append( (ephem.hours(args.source._ra+offset*np.pi/180), args.source._dec) )
            
        unx.append( float(idf.get_info('start_time')) )
        lst.append( obs.sidereal_time() * 12/np.pi )
        pwrX.append( [] )
        pwrY.append( [] )
        pb = ProgressBarPlus(max=len(pnts))
        for pnt in pnts:
            bdy = ephem.FixedBody()
            bdy._ra = pnt[0]
            bdy._dec = pnt[1]
            bdy._epoch = ephem.J2000
            bdy.compute(obs)
            
            beam = beamformer.phase_and_sum(antennas, data,
                                            central_freq=freqs,
                                            azimuth=bdy.az*180/np.pi,
                                            altitude=bdy.alt*180/np.pi)
            beam = (np.abs(beam)**2).mean(axis=1).mean(axis=1)
            beamX, beamY = beam[0], beam[1]
            
            ## Save
            pwrX[-1].append( beamX )
            pwrY[-1].append( beamY )
            
            pb.inc()
            sys.stdout.write(pb.show()+'\r')
            sys.stdout.flush()
            
        ## Done
        idf.close()
        sys.stdout.write(pb.show()+'\r')
        sys.stdout.write('\n')
        sys.stdout.flush()
        
    # Convert to arrays
    unx, lst = np.array(unx), np.array(lst)
    pwrX, pwrY = np.array(pwrX), np.array(pwrY)
    
    # Save for later (needed for debugging)
    outname = "estimateSEFD-%s-%04i%02i%02i.npz" % (os.path.splitext(os.path.basename(ssmif))[0], tTransit.tuple()[0], tTransit.tuple()[1], tTransit.tuple()[2])
    print(f"Saving intermediate data to '{outname}'")
    print(" ")
    np.savez(outname, source=args.source.name, freq=freqs.mean(), 
             unx=unx, lst=lst, pwrX=pwrX, pwrY=pwrY, pnts=pnts)
                
    # Report
    print(f"{args.source.name}:")
    for i in range(lst.size):
        print(f"{ephem.hours(str(lst[i]))}:", pwrX[i,:], pwrY[i,:])
        
    # Plot
    if args.plots:
        fig = plt.figure()
        ax = fig.gca()
        for i in range(pwrX.shape[0]):
            ax.plot(np.arange(pwrX.shape[1]), pwrX[i,:], linestyle='-', marker='+', label=f"XX@{i+1}")
            ax.plot(np.arange(pwrY.shape[1]), pwrY[i,:], linestyle='-', marker='x', label=f"YY@{i+1}")
        ax.set_xlabel('Pointing #')
        ax.set_ylabel('Power [arb.]')
        ax.legend(loc=0)
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="given an SSMIF and a collection of TBX files, use phase-and-sum beamforming to make a basket weave pattern to estimate the system equivalent flux density (SEFD) and pointing error",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('ssmif', type=str,
                        help='station SSMIF')
    parser.add_argument('filename', type=str, nargs='+',
                        help='filename to process')
    parser.add_argument('-s', '--source', type=str, default='CygA',
                        help='source to use')
    parser.add_argument('-r', '--swing-range', type=aph.positive_float, default=12.0,
                        help='+/- swing of each weave track in degrees')
    parser.add_argument('-n', '--nstep', type=aph.positive_int, default=17,
                        help='number of steps per weave track')
    parser.add_argument('-t', '--avg-time', type=aph.positive_or_zero_float, default=0.0, 
                        help='integration time for the beam pointings; 0 = integrate the entire file')
    parser.add_argument('-f', '--freq', type=aph.frequency, default='74MHz',
                        help='frequency to estimate the SEFD at')
    parser.add_argument('-c', '--cal-file', type=str,
                        help='.sc calibration to apply on top of the SSMIF')
    parser.add_argument('-p', '--plots', action='store_true',
                        help='show summary plots at the end')
    args = parser.parse_args()
    main(args)
