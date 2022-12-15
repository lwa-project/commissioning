#!/usr/bin/env python3

# Python2 compatibility
from __future__ import print_function, division
try:
    range = xrange
except NameError:
    pass
    
import os
import sys
import time
import ephem
import numpy
import argparse
from datetime import datetime, timedelta, tzinfo

from lsl import astro
from lsl.common.stations import parse_ssmif
from lsl.reader.ldp import TBWFile
from lsl.misc import beamformer
from lsl.common.dp import fS, SoftwareDP
import lsl.correlator.fx as fxc

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
    
    # Get an observer reader for calculations
    obs = station.get_observer()
    
    # Setup the beamformer gain and delay variables
    course = numpy.zeros(520)
    fine   = numpy.zeros(520)
    gains  = numpy.zeros((260,4))
    gains[:,0] = 1.0
    gains[:,3] = 1.0
    for ant in antennas:
        if ant.combined_status != 33:
            stand = (ant.digitizer - 1) / 2
            gains[stand,:] = 0.0
            
    # Setup the beamformer itself
    dp = SoftwareDP(mode='DRX', filter=7, central_freq=74e6)
    
    # Find the target azimuth/elevation to use
    idf = TBWFile(filenames[0])
    tStart = datetime.utcfromtimestamp(idf.get_info('start_time'))
    idf.close()
    
    obs.date = tStart.strftime("%Y/%m/%d %H:%M:%S")
    tTransit = obs.next_transit(args.source)
    obs.date = tTransit
    args.source.compute(obs)
    targetAz = args.source.az*180/numpy.pi
    targetEl = args.source.alt*180/numpy.pi
    
    # Preliminary report
    print("Working on %i TBW files using SSMIF '%s'" % (len(filenames), os.path.basename(ssmif)))
    print("  Source: '%s'" % args.source.name)
    print("    Transit time: %s" % str(tTransit))
    print("    Transit azimuth: %.2f degrees" % targetAz)
    print("    Transet elevation: %.2f degrees" % targetEl)
    print(" ")
    
    # Loop over input files
    unx, lst, pwrX, pwrY = [], [], [], []
    for filename in filenames:
        ## Get the file reader
        idf = TBWFile(filename)
        
        ## Pull out some metadata and update the observer
        jd = astro.unix_to_utcjd(idf.get_info('start_time'))
        obs.date = ephem.Date(jd - astro.DJD_OFFSET)
        sample_rate = idf.get_info('sample_rate')
        nInts = int(round( idf.get_info('nframe') / (30000.0 * len(antennas) / 2) ))
        transitOffset = (obs.date-tTransit)*86400.0
        
        ## Metadata report
        print("Filename: %s" % os.path.basename(filename))
        print("  Data type:  %s" % type(idf))
        print("  Captures in file: %i (%.3f s)" % (nInts, nInts*30000*400/sample_rate))
        print("  Station: %s" % station.name)
        print("  Date observed: %s" % str(obs.date))
        print("  MJD: %.5f" % (jd-astro.MJD_OFFSET,))
        print("  LST: %s" % str(obs.sidereal_time()))
        print("    %.1f s %s transit" % (abs(transitOffset), 'before' if transitOffset < 0 else 'after'))
        print(" ")
        
        ## Load in the data
        readT, t, data = idf.read(time_in_samples=True)
        
        ## Build up a time array
        t = t + numpy.arange(data.shape[1], dtype=numpy.int64)
        
        ## Update the beamformer delays for the pointing center(s)
        unx.append( idf.get_info('start_time') )
        lst.append( obs.sidereal_time() * 12/numpy.pi )
        pwrX.append( [] )
        pwrY.append( [] )
        
        for offset in (-1, 0, 1):
            ### Compute
            delays = beamformer.calc_delay(antennas, freq=74.0e6, azimuth=targetAz, elevation=targetEl+offset)
            delays *= fS*16
            delays = delays.max() - delays
            ### Decompose into FIFO and FIR
            course = (delays // 16)
            fine   = (delays % 16)
            
            ## Form the beams for both polarizations
            beamX, beamY = dp.form_beam(antennas, t, data, course, fine, gains)
            
            ## Compute the integrated spectra
            ### Convert to int16
            beam = numpy.zeros((2, beamX.size), dtype=numpy.int16)
            beam[0,:] = (numpy.round(beamX)).astype(data.dtype)
            beam[1,:] = (numpy.round(beamY)).astype(data.dtype)
            ### Move into the frequency domain
            freq, spec = fxc.SpecMaster(beam, LFFT=8192, window=fxc.null_window, verbose=False, sample_rate=fS, clip_level=0)
            
            ## Save
            pwrX[-1].append( spec[0,:] )
            pwrY[-1].append( spec[1,:] )
            
        ## Done
        idf.close()
        
    # Convert to arrays
    unx, lst = numpy.array(unx), numpy.array(lst)
    pwrX, pwrY = numpy.array(pwrX), numpy.array(pwrY)
    
    # Save for later (needed for debugging)
    outname = "estimateSEFD-%s-%04i%02i%02i.npz" % (os.path.splitext(os.path.basename(ssmif))[0], tTransit.tuple()[0], tTransit.tuple()[1], tTransit.tuple()[2])
    print("Saving intermediate data to '%s'" % outname)
    print(" ")
    numpy.savez(outname, source=args.source.name, freq=freq, 
                unx=unx, lst=lst, pwrX=pwrX, pwrY=pwrY)
                
    # Report
    print("%s" % (args.source.name,))
    for i in range(lst.size):
        print("%s:  %s  %s" % (str(ephem.hours(str(lst[i]))), pwrX[i,:], pwrY[i,:]))
        
    # Plot
    if args.plots:
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(lst, pwrX, linestyle='-', marker='+')
        ax.plot(lst, pwrY, linestyle='-', marker='x')
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="given an SSMIF and a collection of TBW files, use the SoftwareDP to form beams at the transit point of a source and estimate the system equivalent flux density (SEFD) and pointing error",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('ssmif', type=str,
                        help='station SSMIF')
    parser.add_argument('filename', type=str, nargs='+',
                        help='filename to process')
    parser.add_argument('-s', '--source', type=str, default='CygA',
                        help='source to use')
    parser.add_argument('-p', '--plots', action='store_true',
                        help='show summary plots at the end')
    args = parser.parse_args()
    main(args)
    
