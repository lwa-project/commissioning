#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import time
import ephem
import numpy
import getopt
from datetime import datetime, timedelta, tzinfo

from lsl import astro
from lsl.common.stations import parseSSMIF
from lsl.reader.ldp import LWA1DataFile
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


def usage(exitCode=None):
    print """estimateSEFD.py - Given an SSMIF and a collection of TBW files, use
the SoftwareDP to form beams at the transit point of a source and estimate the
system equivalent flux density (SEFD) and pointing error.

Usage: estimateSEFD.py [OPTIONS] SSMIF tbw [tbw [...]]

Options:
-h, --help             Display this help information
-s, --source           Source to use (default = CygA)
-p, --plots            Show summary plots at the end (default = no)
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseConfig(args):
    config = {}
    # Command line flags - default values
    config['source'] = 'CygA'
    config['showPlots'] = False
    config['args'] = []

    # Read in and process the command line flags
    try:
        opts, arg = getopt.getopt(args, "hs:p", ["help", "source=", "plots"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-s', '--source'):
            config['source'] = value
        elif opt in ('-p', '--plots'):
            config['showPlots'] = True
        else:
            assert False
            
    # Add in arguments
    config['args'] = arg
    
    # Validate the inputs
    ## Find/validate the source
    src = None
    for line in _srcs:
        srcNew = ephem.readdb(line)
        if srcNew.name == config['source']:
            src = srcNew
            break
    if src is None:
        raise RuntimeError("Unknown source '%s'" % config['source'])
    else:
        config['source'] = src
    ## Argument length
    if len(config['args']) < 2:
        raise RuntimeError("Must provide both a SSMIF and stretch file")
        
    # Return configuration
    return config


def main(args):
    # Parse the command line
    config = parseConfig(args)
    
    # Break out the files we need
    ssmif = config['args'][0]
    filenames = config['args'][1:]
    
    # Setup the LWA station information
    station = parseSSMIF(ssmif)
    antennas = station.getAntennas()
    
    # Get an observer reader for calculations
    obs = station.getObserver()
    
    # Setup the beamformer gain and delay variables
    course = numpy.zeros(520)
    fine   = numpy.zeros(520)
    gains  = numpy.zeros((260,4))
    gains[:,0] = 1.0
    gains[:,3] = 1.0
    for ant in antennas:
        if ant.getStatus() != 33:
            stand = (ant.digitizer - 1) / 2
            gains[stand,:] = 0.0
            
    # Setup the beamformer itself
    dp = SoftwareDP(mode='DRX', filter=7, centralFreq=74e6)
    
    # Find the target azimuth/elevation to use
    idf = LWA1DataFile(filenames[0])
    tStart = datetime.utcfromtimestamp(idf.getInfo('tStart'))
    idf.close()
    
    obs.date = tStart.strftime("%Y/%m/%d %H:%M:%S")
    tTransit = obs.next_transit(config['source'])
    obs.date = tTransit
    config['source'].compute(obs)
    targetAz = config['source'].az*180/numpy.pi
    targetEl = config['source'].alt*180/numpy.pi
    
    # Preliminary report
    print "Working on %i TBW files using SSMIF '%s'" % (len(filenames), os.path.basename(ssmif))
    print "  Source: '%s'" % config['source'].name
    print "    Transit time: %s" % str(tTransit)
    print "    Transit azimuth: %.2f degrees" % targetAz
    print "    Transet elevation: %.2f degrees" % targetEl
    print " "
    
    # Loop over input files
    unx, lst, pwrX, pwrY = [], [], [], []
    for filename in filenames:
        ## Get the file reader
        idf = LWA1DataFile(filename)
        
        ## Pull out some metadata and update the observer
        jd = astro.unix_to_utcjd(idf.getInfo('tStart'))
        obs.date = ephem.Date(jd - astro.DJD_OFFSET)
        sampleRate = idf.getInfo('sampleRate')
        nInts = int(round( idf.getInfo('nFrames') / (30000.0 * len(antennas) / 2) ))
        transitOffset = (obs.date-tTransit)*86400.0
        
        ## Metadata report
        print "Filename: %s" % os.path.basename(filename)
        print "  Data type:  %s" % type(idf)
        print "  Captures in file: %i (%.3f s)" % (nInts, nInts*30000*400/sampleRate)
        print "  Station: %s" % station.name
        print "  Date observed: %s" % str(obs.date)
        print "  MJD: %.5f" % (jd-astro.MJD_OFFSET,)
        print "  LST: %s" % str(obs.sidereal_time())
        print "    %.1f s %s transit" % (abs(transitOffset), 'before' if transitOffset < 0 else 'after')
        print " "
        
        ## Load in the data
        readT, t, data = idf.read(timeInSamples=True)
        
        ## Build up a time array
        t = t + numpy.arange(data.shape[1], dtype=numpy.int64)
        
        ## Update the beamformer delays for the pointing center(s)
        unx.append( idf.getInfo('tStart') )
        lst.append( obs.sidereal_time() * 12/numpy.pi )
        pwrX.append( [] )
        pwrY.append( [] )
        
        for offset in (-1, 0, 1):
            ### Compute
            delays = beamformer.calcDelay(antennas, freq=74.0e6, azimuth=targetAz, elevation=targetEl+offset)
            delays *= fS*16
            delays = delays.max() - delays
            ### Decompose into FIFO and FIR
            course = (delays // 16)
            fine   = (delays % 16)
            
            ## Form the beams for both polarizations
            beamX, beamY = dp.formBeam(antennas, t, data, course, fine, gains)
            
            ## Compute the integrated spectra
            ### Convert to int16
            beam = numpy.zeros((2, beamX.size), dtype=numpy.int16)
            beam[0,:] = (numpy.round(beamX)).astype(data.dtype)
            beam[1,:] = (numpy.round(beamY)).astype(data.dtype)
            ### Move into the frequency domain
            freq, spec = fxc.SpecMaster(beam, LFFT=8192, window=fxc.noWindow, verbose=False, SampleRate=fS, ClipLevel=0)
            
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
    print "Saving intermediate data to '%s'" % outname
    print " "
    numpy.savez(outname, source=config['source'].name, freq=freq, 
                unx=unx, lst=lst, pwrX=pwrX, pwrY=pwrY)
                
    # Report
    print "%s" % (config['source'].name,)
    for i in xrange(lst.size):
        print "%s:  %s  %s" % (str(ephem.hours(str(lst[i]))), pwrX[i,:], pwrY[i,:])
        
    # Plot
    if config['showPlots']:
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(lst, pwrX, linestyle='-', marker='+')
        ax.plot(lst, pwrY, linestyle='-', marker='x')
        plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
    