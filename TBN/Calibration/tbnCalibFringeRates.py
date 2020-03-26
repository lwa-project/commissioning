#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import ephem
import numpy
import getopt

from datetime import datetime

from lsl.astro import unix_to_utcjd, utcjd_to_unix
from lsl.common.stations import parse_ssmif, lwa1
from lsl.correlator.uvutil import compute_uvw
from lsl.statistics import robust
from lsl.common.progress import ProgressBar

import lsl.sim.vis as simVis

from matplotlib import pyplot as plt

# List of bright radio sources and pulsars in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",]

def usage(exitCode=None):
    print """tbnCalibFringeRates.py - Calculate fringe rates for a few bright sources
for baselines with all of the outriggers.

Usage: tbnCalibFringeRates.py [OPTIONS] [YYYY/MM/DD HH:MM:SS.S]

Options:
-h, --help            Display this help information
-m, --metadata        Name of SSMIF file to use for mappings
-f, --frequency       Frequency in MHz (default = 74)
-e, --elevation-cut   Source elevation cut (default = 10 degrees)
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['SSMIF'] = None
    config['freq'] = 74e6
    config['elevCut'] = 10.0
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hm:f:e:", ["help", "metadata=", "frequency=", "elevation-cut="])
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
        elif opt in ('-f', '--frequency'):
            config['freq'] = float(value)*1e6
        elif opt in ('-e', '--elevation-cut'):
            config['elevCut'] = float(value)
        else:
            assert False
            
    # Convert the elevation cut to radians
    config['elevCut'] *= numpy.pi/180.0	
    
    # Add in arguments
    config['args'] = args
    
    # Return configuration
    return config


def getFringeRate(antenna1, antenna2, observer, src, freq):
    """
    Get the fringe rate for a baseline formed by antenna1-antenna2 for source src
    as observed by an observer.
    """
    
    # Update the source position
    src.compute(observer)
    
    # Calculate the hour angle
    HA = (float(observer.sidereal_time()) - float(src.ra))*12.0/numpy.pi
    dec = float(src.dec)*180/numpy.pi
    
    # Get the u,v,w coordinates
    uvw = compute_uvw([antenna1, antenna2], HA=HA, dec=dec, freq=freq)
    
    return -(2*numpy.pi/86164.0905)*uvw[0,0,0]*numpy.cos(src.dec)


def main(args):
    config = parseOptions(args)
    
    #
    # Gather the station meta-data from its various sources
    #
    if config['SSMIF'] is not None:
        site = parse_ssmif(config['SSMIF'])
    else:
        site = lwa1
    observer = site.get_observer()
    antennas = site.antennas
    nAnts = len(antennas)
    
    #
    # Grab the time
    #
    if len(config['args']) == 2:
        config['args'][0] = config['args'][0].replace('-', '/')
        year, month, day = config['args'][0].split('/', 2)
        year = int(year)
        month = int(month)
        day = int(day)
        
        hour, minute, second = config['args'][1].split(':', 2)
        hour = int(hour)
        minute = int(minute)
        second = int(second)
        
        tNow = datetime(year, month, day, hour, minute, second)
        
    else:
        tNow = datetime.utcnow()
    observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S")
    print "Current time is %s" % tNow.strftime("%Y/%m/%d %H:%M:%S UTC")
    print "Current LST at %s is %s" % (lwa1.name, observer.sidereal_time())
    print " "
    
    #
    # Load the sources source
    #
    srcs = [ephem.Sun(),]
    for line in _srcs:
        srcs.append( ephem.readdb(line) )
        
    #
    # Find the outriggers
    #
    outriggers = []
    inside = None
    best = 1e6
    for a in antennas:
        d = numpy.sqrt( a.stand.x**2 + a.stand.y**2 + a.stand.z**2 )
        if d < best:
            if a.pol == 0:
                best = d
                inside = a
        if d < 150:
            continue
            
        if a.pol == 0:
            outriggers.append(a)
    print "Outriggers:"
    for outrigger in outriggers:
        print " %s" % str(outrigger.stand)
    print " "
    print "Array Antenna:"
    print " %s" % str(inside.stand)
    print " "
    
    #
    # Calculate the source positions to find what is up
    #
    print "Visible Sources:"
    for src in srcs:
        src.compute(observer)
        if src.alt > config['elevCut']:
            print "  %s at %s degrees elevation" % (src.name, src.alt)
    print " "
    
    #
    # Calculate the fringe rates - At the specified time
    #
    allRates = {}
    for outrigger in outriggers:
        for src in srcs:
            src.compute(observer)
            if src.alt > config['elevCut']:
                fRate = getFringeRate(inside, outrigger, observer, src, config['freq'])
                try:
                    allRates[src.name][outrigger.stand.id] = fRate
                except KeyError:
                    allRates[src.name] = {outrigger.stand.id: fRate}
                    
    #
    # Calculate the fringe rates - Two hours centered on the specified time
    #
    allRates2 = {}
    for outrigger in outriggers:
        observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S")
        observer.date = observer.date - 1.0/24.0
        
        for tRel in numpy.linspace(-1.0, 1.0, 61):
            observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S")
            observer.date = observer.date + tRel/24.0
            
            for src in srcs:
                src.compute(observer)
                if src.alt > config['elevCut']:
                    fRate = getFringeRate(antennas[0], outrigger, observer, src, config['freq'])
                    try:
                        allRates2[src.name][outrigger.stand.id].append(fRate)
                    except KeyError:
                        try:
                            allRates2[src.name][outrigger.stand.id] = [fRate,]
                        except KeyError:
                            allRates2[src.name] = {}
                            allRates2[src.name][outrigger.stand.id] = [fRate,]
                            
    #
    # Report
    #
    srcNames = allRates.keys()
    stands = [outrigger.stand.id for outrigger in outriggers]
    print "%-4s  %-s" % ("Src", ''.join(["#%-10i  " % stand for stand in stands]))
    print "=" * (4+2+(11+2)*len(stands))
    for src in srcNames:
        line = "%-4s  " % src
        for stand in stands:
            line = "%s%-+7.3f mHz  " % (line, allRates[src][stand]*1000.0)
        print line
        
    #
    # Plot
    #
    fig = plt.figure()
    axs = []
    for i,stand in enumerate(stands):
        axs.append( fig.add_subplot(2, 3, i+1) )
        axs[-1].set_title('Stand #%i - Stand #%i' % (stand, inside.stand.id))
        
    for src in srcNames:
        if src not in ('CygA', 'CasA', 'Sun', 'VirA', 'TauA'):
            continue
        for i,stand in enumerate(stands):
            data = numpy.array(allRates2[src][stand])
            axs[i].hist(data*1000.0, bins=5, label=src)
            
    for i,stand in enumerate(stands):
        fRate = 0.0
        ylim = axs[i].get_ylim()
        axs[i].plot([fRate, fRate], ylim, label='DC')

        axs[i].set_xlim((-10, 10))
        axs[i].set_xlabel('Rate [mHz]')
        axs[i].set_ylim(ylim)
        axs[i].set_ylabel('Number of Minutes')
        axs[i].legend(loc=0)
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
    