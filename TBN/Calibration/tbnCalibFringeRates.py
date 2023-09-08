#!/usr/bin/env python3

import os
import sys
import ephem
import numpy
import argparse

from datetime import datetime

from lsl.astro import unix_to_utcjd, utcjd_to_unix
from lsl.common.stations import parse_ssmif, lwa1
from lsl.correlator.uvutils import compute_uvw
from lsl.statistics import robust
from lsl.common.progress import ProgressBar
from lsl.misc import parser as aph

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
    #
    # Gather the station meta-data from its various sources
    #
    if args.metadata is not None:
        site = parse_ssmif(args.metadata)
    else:
        site = lwa1
    observer = site.get_observer()
    antennas = site.antennas
    nAnts = len(antennas)
    
    #
    # Grab the time
    #
    if args.date is not None and args.time is not None:
        year, month, day = args.date.split('/', 2)
        year = int(year)
        month = int(month)
        day = int(day)
        
        hour, minute, second = args.time.split(':', 2)
        hour = int(hour)
        minute = int(minute)
        second = int(float(second))
        
        tNow = datetime(year, month, day, hour, minute, second)
        
    else:
        tNow = datetime.utcnow()
    observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S")
    print("Current time is %s" % tNow.strftime("%Y/%m/%d %H:%M:%S UTC"))
    print("Current LST at %s is %s" % (lwa1.name, observer.sidereal_time()))
    print(" ")
    
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
    print("Outriggers:")
    for outrigger in outriggers:
        print(" %s" % str(outrigger.stand))
    print(" ")
    print("Array Antenna:")
    print(" %s" % str(inside.stand))
    print(" ")
    
    #
    # Calculate the source positions to find what is up
    #
    print("Visible Sources:")
    for src in srcs:
        src.compute(observer)
        if src.alt > args.elevation_cut:
            print("  %s at %s degrees elevation" % (src.name, src.alt))
    print(" ")
    
    #
    # Calculate the fringe rates - At the specified time
    #
    allRates = {}
    for outrigger in outriggers:
        for src in srcs:
            src.compute(observer)
            if src.alt > args.elevation_cut:
                fRate = getFringeRate(inside, outrigger, observer, src, args.frequency)
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
                if src.alt > args.elevation_cut:
                    fRate = getFringeRate(antennas[0], outrigger, observer, src, args.frequency)
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
    print("%-4s  %-s" % ("Src", ''.join(["#%-10i  " % stand for stand in stands])))
    print("=" * (4+2+(11+2)*len(stands)))
    for src in srcNames:
        line = "%-4s  " % src
        for stand in stands:
            line = "%s%-+7.3f mHz  " % (line, allRates[src][stand]*1000.0)
        print(line)
        
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
    parser = argparse.ArgumentParser(
        description="calculate fringe rates for a few bright sources for baselines with all of the outriggers",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('date', type=aph.date, nargs='?',
                        help='UTC date to compute rates in YYYY/MM/DD')
    parser.add_argument('time', type=aph.time, nargs='?',
                        help='UTC time to compute rates in HH:MM:SS')
    parser.add_argument('-m', '--metadata', type=str,
                        help='name of SSMIF file to use for mappings')
    parser.add_argument('-f', '--frequency', type=aph.positive_float, default=74.0,
                        help='frequency in MHz')
    parser.add_argument('-e', '--elevation-cut', type=aph.positive_or_zero_float, default=10.0,
                        help='source elevation cut in degrees')
    args = parser.parse_args()
    args.frequency *= 1e6
    args.elevation_cut *= numpy.pi/180
    main(args)
    
