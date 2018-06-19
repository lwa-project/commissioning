#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to generate delay-and-sum beam forming coefficients as well as a 
BAM script to move all beams with ~4 minute steps.

Usage:
trackSource <SSMIF> <source_name> <start date> <start time> <duration in hr>

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import pytz
import ephem
import numpy
from datetime import datetime, timedelta

from lsl.common import stations
from lsl.misc.beamformer import calcDelay


# Time zones
_UTC = pytz.utc
_MST = pytz.timezone('US/Mountain')


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",
         "B0329+54,f|J,03:32:59.37,+54:34:43.6,1",
         "B0809+74,f|J,08:14:59.44,+74:29:05.8,1", 
         "B0950+08,f|J,09:53:09.31,+07:55:35.8,1",
         "B1133+16,f|J,11:36:03.25,+15:51:04.5,1",
         "B1919+21,f|J,19:21:44.80,+21:53:01.8,1",
         "J2145-0750,f|J,21:45:50.47,-07:50:18.3,1",
         "J2339-0533,f|J,23:39:38.75,-05:33:05.3,1"]


# Source tracking step time in minutes
tStep = 4.0

# Observing frequency in Hz
centralFreq = 74.03e6

# Beams to use
beamsToUse = (1, 2, 4)


def main(args):
    # Divy up the command line arguments
    filename = args[0]
    source = args[1]
    startDate = args[2]
    startTime = args[3]
    duration  = float(args[4])
    
    year, month, day = startDate.split('/', 2)
    year = int(year)
    month = int(month)
    day = int(day)
    
    hour, minute, second = startTime.split(':', 2)
    hour = int(hour)
    minute = int(minute)
    second = int(second)
        
    tStart = _MST.localize(datetime(year, month, day, hour, minute, second))
    tStart = tStart.astimezone(_UTC)
    
    # Load the SSMIF
    station = stations.parseSSMIF(filename)

    # Gather the necessary information to figure out where things are
    observer = station.getObserver()
    antennas = station.getAntennas()

    # Find the "good" antennas to use
    digs    = numpy.array([ant.digitizer  for ant in antennas])
    ants    = numpy.array([ant.id         for ant in antennas])
    stands  = numpy.array([ant.stand.id   for ant in antennas])
    pols    = numpy.array([ant.pol        for ant in antennas])
    antStat = numpy.array([ant.status     for ant in antennas])
    feeStat = numpy.array([ant.fee.status for ant in antennas])

    badStands = numpy.where( antStat != 3 )[0]
    badFees   = numpy.where( feeStat != 3 )[0]
    bad = numpy.where( (stands > 256) | (antStat != 3) | (feeStat != 3) )[0]
    ## print "Number of bad stands:   %3i" % len(badStands)
    ## print "Number of bad FEEs:     %3i" % len(badFees)
    ## print "---------------------------"
    ## print "Total number bad inuts: %3i" % len(bad)
    ## print " "
    
    # Build the source list
    srcs = [ephem.Sun(), ephem.Jupiter(),]
    for line in _srcs:
        srcs.append( ephem.readdb(line) )
        
    # Identify the source to track
    refSource  = None
    for i in xrange(len(srcs)):
        if srcs[i].name.lower() == source.lower():
            refSource = srcs[i]
            source = refSource.name
    
    # Make sure we have a source to track
    if refSource is None:
        print "Unknown source '%s', quitting" % source
        sys.exit(1)
    
    print """#!/bin/bash
    
#
# Source tracking script for %s starting at %s
# -> tuning frequency is %.3f Hz
# -> track duration is %.3f hours
# -> update interval is %.3f minutes
#

""" % (source, tStart.astimezone(_MST), centralFreq, duration, tStep)
    
    # Create the DFT files and build the script
    nSteps = int(numpy.ceil(duration * 60 / 4))
    stepSize = timedelta(0, int(tStep*60), int((tStep*60*1000000) % 1000000))
    for s in xrange(nSteps):
        # Compute the source location half-way into the step
        tBeam = tStart + timedelta(0, int(tStep*60/2), int((tStep*60/2*1000000) % 1000000))
        observer.date = tBeam.strftime("%Y/%m/%d %H:%M:%S")
        refSource.compute(observer)
        
        pointingAz = refSource.az  * 180.0 / numpy.pi
        pointingEl = refSource.alt * 180.0 / numpy.pi
        
        # Compute the delays
        delays = calcDelay(antennas, freq=centralFreq, azimuth=pointingAz, elevation=pointingEl)
        delays *= 1e9
        delays = delays.max() - delays
        
        # Save - delays
        import delay
        dftBase = 'delay_beam_%s_%03i_%iMHz' % (source, (s+1), centralFreq/1e6,)
        junk = delay.list2delayfile('.', dftBase, delays)

        # Compute gains
        gains = [[1.0, 0.0, 0.0, 1.0]]*260 # initialize gain list
        for d in digs[bad]:
            # Digitizers start at 1, list indicies at 0
            i = d - 1
            gains[i/2] = [0,0,0,0]

        # Save - gains
        import gain
        gftBase = 'delay_beam_%s_%03i_%iMHz' % (source, (s+1), centralFreq/1e6,)
        junk = gain.list2gainfile('.', gftBase, gains)
        
        # Output script command - step start
        print """
#
# Begin step #%i at %s
# -> %s at %.3f az, %.3f el
#
tNow=`date -u +%%s `
tNow=$(($tNow*1))

## Wait for the right time
while [ $tNow -lt %s ]; do
    sleep 5
    tNow=`date -u +%%s `
    tNow=$(($tNow*1))
done

## Send BAM commands
tString=`date `
echo "Sending BAM commands for step #%i at $tString"
""" % ((s+1), tStart.astimezone(_MST), source, pointingAz, pointingEl, tStart.astimezone(_MST).strftime('%s'), (s+1))

        # Output script command - BAM commands
        for beam in beamsToUse:
            print """/home/joecraig/MCS/exec/mesix DP_ BAM "%i %s.df %s.gf 1"
sleep 1""" % (beam, dftBase, gftBase)

        # Output script command - step end
        print """
#
# End step #%i
#
""" % ((s+1),)
        
        # Update time
        tStart = tStart + stepSize


if __name__ == "__main__":
    main(sys.argv[1:])
    
