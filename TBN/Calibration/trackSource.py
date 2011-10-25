#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to generate phase-and-sum beam forming coefficients as well as a 
BAM script to move all beams with ~4 minute steps.

Usage:
trackSource <reference_NPZ> <source_name> <start date> <start time> <duration in hr>

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

from scipy.signal import triang

from lsl.common.constants import c as vLight
from lsl.common.stations import lwa1
from lsl.correlator.uvUtils import computeUVW


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
         "J2145-0750,f|J,21:45:50.47,-07:50:18.3,1",]


# Source tracking step time in minutes
tStep = 4.0


def getGeoDelay(antenna, az, el, Degrees=False):
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
	
	xyz = numpy.array([antenna.stand.x, antenna.stand.y, antenna.stand.z])
	return numpy.dot(source, xyz) / vLight


def main(args):
	# Gather the necessary information to figure out where things are
	observer = lwa1.getObserver()
	antennas = lwa1.getAntennas()
	
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
	
	# Load the data
	dataDict = numpy.load(filename)
	## Frequency
	centralFreq = dataDict['centralFreq']
	## Integration time
	tInt = dataDict['tInt']
	## Start times of the integrations
	times = dataDict['times']
	## The visiblity data
	phase = dataDict['simpleVis']
	
	# Build the source list
	beginDate = datetime.utcfromtimestamp(times[0])
	observer.date = beginDate.strftime("%Y/%m/%d %H:%M:%S")
	srcs = [ephem.Sun(),]
	for line in _srcs:
		srcs.append( ephem.readdb(line) )
		
	# Identify the location of the reference source (the Sun in this case) and the
	# reference source
	az = -99
	el = -99
	refSource  = None
	for i in xrange(len(srcs)):
		srcs[i].compute(observer)
			
		if srcs[i].name == 'Sun':
			az = srcs[i].az  * 180.0/numpy.pi
			el = srcs[i].alt * 180.0/numpy.pi
			
		if srcs[i].name.lower() == source.lower():
			refSource = srcs[i]
			source = refSource.name
	
	# Make sure we have a source to track
	if refSource is None:
		print "Unknown source '%s', quitting" % source
		sys.exit(1)
	
	# Generate geometric delay coefficients
	aln = []
	for i in xrange(phase.shape[1]):
		gd = getGeoDelay(antennas[i], az, el, Degrees=True)
		aln.append( numpy.exp(2j*numpy.pi*centralFreq*gd) )
	aln = numpy.array(aln)
	
	# Build the c^l_n values from Steve's "Fun with TBN" memo (Eqn. 10)
	cln = numpy.zeros(phase.shape, dtype=numpy.complex128)
	for i in xrange(cln.shape[1]):
		if i % 2 == 0:
			cln[:,i] = phase[:,i] / phase[:,0]
		else:
			cln[:,i] = phase[:,i] / phase[:,1]
	cln /= aln
	
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
		# Compute the source location
		observer.date = tStart.strftime("%Y/%m/%d %H:%M:%S")
		refSource.compute(observer)
		
		pointingAz = refSource.az
		pointingEl = refSource.alt
		
		# Compute the geometric delay for the requested pointing
		alnPointing = []
		for i in xrange(phase.shape[1]):
			gd = getGeoDelay(antennas[i], pointingAz, pointingEl, Degrees=False)
			alnPointing.append( numpy.exp(2j*numpy.pi*centralFreq*gd) )
		alnPointing = numpy.array(alnPointing)
	
		# Calculate the beamforming coefficients
		blnPointing = (cln*alnPointing).conj() / numpy.abs(cln*alnPointing)
		
		# Intepret these purely as delays
		delays = numpy.angle(blnPointing) / (2*numpy.pi*centralFreq)
		delays = delays.max() - delays
		
		# Save
		import gain
		import delay
		dftBase = 'phased_beam_%s_%03i_%iMHz' % (source, (s+1), centralFreq/1e6,)
		junk = delay.list2delayfile('.', dftBase, delays[0,:]*1e9)
		
		# Output script command
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
/home/joecraig/MCS/exec/mesix DP_ BAM "1 %s.df beams_56MHz_180az_77el_100bg.gf 1"
sleep 1
/home/joecraig/MCS/exec/mesix DP_ BAM "2 %s.df beams_56MHz_180az_77el_100bg.gf 1"
sleep 1
/home/joecraig/MCS/exec/mesix DP_ BAM "4 %s.df beams_56MHz_180az_77el_100bg.gf 1"
sleep 1
#
# End step #%i
#

""" % ((s+1), tStart.astimezone(_MST), source, pointingAz*180/numpy.pi, pointingEl*180/numpy.pi, tStart.astimezone(_MST).strftime('%s'), (s+1), dftBase, dftBase, dftBase, (s+1))
		
		# Update time
		tStart = tStart + stepSize


if __name__ == "__main__":
	main(sys.argv[1:])
	
