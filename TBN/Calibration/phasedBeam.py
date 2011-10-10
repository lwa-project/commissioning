#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import ephem
import numpy
from datetime import datetime, timedelta

from scipy.signal import triang

from lsl.common.constants import c as vLight
from lsl.common.stations import lwa1
from lsl.correlator.uvUtils import computeUVW


# List of bright radio sources in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",]


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
	pointingAz = float(args[1])
	pointingEl = float(args[2])
	
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
	
	print "Central frequency: %.3f Hz" % centralFreq
	
	# Build the source list
	beginDate = datetime.utcfromtimestamp(times[0])
	observer.date = beginDate.strftime("%Y/%m/%d %H:%M:%S")
	srcs = [ephem.Sun(),]
	for line in _srcs:
		srcs.append( ephem.readdb(line) )
	
	# Identify the location of the reference source (the Sun in this case)
	az = -99
	el = -99
	for i in xrange(len(srcs)):
		srcs[i].compute(observer)
			
		if srcs[i].name == 'Sun':
			az = srcs[i].az  * 180.0/numpy.pi
			el = srcs[i].alt * 180.0/numpy.pi
	
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
	
	# Compute the geometric delay for the requested pointing
	alnPointing = []
	for i in xrange(phase.shape[1]):
		gd = getGeoDelay(antennas[i], pointingAz, pointingEl, Degrees=True)
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
	dftBase = 'phased_beam_%.2faz_%.2fel_%iMHz' % (pointingAz, pointingEl, centralFreq/1e6,)
	junk = delay.list2delayfile('.', dftBase, delays[0,:]*1e9)


if __name__ == "__main__":
	main(sys.argv[1:])
	
