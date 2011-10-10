#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import ephem
import numpy
from datetime import datetime, timedelta

from scipy.signal import triang

from lsl.common.constants import c
from lsl.common.stations import lwa1
from lsl.correlator.uvUtils import computeUVW


# List of bright radio sources and pulsars in PyEphem format
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
	return numpy.dot(source, xyz) / c


def main(args):
	observer = lwa1.getObserver()
	antennas = lwa1.getAntennas()
	
	filename = args[0]
	pointingAz = float(args[1])
	pointingEl = float(args[2])
	dataDict = numpy.load(filename)
	
	ref  = dataDict['ref']
	refX = dataDict['refX'] - 1
	refY = dataDict['refY'] - 1
	centralFreq = dataDict['centralFreq']
	tInt = dataDict['tInt']
	times = dataDict['times']
	phase = dataDict['simpleVis']
	
	print "Central frequency: %.3f Hz" % centralFreq
	
	# Fix for the crappy file without times
	from datetime import datetime
	if times[0] == 0:
		firstTime = datetime(2011, 9, 22, 22, 42, 06)
		firstTime = int(firstTime.strftime("%s"))
		times = [firstTime+tInt*i for i in xrange(len(times))]
	
	beginDate = datetime.utcfromtimestamp(times[0])
	observer.date = beginDate.strftime("%Y/%m/%d %H:%M:%S")
	srcs = [ephem.Sun(),]
	for line in _srcs:
		srcs.append( ephem.readdb(line) )
	
	az = -99
	el = -99
	for i in xrange(len(srcs)):
		srcs[i].compute(observer)
		
		if srcs[i].alt > 0:
			print "source %s: alt %.1f degrees, az %.1f degrees" % (srcs[i].name, srcs[i].alt*180/numpy.pi, srcs[i].az*180/numpy.pi)
			
		if srcs[i].name == 'Sun':
			az = srcs[i].az  * 180.0/numpy.pi
			el = srcs[i].alt * 180.0/numpy.pi
			
	aln = []
	for i in xrange(phase.shape[1]):
		gd = getGeoDelay(antennas[i], az, el, Degrees=True)
		aln.append( numpy.exp(2j*numpy.pi*centralFreq*gd) )
	aln = numpy.array(aln)
	
	cln = numpy.zeros(phase.shape, dtype=numpy.complex128)
	for i in xrange(cln.shape[1]):
		if i % 2 == 0:
			cln[:,i] = phase[:,i] / phase[:,0]
		else:
			cln[:,i] = phase[:,i] / phase[:,1]
	cln /= aln
	delay = numpy.angle(cln) / (2*numpy.pi*centralFreq)
	
	alnPointing = []
	for i in xrange(phase.shape[1]):
		gd = getGeoDelay(antennas[i], pointingAz, pointingEl, Degrees=True)
		alnPointing.append( numpy.exp(2j*numpy.pi*centralFreq*gd) )
	alnPointing = numpy.array(alnPointing)
	blnPointing = (cln*alnPointing).conj() / numpy.abs(cln*alnPointing)
	
	delays = numpy.angle(blnPointing) / (2*numpy.pi*centralFreq)
	delays = delays.max() - delays
	
	import gain
	import delay
	dftBase = 'phased_beam_%.2faz,_%.2fel_%iMHz' % (pointingAz, poinintEl, centralFreq/1e6,)
	junk = delay.list2delayfile('.', dftBase, delays[0,:]*1e9)


if __name__ == "__main__":
	main(sys.argv[1:])
	
