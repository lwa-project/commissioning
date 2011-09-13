#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simulate fringes for a dipole-dipole data set using the lsl.sim.vis.buildSimData()
function and the bright sources listed in lsl.sim.vis.srcs.

Usage:
  simulateFringesBright.py freq1 freq2 stand1 stand2 reference_file
where freq1 and freq1 are in MHz.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy

import aipy

from lsl.reader import drx
from lsl.common.dp import fS
from lsl.common import stations
from lsl.common.paths import data as dataPath
from lsl.astro import unix_to_utcjd
from lsl.sim import vis as simVis

from matplotlib import pyplot as plt

def main(args):
	cFreq1 = float(args[0])*1e6
	cFreq2 = float(args[1])*1e6
	stand1 = int(args[2])
	stand2 = int(args[3])
	filename = args[4]

	LFFT = 2048
	
	# Build up the station
	site = stations.lwa1
	
	rawAntennas = site.getAntennas()
	
	antennas = []
	for ant in rawAntennas:
		if ant.stand.id == stand1 and ant.pol == 0:
			antennas.append(ant)
	for ant in rawAntennas:
		if ant.stand.id == stand2 and ant.pol == 0:
			antennas.append(ant)
	# Get file parameters
	fh = open(filename, "rb")
	nFramesFile = os.path.getsize(filename) / drx.FrameSize
	junkFrame = drx.readFrame(fh)
	fh.close()

	tStart = junkFrame.getTime()
	srate = junkFrame.getSampleRate()

	# Get frequencies
	freq1 = numpy.fft.fftfreq(LFFT, d=1/srate)
	freq1 = cFreq1 + (numpy.fft.fftshift(freq1))[1:]
	freq2 = numpy.fft.fftfreq(LFFT, d=1/srate)
	freq2 = cFreq2 + (numpy.fft.fftshift(freq2))[1:]

	# Create the simulated array
	refJD = unix_to_utcjd(tStart)
	aa1 = simVis.buildSimArray(site, antennas, freq1/1e9, jd=refJD)
	aa2 = simVis.buildSimArray(site, antennas, freq2/1e9, jd=refJD)

	# Build the model times and range.  This is hard coded for now.
	tGap = 4.0
	jdList = []
	dTimes = []
	for i in xrange(899):
		tNow = tStart + i*tGap
		jdNow = unix_to_utcjd(tNow)
		print i, tNow, jdNow
		jdList.append(jdNow)
		dTimes.append( i*tGap )
		
	# Acutally run the simulations
	simDict1 = simVis.buildSimData(aa1, simVis.srcs, jd=jdList, pols=['xx',], verbose=True)
	simDict2 = simVis.buildSimData(aa2, simVis.srcs, jd=jdList, pols=['xx',], verbose=True)

	# Plot
	fig = plt.figure()
	ax1 = fig.add_subplot(2, 1, 1)
	ax2 = fig.add_subplot(2, 1, 2)

	vis1 = []
	for vis in simDict1['vis']['xx']:
		vis1.append( vis )
	vis2 = []
	for vis in simDict2['vis']['xx']:
		vis2.append( vis )
	
	vis1 = numpy.array(vis1)
	vis1 = numpy.ma.array(vis1, mask=~numpy.isfinite(vis1))
	vis2 = numpy.array(vis2)
	vis2 = numpy.ma.array(vis2, mask=~numpy.isfinite(vis2))

	data = numpy.abs(vis1)
	data = data.ravel()
	data.sort()
	vmax1 = data[int(round(0.90*len(data)))]
	
	data = numpy.abs(vis2)
	data = data.ravel()
	data.sort()
	vmax2 = data[int(round(0.90*len(data)))]

	ax1.imshow(numpy.abs(vis1), extent=(simDict1['freq'][0]/1e6, simDict1['freq'][-1]/1e6, dTimes[0], dTimes[-1]), 
			origin='lower', vmin=0.0, vmax=vmax1)
	ax2.imshow(numpy.abs(vis2), extent=(simDict2['freq'][0]/1e6, simDict2['freq'][-1]/1e6, dTimes[0], dTimes[-1]), 
			origin='lower', vmin=0.0, vmax=vmax2)

	ax1.axis('auto')
	ax2.axis('auto')

	ax1.set_xlabel('Frequency [MHz]')
	ax2.set_xlabel('Frequency [MHz]')
	ax1.set_ylabel('Elapsed Time [s]')
	ax2.set_ylabel('Elapsed Time [s]')

	plt.show()

if __name__ == '__main__':
	main(sys.argv[1:])

