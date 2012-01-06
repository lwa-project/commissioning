#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to fringe special DRX files that have a beam X pol. and a dipole 
on Y pol.  The visibilites are written to a NPZ file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy

from lsl.statistics import robust

from lsl.reader import drx
from lsl.common.dp import fS
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.correlator import fx as fxc

from matplotlib import pyplot as plt

def main(args):
	stand1 = 0
	stand2 = int(args[0])
	filename = args[1]
	
	# Build up the station
	site = stations.lwa1
	
	rawAntennas = site.getAntennas()

	antennas = []
	xyz = numpy.zeros((len(rawAntennas),3))
	i = 0
	for ant in rawAntennas:
		xyz[i,0] = ant.stand.x
		xyz[i,1] = ant.stand.y
		xyz[i,2] = ant.stand.z
		i += 1

	arrayX = xyz[:,0].mean()
	arrayY = xyz[:,1].mean()
	arrayZ = xyz[:,2].mean()

	beamStand   = stations.Stand(0, arrayX, arrayY, arrayZ)
	beamFEE     = stations.FEE('Beam', 0, gain1=0, gain2=0, status=3)
	beamCable   = stations.Cable('Beam', 0, vf=1.0)
	beamAntenna = stations.Antenna(0, stand=beamStand, pol=0, theta=0, phi=0, status=3)
	beamAntenna.fee = beamFEE
	beamAntenna.feePort = 1
	beamAntenna.cable = beamCable

	antennas.append(beamAntenna)
	
	for ant in rawAntennas:
		if ant.stand.id == stand2 and ant.pol == 0:
			antennas.append(ant)
	
	fh = open(filename, "rb")
	nFramesFile = os.path.getsize(filename) / drx.FrameSize
	junkFrame = drx.readFrame(fh)
	fh.seek(0)
	
	beam, tune, pol = junkFrame.parseID()
	srate = junkFrame.getSampleRate()
	
	tunepols = drx.getFramesPerObs(fh)
	tunepols = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepols
		
	tnom = junkFrame.header.timeOffset
	tStart = junkFrame.getTime()
	
	# Get the DRX frequencies
	cFreq1 = 0.0
	cFreq2 = 0.0
	for i in xrange(4):
		junkFrame = drx.readFrame(fh)
		b,t,p = junkFrame.parseID()
		if p == 0 and t == 1:
			try:
				cFreq1 = junkFrame.getCentralFreq()
			except AttributeError:
				from lsl.common.dp import fS
				cFreq1 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
		elif p == 0 and t == 2:
			try:
				cFreq2 = junkFrame.getCentralFreq()
			except AttributeError:
				from lsl.common.dp import fS
				cFreq2 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
		else:
			pass
	fh.seek(-4*drx.FrameSize, 1)
	
	# Align the files as close as possible by the time tags and then make sure that
	# the first frame processed is from tuning 1, pol 0.
	junkFrame = drx.readFrame(fh)
	beam, tune, pol = junkFrame.parseID()
	pair = 2*(tune-1) + pol
	j = 0
	while pair != 0:
		junkFrame = drx.readFrame(fh)
		beam, tune, pol = junkFrame.parseID()
		pair = 2*(tune-1) + pol
		j += 1
	fh.seek(-drx.FrameSize, 1)
	print "Shifted beam %i data by %i frames (%.4f s)" % (beam, j, j*4096/srate/4)
	
	# Set integration time
	tInt = 4.0
	nFrames = int(round(tInt*srate/4096))
	tInt = nFrames*4096/srate
	
	# Read in some data
	tFile = nFramesFile / 4 * 4096 / srate
	
	filePos = fh.tell()
		
	# Read in the first 100 frames for each tuning/polarization
	count = {0:0, 1:0, 2:0, 3:0}
	data = numpy.zeros((4, 4096*100), dtype=numpy.csingle)
	for i in xrange(4*100):
		try:
			cFrame = drx.readFrame(fh, Verbose=False)
		except errors.eofError:
			break
		except errors.syncError:
			continue
		
		beam,tune,pol = cFrame.parseID()
		aStand = 2*(tune-1) + pol
		
		data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
		count[aStand] +=  1
	
	# Go back to where we started
	fh.seek(filePos)
	
	# Compute the robust mean and standard deviation for I and Q for each
	# tuning/polarization
	meanI = []
	meanQ = []
	stdsI = []
	stdsQ = []
	for i in xrange(4):
		meanI.append( robust.mean(data[i,:].real) )
		meanQ.append( robust.mean(data[i,:].imag) )
		
		stdsI.append( robust.std(data[i,:].real) )
		stdsQ.append( robust.std(data[i,:].imag) )
	
	# Report
	print "Statistics:"
	for i in xrange(4):
		print " Mean %i: %.3f + %.3f j" % (i+1, meanI[i], meanQ[i])
		print " Std  %i: %.3f + %.3f j" % (i+1, stdsI[i], stdsQ[i])
	
	# Come up with the clip levels based on 4 sigma
	clip1 = (meanI[0] + meanI[1] + meanQ[0] + meanQ[1]) / 4.0
	clip2 = (meanI[2] + meanI[3] + meanQ[2] + meanQ[3]) / 4.0
	
	clip1 += 5*(stdsI[0] + stdsI[1] + stdsQ[0] + stdsQ[1]) / 4.0
	clip2 += 5*(stdsI[2] + stdsI[3] + stdsQ[2] + stdsQ[3]) / 4.0
	
	clip1 = int(round(clip1))
	clip2 = int(round(clip2))
	
	# Report again
	print "Clip Levels:"
	print " Tuning 1: %i" % clip1
	print " Tuning 2: %i" % clip2
	
	print filename, tFile, tInt, int(tFile/tInt)
	nChunks = int(tFile/tInt)
	for i in xrange(nChunks):
		junkFrame = drx.readFrame(fh)
		tStart = junkFrame.getTime()
		fh.seek(-drx.FrameSize, 1)

		count1 = [0,0]
		data1 = numpy.zeros((2, 4096*nFrames), dtype=numpy.complex64)
		count2 = [0,0]
		data2 = numpy.zeros((2, 4096*nFrames), dtype=numpy.complex64)
		for j in xrange(nFrames):
			for k in xrange(4):
				cFrame = drx.readFrame(fh)
				beam, tune, pol = cFrame.parseID()
				pair = 2*(tune-1) + pol

				if tune == 1:
					data1[pol, count1[pol]*4096:(count1[pol]+1)*4096] = cFrame.data.iq
					count1[pol] += 1
				else:
					data2[pol, count2[pol]*4096:(count2[pol]+1)*4096] = cFrame.data.iq
					count2[pol] += 1
					
		# Correlate
		LFFT = 2048
		for ant in antennas:
			print i+1, tStart, tInt, str(ant)
	
		blList1, freq1, vis1 = fxc.FXMaster(data1, antennas, LFFT=LFFT, Overlap=1, IncludeAuto=True, verbose=False, SampleRate=srate, CentralFreq=cFreq1, Pol='XX', ReturnBaselines=True, GainCorrect=False, ClipLevel=clip1)
		
		blList2, freq2, vis2 = fxc.FXMaster(data2, antennas, LFFT=LFFT, Overlap=1, IncludeAuto=True, verbose=False, SampleRate=srate, CentralFreq=cFreq2, Pol='XX', ReturnBaselines=True, GainCorrect=False, ClipLevel=clip2)
	
		if nChunks != 1:
			outfile = filename.replace('.dat', '-vis-%04i.npz' % (i+1))
		else:
			outfile = filename.replace('.dat', '-vis.npz')
		numpy.savez(outfile, freq1=freq1, vis1=vis1, freq2=freq2, vis2=vis2, tStart=tStart, tInt=tInt, stands=numpy.array([stand1, stand2]))

		del data1
		del data2

	# Plot
	fig = plt.figure()
	i = 0
	for bl, vi in zip(blList1, vis1):
		ax = fig.add_subplot(4, 3, i+1)
		ax.plot(freq1/1e6, numpy.unwrap(numpy.angle(vi)))
		ax.set_title('Stand %i - Stand %i' % (bl[0].stand.id, bl[1].stand.id))
		ax = fig.add_subplot(4, 3, i+4)
		ax.plot(freq1/1e6, numpy.abs(vi))
		i += 1

		coeff = numpy.polyfit(freq1, numpy.unwrap(numpy.angle(vi)), 1)
		print coeff[0]/2/numpy.pi*1e9, coeff[1]*180/numpy.pi
		
	i = 6
	for bl, vi in zip(blList2, vis2):
		ax = fig.add_subplot(4, 3, i+1)
		ax.plot(freq2/1e6, numpy.unwrap(numpy.angle(vi)))
		ax.set_title('Stand %i - Stand %i' % (bl[0].stand.id, bl[1].stand.id))
		ax = fig.add_subplot(4, 3, i+4)
		ax.plot(freq2/1e6, numpy.abs(vi))
		i += 1

		coeff = numpy.polyfit(freq2, numpy.unwrap(numpy.angle(vi)), 1)
		print coeff[0]/2/numpy.pi*1e9, coeff[1]*180/numpy.pi

	#plt.show()
	

if __name__ == "__main__":
	main(sys.argv[1:])
	
