#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script takes a DROSv2  spectrometer data and checks the 
power and clip levels.

Usage:
  specFileCheck.py spectrometer_file

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import numpy
import getopt
from datetime import datetime

from lsl.reader import drspec
from lsl.common import progress
from optparse import OptionParser

if __name__ == '__main__':
	parser = OptionParser(usage="specFileCheck.py spectrometer_file", description="Inspect the data in a DR spectrometer file")
	(options, args) = parser.parse_args()
	
	fh = open(args[0], 'rb')
	# Interogate the file to figure out what frames sizes to expect, now many 
	# frames there are, and what the transform length is
	nFrames = os.path.getsize(args[0]) / drspec.getFrameSize(fh)
	nChunks = nFrames
	LFFT = drspec.getTransformSize(fh)
	
	# Read in the first frame to figure out the DP information
	cPos = fh.tell()
	junkFrame = drspec.readFrame(fh)
	fh.seek(cPos)
	
	beam = junkFrame.parseID()
	centralFreq1 = junkFrame.getCentralFreq(1)
	centralFreq2 = junkFrame.getCentralFreq(2)
	srate = junkFrame.getSampleRate()
	tInt = junkFrame.header.nInts*LFFT/srate
	decFactor = junkFrame.header.decimation
	beginDate = datetime.utcfromtimestamp(junkFrame.getTime())
	
	products = junkFrame.getDataProducts()
	nProducts = len(products)
	
	# Report
	print "Filename: %s" % args[0]
	print "Date of First Frame: %s" % beginDate
	print "Beam: %i" % beam
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Data Products: %s" % (', '.join(products),)
	print "Frames: %i (%.3f s)" % (nFrames, nFrames*tInt)
	print "---"
	print "Transform Length: %i" % LFFT
	print "Integration: %.3f s" % tInt
	print "Decimation Factor: %.3f" % decFactor
	
	percentages = numpy.zeros((4,nFrames),numpy.dtype('float'))
	power = numpy.zeros((2*nProducts,nFrames),numpy.dtype('float'))
	
	try:
		pbar = progress.ProgressBarPlus(max=nFrames)
	except:
		pbar = progress.ProgressBar(max=nFrames)
		
	for i in xrange(nFrames):
		frame = drspec.readFrame(fh)
		
		data = []
		for j,tuning in enumerate((0, 1)):
			for k,product in enumerate(products):
				data.append( getattr(frame.data, '%s%i' % (product, tuning)).mean() )
		for j in range(nProducts*2):
			power[j][i] = data[j]
		for j in range(4):
			percentages[j][i] = int(frame.data.saturations[j])/(srate*tInt)
			
		pbar.inc(1)
		if i % 10 == 0:
			sys.stdout.write(pbar.show()+'\r')
			sys.stdout.flush()
			
	sys.stdout.write(pbar.show()+'\n')
	sys.stdout.flush()
	
	print "Average Saturation Counts:"
	print "Tuning 1 XX: %6.2f%% +/- %6.2f%%" % (percentages[0].mean()*100,  percentages[0].std()*100)
	print "Tuning 1 YY: %6.2f%% +/- %6.2f%%" % (percentages[1].mean()*100,  percentages[1].std()*100)
	print "Tuning 2 XX: %6.2f%% +/- %6.2f%%" % (percentages[2].mean()*100,  percentages[2].std()*100)
	print "Tuning 2 YY: %6.2f%% +/- %6.2f%%" % (percentages[3].mean()*100,  percentages[3].std()*100)
	print "---"
	
	print "Average Power Levels:"
	for j,tuning in enumerate((0, 1)):
		for k,product in enumerate(products):
			print "Tuning %i %s: %5.2f" % (tuning+1, product, power[j*nProducts+k].mean())
			
	fh.close()