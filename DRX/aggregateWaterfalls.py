#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a collection of NPZ files created by drxWatefall, average them down to a single 
spectrum and save them to a common NPZ file.

Usage:
  aggregateWaterfalls.py npz_file_1 npz_file_2 npz_file_3 ...

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import time
import numpy
from datetime import datetime


def main(args):
	# Gather together all of the files we should be working on
	filenames = args
	filenames.sort()
	
	print "Given %i NPZ files" % len(filenames)
	
	# Read in the first file to figure out what we are dealing with
	dataDict = numpy.load(filenames[0])
	freq = dataDict['freq']
	spec = dataDict['spec']
	try:
		srate = dataDict['srate']
	except KeyError:
		srate = 19.6e6
	standMapper = dataDict['standMapper']
	
	# Compute the actual integration time after processing
	tIntOriginal = dataDict['tInt']
	tIntActual = spec.shape[0] * dataDict['tInt']
	
	# Report
	print "FFT channels: %i" % freq.size
	print "Total time per file: %.3f s" % tIntActual
	
	# Loop over the files and average
	times = numpy.zeros(len(filenames))
	spec = numpy.zeros((len(filenames), spec.shape[1], freq.size))
	for i,filename in enumerate(filenames):
		dataDict = numpy.load(filename)
		
		cSpec = dataDict['spec']
		cTime = dataDict['times']
		
		times[i] = cTime[0]
		spec[i,:,:] = cSpec.mean(axis=0)
		
	# Compute the "fake" integration time that is needed by plotWaterfall to 
	# map points correctly on a click.
	tInt = times[1] - times[0]
	
	# Save
	outname = 'aggregated-waterfall.npz'
	numpy.savez(outname, freq=freq, times=times, spec=spec, tInt=tInt, tIntActual=tIntActual, tIntOriginal=tIntOriginal, srate=srate, standMapper=standMapper, filenames=filenames)


if __name__ == "__main__":
	main(sys.argv[1:])


