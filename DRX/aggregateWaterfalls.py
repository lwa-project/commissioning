#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a collection of NPZ files created by drxWatefall, average them down to a single 
spectrum and save them to a common NPZ file.  The output filename is always 
'aggregated-waterfall.npz' and the file is saved in the current directory.

Usage:
  aggregateWaterfalls.py [OPTIONS] npz_file_1 npz_file_2 npz_file_3 ...

Options:
  -h, --help       Show this help message.
  -m, --maximum    Save the spectra associated with the maximum power in
                   each beam/tuning/polarization.

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
	# Figure out if we are in average mode (default) or maximum
	# mode (with a -m)
	if args[0] in ('-m', '--maximum'):
		mode = 'max'
		args = args[1:]
		print "Working in maximum power mode"
		
	elif args[0] in ('-h', '--help'):
		print """Given a collection of NPZ files created by drxWatefall, average them down to 
a single spectrum and save them to a common NPZ file.  The output filename is 
always 'aggregated-waterfall.npz' and the file is saved in the current 
directory.

Usage:
  aggregateWaterfalls.py [OPTIONS] npz_file_1 npz_file_2 npz_file_3 ...

Options:
  -h, --help       Show this help message.
  -m, --maximum    Save the spectra associated with the maximum power in
                   each beam/tuning/polarization."""
		
		sys.exit(0)
		
	else:
		mode = 'avg'
		print "Working in averaging mode"
	
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
	if mode == 'avg':
		tIntActual = spec.shape[0] * dataDict['tInt']
	else:
		tIntActual = tIntOriginal
	
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
		if mode == 'avg':
			spec[i,:,:] = cSpec.mean(axis=0)
		else:
			power = cSpec.sum(axis=2)
			for j in xrange(cSpec.shape[1]):
				localMax = numpy.where( power[:,j] == power[:,j].max() )[0][0]
				spec[i,j,:] = cSpec[localMax,j,:]
		
	# Compute the "fake" integration time that is needed by plotWaterfall to 
	# map points correctly on a click.
	tInt = times[1] - times[0]
	
	# Save
	outname = 'aggregated-waterfall.npz'
	numpy.savez(outname, freq=freq, times=times, spec=spec, tInt=tInt, tIntActual=tIntActual, tIntOriginal=tIntOriginal, srate=srate, standMapper=standMapper, filenames=filenames)


if __name__ == "__main__":
	main(sys.argv[1:])


