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

from scipy.stats import scoreatpercentile

from lsl.statistics import kurtosis


def main(args):
	mode = 'avg'
	flagRFI = False

	# Figure out if we are in average mode (default) or maximum
	# mode (with a -m)
	if args[0] in ('-m', '--maximum'):
		mode = 'max'
		args = args[1:]
		print "Working in maximum power mode"

	elif args[0] in ('-f' '--flag-rfi'):
		flagRFI = True
		args = args[1:]
		
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
	try:
		freq1 = dataDict['freq1']
		freq2 = dataDict['freq2']
	except KeyError:
		freq1 = freq
		freq2 = freq
	
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
	p50  = numpy.zeros((len(filenames), spec.shape[1], freq.size))
	p90  = numpy.zeros((len(filenames), spec.shape[1], freq.size))
	p99  = numpy.zeros((len(filenames), spec.shape[1], freq.size))
	p100 = numpy.zeros((len(filenames), spec.shape[1], freq.size))
	mask = numpy.zeros((len(filenames), spec.shape[1], freq.size), dtype=numpy.bool)
	for i,filename in enumerate(filenames):
		print '%i of %i' % (i, len(filenames))
		dataDict = numpy.load(filename)
		
		cSpec = dataDict['spec']
		cTime = dataDict['times']
		
		times[i] = cTime[0]
		if mode == 'avg':
			if flagRFI:
				# S-K RFI identification
				kurtosisCut = 4 
	
				N = srate/(freq.size+1) * tIntOriginal
				kurtosis = numpy.zeros((spec.shape[1], spec.shape[2]))
		
				for k in xrange(spec.shape[1]):
					for j in xrange(spec.shape[2]):
						channel = cSpec[:,k,j]
						kurtosis[k,j] = kurtosis.spectralPower(channel, N=N)
		
				kMean = 1.0
				kStd  = kurtosis.std(cSpec.shape[0], N)
				
				for k in xrange(spec.shape[1]):
					bad = numpy.where( numpy.abs(kurtosis[k,:] - kMean) >= kurtosisCut*kStd )[0]
					for b in bad:
						try:
							for j in xrange(b-2, b+3):
								mask[i,k,j] = True
						except IndexError:
							pass

			spec[i,:,:] = cSpec.mean(axis=0)
			p50[i,:,:] = numpy.median(cSpec, axis=0)
			p90[i,:,:] = scoreatpercentile(cSpec, 90)
			p99[i,:,:] = scoreatpercentile(cSpec, 99)
			p100[i,:,:] = numpy.max(cSpec, axis=0)
					
		else:
			power = cSpec.sum(axis=2)
			for j in xrange(cSpec.shape[1]):
				localMax = numpy.where( power[:,j] == power[:,j].max() )[0][0]
				spec[i,j,:] = cSpec[localMax,j,:]
		
	# Compute the "fake" integration time that is needed by plotWaterfall to 
	# map points correctly on a click.
	tInt = times[1] - times[0]
	
	# Convert filenames to absolute paths
	filenames = [os.path.abspath(f) for f in filenames]
	
	p50 = numpy.median(p50, axis=0)
	p90 = scoreatpercentile(p90, 90)
	p99 = scoreatpercentile(p99, 99)
	p100 = numpy.max(p100, axis=0)
	
	# Save
	outname = 'aggregated-waterfall.npz'
	numpy.savez(outname, freq=freq, freq1=freq1, freq2=freq2, times=times, spec=spec, mask=mask, 
				tInt=tInt, tIntActual=tIntActual, tIntOriginal=tIntOriginal, 
				p50=p50, p90=p90, p99=p99, p100=p100, 
				srate=srate, standMapper=standMapper, filenames=filenames)


if __name__ == "__main__":
	main(sys.argv[1:])


