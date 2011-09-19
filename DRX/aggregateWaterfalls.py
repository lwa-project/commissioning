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


def spectralKurtosis(x, N=1):
	"""
	Compute the spectral kurtosis for a set of power measurments averaged
	over N FFTs.  For a distribution consistent with Gaussian noise, this
	value should be ~1.
	"""
	
	M = len(x)
	
	k = M*(x**2).sum()/(x.sum())**2 - 1.0
	k *= (M*N+1)/(M-1)
	
	return k


def skStd(M, N=1):
	"""
	Return the expected standard deviation of the spectral kurtosis for M points 
	each composed of N measurments.
	
	.. note::
		In the future (LSL 0.5), this will be lsl.statistics.kurtosis.std()
	"""
	
	return numpy.sqrt( skVar(M, N) )


def skVar(M, N=1):
	"""
	Return the expected variance (second central moment) of the spectral kurtosis 
	for M points each composed of N measurments.
	.. note::
		In the future (LSL 0.5), this will be lsl.statistics.kurtosis.var()
	"""

	return 2.0*N*(N+1)*M**2/ float( (M-1)*(M*N+3)*(M*N+2) )


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
	mask = numpy.zeros((len(filenames), spec.shape[1], freq.size), dtype=numpy.bool)
	for i,filename in enumerate(filenames):
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
						kurtosis[k,j] = spectralKurtosis(channel, N=N)
		
				kMean = 1.0
				kStd  = skStd(cSpec.shape[0], N)
				
				for k in xrange(spec.shape[1]):
					bad = numpy.where( numpy.abs(kurtosis[k,:] - kMean) >= kurtosisCut*kStd )[0]
					for b in bad:
						try:
							for j in xrange(b-2, b+3):
								mask[i,k,j] = True
						except IndexError:
							pass

				spec[i,:,:] = cSpec.mean(axis=0)
			else:
				spec[i,:,:] = cSpec.mean(axis=0)
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
	
	# Save
	outname = 'aggregated-waterfall.npz'
	numpy.savez(outname, freq=freq, freq1=freq1, freq2=freq2, times=times, spec=spec, mask=mask, tInt=tInt, tIntActual=tIntActual, tIntOriginal=tIntOriginal, srate=srate, standMapper=standMapper, filenames=filenames)


if __name__ == "__main__":
	main(sys.argv[1:])


