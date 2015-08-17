#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given an HDF5 file, apply incoherent dedispersion to the data at the specified
DM and save the results to a new file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import time
import numpy
import getopt
import shutil
from datetime import datetime

from lsl.misc.dedispersion import delay, incoherent


def usage(exitCode=None):
	print """dedisperseHDF.py - Read in DRX/HDF5 waterfall file and apply incoherent 
dedispersion to the data.

Usage: dedisperseHDF.py [OPTIONS] DM file

Options:
-h, --help                Display this help information
-c, --correct-time        Correct the timestamps relative to the infinite 
                          frequency case (default = no)
-f, --force               Force overwritting of existing HDF5 files
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['correctTime'] = False
	config['force'] = False
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hcf", ["help", "correct-time", "force"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
		
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-c', '--correct-time'):
			config['correctTime'] = True
		elif opt in ('-f', '--force'):
			config['force'] = True
		else:
			assert False
			
	# Add in arguments
	config['args'] = args
	
	# Validate
	if len(config['args']) != 2:
		raise RuntimeError("Must specify both a DM in pc/cm^3 and a HDF5 file to process")
		
	# Return configuration
	return config


def main(args):
	config = parseOptions(args)
	
	dm = float(config['args'][0])
	filename = config['args'][1]
	
	# Ready the outut filename
	outname = os.path.basename(filename)
	outname = os.path.splitext(outname)[0]
	outname = "%s-DM%.4f.hdf5" % (outname, dm)
	
	if os.path.exists(outname):
		if not config['force']:
			yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
		else:
			yn = 'y'
			
		if yn not in ('n', 'N'):
			shutil.copy(filename, outname)
		else:
			raise RuntimeError("Output file '%s' already exists" % outname)
	else:
		shutil.copy(filename, outname)
		
	# Open the file and loop over the observation sections
	h = h5py.File(outname, mode='a')
	for obsName in h.keys():
		obs = h.get(obsName, None)
		
		# Load in the information we need to calculate the pseudo-spectral kurtosis (pSK)
		tInt = obs.attrs['tInt']
		LFFT = obs.attrs['LFFT']
		srate = obs.attrs['sampleRate']
		
		print "Staring Observation #%i" % int(obsName.replace('Observation', ''))
		print "  Sample Rate: %.1f Hz" % srate
		print "  LFFT: %i" % LFFT
		print "  tInt: %.3f s" % tInt
		print "  DM: %.3f pc/cm^3" % dm
		
		time = obs.get('time', None)[:]
		tuning1 = obs.get('Tuning1', None)
		tuning2 = obs.get('Tuning2', None)
		
		# Get the frequencies
		freq1 = tuning1.get('freq', None)
		freq2 = tuning2.get('freq', None)
		
		# Get the supplemental information
		baseMask1 = tuning1.get('Mask', None)
		baseMask2 = tuning2.get('Mask', None)
		baseSK1 = tuning1.get('SpectralKurtosis', None)
		baseSK2 = tuning2.get('SpectralKurtosis', None)
		
		# Update the time
		if config['correctTime']:
			maxFreq = max( [max(freq1), max(freq2)] )
			infDelay = delay(numpy.array([maxFreq, numpy.inf]), dm)
			print "  Infinite Time Delay: %.3f s" % max(infDelay)
			
			time = time - max(infDelay)
			obs['time'][:] = time
			
		# Loop over data products
		for dp in ('XX', 'YY', 'XY', 'YX', 'I', 'Q', 'U', 'V'):
			data1 = tuning1.get(dp, None)
			data2 = tuning2.get(dp, None)
			if data1 is None:
				continue
				
			if baseMask1 is not None:
				mask1 = baseMask1.get(dp, None)
				mask2 = baseMask2.get(dp, None)
			else:
				mask1 = None
				mask2 = None
			if baseSK1 is not None:
				sk1 = baseSK1.get(dp, None)
				sk2 = baseSK2.get(dp, None)
			else:
				sk1 = None
				sk2 = None
				
			## Combine the frequency segments so that the dedispersion is applied relative 
			## to the same frequency for both tunings
			freqC = numpy.zeros(freq1.shape[0]+freq2.shape[0], dtype=freq1.dtype)
			freqC[:freq1.shape[0]] = freq1[:]
			freqC[freq1.shape[0]:] = freq2[:]
			
			dataC = numpy.zeros((data1.shape[0], freqC.size), dtype=data1.dtype)
			dataC[:,:freq1.shape[0]] = data1[:,:]
			dataC[:,freq1.shape[0]:] = data2[:,:]
			
			if mask1 is not None:
				maskC = numpy.zeros((mask1.shape[0], freqC.size), dtype=mask1.dtype)
				maskC[:,:freq1.shape[0]] = mask1[:,:]
				maskC[:,freq1.shape[0]:] = mask2[:,:]
			if sk1 is not None:
				skC = numpy.zeros((sk1.shape[0], freqC.size), dtype=sk1.dtype)
				skC[:,:freq1.shape[0]] = sk1[:,:]
				skC[:,freq1.shape[0]:] = sk2[:,:]
				
			## Dedisperse
			try:
				dataCD = incoherent(freqC, dataC, tInt, dm, boundary='fill', fill_value=numpy.nan)
			except TypeError:
				dataCD = incoherent(freqC, dataC, tInt, dm)
				
			if mask1 is not None:
				try:
					maskCD = incoherent(freqC, maskC, tInt, dm, boundary='fill', fill_value=numpy.nan)
				except TypeError:
					maskCD = incoherent(freqC, maskC, tInt, dm)
					
				maskCD[numpy.where( ~numpy.isfinite(dataCD) )] = True
				maskCD = maskCD.astype(mask1.dtype)
			if sk1 is not None:
				try:
					skCD = incoherent(freqC, skC, tInt, dm, boundary='fill', fill_value=numpy.nan)
				except TypeError:
					skCD = incoherent(freqC, skC, tInt, dm)
					
				skCD[numpy.where( ~numpy.isfinite(dataCD) )] = 0.0
				skCD = skCD.astype(sk1.dtype)
				
			## Update
			data1[:,:] = dataCD[:,:freq1.shape[0]]
			data2[:,:] = dataCD[:,freq1.shape[0]:]
			if mask1 is not None:
				mask1[:,:] = maskCD[:,:freq1.shape[0]]
				mask2[:,:] = maskCD[:,freq1.shape[0]:]
			if sk1 is not None:
				sk1[:,:] = skCD[:,:freq1.shape[0]]
				sk2[:,:] = skCD[:,freq1.shape[0]:]
				
	# Done!
	h.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	
