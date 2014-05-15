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

from lsl.misc.dedispersion import incoherent


def usage(exitCode=None):
	print """dedisperseHDF.py - Read in DRX/HDF5 waterfall file and apply incoherent 
dedispersion to the data.

Usage: dedisperseHDF.py [OPTIONS] DM file

Options:
-h, --help                Display this help information
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "h", ["help",])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
		
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
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
		yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
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
		
		# Loop over data products
		baseMask1 = tuning1.get('Mask', None)
		baseMask2 = tuning2.get('Mask', None)
		
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
				
			## Dedisperse
			dataCD = incoherent(freqC, dataC, tInt, dm)
			if mask1 is not None:
				maskCD = incoherent(freqC, maskC, tInt, dm)
				maskCD = maskCD.astype(mask1.dtype)
				
			## Update
			data1[:,:] = dataCD[:,:freq1.shape[0]]
			data2[:,:] = dataCD[:,freq1.shape[0]:]
			if mask1 is not None:
				mask1[:,:] = maskCD[:,:freq1.shape[0]]
				mask2[:,:] = maskCD[:,freq1.shape[0]:]
				
	# Done!
	h.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	
