#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
calculateSK.py - Read in a DRX/HDF5 wataerfall file and calculate the pseudo-
spectral kurtosis (pSK) for each observation (XX and YY or I).  The pSK values
are "pseudo" since the spectra contained within the HDF5 file are averaged 
over N FFTs where N > 1.  The pSK value for each channel is calculated by 
examining M consecutive spectra.

Note:  Although the pSK method works reasonable well for short integration times
(<~0.3 s) the method may break down for much longer integration times.

Usage:
./calculateSK.py [OPTIONS] file

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

from lsl.statistics import robust, kurtosis


def usage(exitCode=None):
	print """calculateSK.py - Read in DRX/HDF5 waterfall file and calculate the pseudo-
spectral kurtosis (pSK).

Usage: calculateSK.py [OPTIONS] file

Options:
-h, --help                Display this help information
-d, --duration            pSK update interval (default = 10s)
-n, --no-update           Do not add the pSK information to the HDF5 file
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['duration'] = 10.0
	config['update'] = True
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hd:n", ["help", "duration=", "no-update"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
		
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-d', '--duration'):
			config['duration'] = float(value)
		elif opt in ('-n', '--no-update'):
			config['update'] = False
		else:
			assert False
			
	# Add in arguments
	config['args'] = args
	
	# Return configuration
	return config


def main(args):
	config = parseOptions(args)
	
	filename = config['args'][0]
	
	# Open the file and loop over the observation sections
	h = h5py.File(filename, mode='a' if config['update'] else 'r')
	for obsName in h.keys():
		obs = h.get(obsName, None)
		
		# Load in the information we need to calculate the pseudo-spectral kurtosis (pSK)
		tInt = obs.attrs['tInt']
		LFFT = obs.attrs['LFFT']
		srate = obs.attrs['sampleRate']
		
		skN = int(tInt*srate / LFFT)
		chunkSize = int(round(config['duration']/tInt))
		
		print "Staring Observation #%i" % int(obsName.replace('Observation', ''))
		print "  Sample Rate: %.1f Hz" % srate
		print "  LFFT: %i" % LFFT
		print "  tInt: %.3f s" % tInt
		print "  Integrations per spectrum: %i" % skN
		
		time = obs.get('time', None)[:]
		tuning1 = obs.get('Tuning1', None)
		tuning2 = obs.get('Tuning2', None)
		
		# What's this?  Stokes I is the sum of XX and YY so there are 
		# effectively twice as many FFTs per spectrum relative to XX and
		# YY
		nAdjust = {'XX': 1, 'YY': 1, 'I': 2}
		
		# Loop over data products
		for dp in ('XX', 'YY', 'I'):
			data1 = tuning1.get(dp, None)
			data2 = tuning2.get(dp, None)
			if data1 is None:
				continue
			
			# Loop over chunks for computing pSK
			sk1 = numpy.zeros_like(data1)
			sk2 = numpy.zeros_like(data2)
			for i in xrange(time.size/chunkSize+1):
				start = i*chunkSize
				stop = start + chunkSize
				if stop >= time.size:
					stop = time.size - 1
				if start == stop:
					continue
					
				## Compute the pSK for each channel in both tunings
				section1 = data1[start:stop,:]
				section2 = data2[start:stop,:]
				for j in xrange(section1.shape[1]):
					sk1[start:stop,j] = kurtosis.spectralPower(section1[:,j], N=skN*nAdjust[dp])
					sk2[start:stop,j] = kurtosis.spectralPower(section2[:,j], N=skN*nAdjust[dp])
					
			# Report
			print "  => %s-1 SK Mean: %.3f" % (dp, robust.mean(sk1))
			print "     %s-1 SK Std. Dev.: %.3f" % (dp, robust.std(sk1))
			print "     %s-2 SK Mean: %.3f" % (dp, robust.mean(sk2))
			print "     %s-2 SK Std. Dev.: %.3f" % (dp, robust.std(sk2))
			
			# Save the pSK information to the HDF5 file if we need to
			if config['update']:
				h.attrs['FileLastUpdated'] = datetime.utcnow().strftime("UTC %Y/%m/%d %H:%M:%S")
				
				## Tuning 1
				sks1 = tuning1.get('SpectralKurtosis', None)
				if sks1 is None:
					sks1 = tuning1.create_group('SpectralKurtosis')
				try:
					sks1.create_dataset(dp, sk1.shape, 'f4')
				except:
					pass
				sks1[dp][:,:] = sk1
				sks1.attrs['N'] = skN*nAdjust[dp]
				sks1.attrs['M'] = chunkSize
				
				## Tuning 2
				sks2 = tuning2.get('SpectralKurtosis', None)
				if sks2 is None:
					sks2 = tuning2.create_group('SpectralKurtosis')
				try:
					sks2.create_dataset(dp, sk2.shape, 'f4')
				except:
					pass
				sks2[dp][:,:] = sk2
				sks2.attrs['N'] = skN*nAdjust[dp]
				sks2.attrs['M'] = chunkSize
				
	# Done!
	h.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	