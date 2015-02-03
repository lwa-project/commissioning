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
-g, --generate-mask       Added the pSK information to the file as a mask rather
                          than pSK values.  The masking parameters are 
                          controlled by the -t/--threshold, -f/--fill, and 
                          -m/--merge flags.
-t, --threshold           Masking threshold, in sigma, if the -g/--generate-mask 
                          flag is set (default = 4.0)
-f, --fill                Fill in the secondary polarizations (XY, YX; Q, U, V; 
                          RL, LR) using the mask from the primaries (XX, YY; I; 
                          RR, LL) (default = no)
-m, --merge               Merge the new mask table with the existing one, if it
                          exists (default = replace)
                          
Note: The -g/--generate-mask flag overrides the -n/--no-update flag.
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
	config['mask'] = False
	config['threshold'] = 4.0
	config['fill'] = False
	config['merge'] = False
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hd:ngt:fm", ["help", "duration=", "no-update", "generate-mask", "threshold=", "fill", "merge"])
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
		elif opt in ('-g', '--generate-mask'):
			config['update'] = True
			config['mask'] = True
		elif opt in ('-t', '--threshold'):
			config['threshold'] = float(value)
		elif opt in ('-f', '--fill'):
			config['fill'] = True
		elif opt in ('-m', '--merge'):
			config['merge'] = True
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
		
		# Get all of the avaliable data products
		dataProducts = list(tuning1)
		for toRemove in ('Mask', 'Saturation', 'SpectralKurtosis', 'freq'):
			try:
				del dataProducts[dataProducts.index(toRemove)]
			except ValueError:
				pass
				
		# What's this?  Stokes I is the sum of XX and YY so there are 
		# effectively twice as many FFTs per spectrum relative to XX and
		# YY
		nAdjust = {'XX': 1, 'YY': 1, 'I': 2, 'RR': 1, 'LL': 1}
		
		# Loop over data products - primary
		for dp in ('XX', 'YY', 'I', 'RR', 'LL'):
			data1 = tuning1.get(dp, None)
			data2 = tuning2.get(dp, None)
			if data1 is None:
				continue
				
			# Loop over chunks for computing pSK
			sk1 = numpy.zeros_like(data1)
			sk2 = numpy.zeros_like(data2)
			
			section1 = numpy.empty((chunkSize,data1.shape[1]), dtype=data1.dtype)
			section2 = numpy.empty((chunkSize,data2.shape[1]), dtype=data2.dtype)
			for i in xrange(time.size/chunkSize+1):
				start = i*chunkSize
				stop = start + chunkSize
				if stop >= time.size:
					stop = time.size - 1
				if start == stop:
					continue
					
				## Compute the pSK for each channel in both tunings
				selection = numpy.s_[start:stop, :]
				try:
					data1.read_direct(section1, selection)
					data2.read_direct(section2, selection)
				except TypeError:
					section1 = data1[start:stop,:]
					section2 = data2[start:stop,:]
				for j in xrange(section1.shape[1]):
					sk1[start:stop,j] = kurtosis.spectralPower(section1[:,j], N=skN*nAdjust[dp])
					sk2[start:stop,j] = kurtosis.spectralPower(section2[:,j], N=skN*nAdjust[dp])
					
			# Report
			print "  => %s-1 SK Mean: %.3f" % (dp, numpy.mean(sk1))
			print "     %s-1 SK Std. Dev.: %.3f" % (dp, numpy.std(sk1))
			print "     %s-2 SK Mean: %.3f" % (dp, numpy.mean(sk2))
			print "     %s-2 SK Std. Dev.: %.3f" % (dp, numpy.std(sk2))
			
			# Save the pSK information to the HDF5 file if we need to
			if config['update']:
				h.attrs['FileLastUpdated'] = datetime.utcnow().strftime("UTC %Y/%m/%d %H:%M:%S")
				
				if config['mask']:
					## Calculate the expected pSK limits for the threshold
					kLow, kHigh = kurtosis.getLimits(config['threshold'], chunkSize, N=skN*nAdjust[dp])
					
					## Generate the mask arrays
					maskTable1 = numpy.where( (sk1 < kLow) | (sk1 > kHigh), True, False )
					maskTable2 = numpy.where( (sk2 < kLow) | (sk2 > kHigh), True, False )
					
					## Report
					print "       => %s-1 Mask Fraction: %.1f%%" % (dp, 100.0*maskTable1.sum()/maskTable1.size,)
					print "          %s-2 Mask Fraction: %.1f%%" % (dp, 100.0*maskTable2.sum()/maskTable2.size,)
					
					## Pull out the Mask group from the HDF5 file
					### Tuning 1
					mask1 = tuning1.get('Mask', None)
					if mask1 is None:
						mask1 = tuning1.create_group('Mask')
						for p in dataProducts:
							mask1.create_dataset(p, tuning1[p].shape, 'bool')
					mask1DP = mask1.get(dp, None)
					if config['merge']:
						mask1DP[:,:] |= maskTable1.astype(numpy.bool)
					else:
						mask1DP[:,:] = maskTable1.astype(numpy.bool)
					
					### Tuning 2
					mask2 = tuning2.get('Mask', None)
					if mask2 is None:
						mask2 = tuning2.create_group('Mask')
						for p in dataProducts:
							mask2.create_dataset(p, tuning2[p].shape, 'bool')
					mask2DP = mask2.get(dp, None)
					if config['merge']:
						mask2DP[:,:] |= maskTable2.astype(numpy.bool)
					else:
						mask2DP[:,:] = maskTable2.astype(numpy.bool)
						
				else:
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
					
		if config['mask'] and config['fill']:	
			# Loop over data products - secondary
			for dp in dataProducts:
				## Jump over the primary polarizations
				if dp in ('XX', 'YY', 'I', 'RR', 'LL'):
					continue
					
				if dp in ('Q', 'U', 'V'):
					## Case 1 - Flag Stokes Q, U, and V off I
					
					## Pull out the Mask group from the HDF5 file
					### Tuning 1
					mask1 = tuning1.get('Mask', None)
					mask1DP = mask1.get(dp, None)
					mask1Base = mask1.get('I', None)
					if config['merge']:
						mask1DP[:,:] |= mask1Base[:,:]
					else:
						mask1DP[:,:] = mask1Base[:,:]
						
					### Tuning 2
					mask2 = tuning2.get('Mask', None)
					mask2DP = mask2.get(dp, None)
					mask2Base = mask2.get('I', None)
					if config['merge']:
						mask2DP[:,:] |= mask2Base[:,:]
					else:
						mask2DP[:,:] = mask2Base[:,:]
						
				elif dp in ('XY', 'YX'):
					## Case 2 - Flag XY and YX off XX and YY
					
					## Pull out the Mask group from the HDF5 file
					### Tuning 1
					mask1 = tuning1.get('Mask', None)
					mask1DP = mask1.get(dp, None)
					mask1Base1 = mask1.get('XX', None)
					mask1Base2 = mask1.get('YY', None)
					if config['merge']:
						mask1DP[:,:] |= (mask1Base1[:,:] | mask1Base2[:,:])
					else:
						mask1DP[:,:] = (mask1Base1[:,:] | mask1Base2[:,:])
						
					### Tuning 2
					mask2 = tuning2.get('Mask', None)
					mask2DP = mask2.get(dp, None)
					mask2Base1 = mask2.get('XX', None)
					mask2Base2 = mask2.get('YY', None)
					if config['merge']:
						mask2DP[:,:] |= (mask2Base1[:,:] | mask2Base2[:,:])
					else:
						mask2DP[:,:] = (mask2Base1[:,:] | mask2Base2[:,:])
						
				elif dp in ('RL', 'LR'):
					## Case 3 - Flag RL and LR off RR and LL
					
					## Pull out the Mask group from the HDF5 file
					### Tuning 1
					mask1 = tuning1.get('Mask', None)
					mask1DP = mask1.get(dp, None)
					mask1Base1 = mask1.get('RR', None)
					mask1Base2 = mask1.get('LL', None)
					if config['merge']:
						mask1DP[:,:] |= (mask1Base1[:,:] | mask1Base2[:,:])
					else:
						mask1DP[:,:] = (mask1Base1[:,:] | mask1Base2[:,:])
						
					### Tuning 2
					mask2 = tuning2.get('Mask', None)
					mask2DP = mask2.get(dp, None)
					mask2Base1 = mask2.get('RR', None)
					mask2Base2 = mask2.get('LL', None)
					if config['merge']:
						mask2DP[:,:] |= (mask2Base1[:,:] | mask2Base2[:,:])
					else:
						mask2DP[:,:] = (mask2Base1[:,:] | mask2Base2[:,:])
						
	# Done!
	h.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	
