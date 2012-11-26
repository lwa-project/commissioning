#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Small script to read in a DR spectrometer binary data file and create a HDF5 in 
the image of hdfWaterfall.py that can be plotted with plotHDF.py

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


def usage(exitCode=None):
	print """drspec2hdf.py - Convert a DR spectrometer file to a HDF5 similar to what
hdfWaterfall.py generates.

Usage: drspec2hdf.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-t, --integrate-time        Resample the time resolution
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['intTime'] = None
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "ht:", ["help", "integrate-time="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-q', '--quiet'):
			config['verbose'] = False
		elif opt in ('-t', '--integrate-time'):
			config['intTime'] = float(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)

	filename = config['args'][0]
	fh = open(filename, 'rb')

	# Interogate the file to figure out what frames sizes to expect, now many 
	# frames there are, and what the transform length is
	nFrames = os.path.getsize(filename) / drspec.getFrameSize(fh)
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
	beginDate = datetime.utcfromtimestamp(junkFrame.getTime())
	
	products = junkFrame.getDataProducts()
	nProducts = len(products)
	
	# Conversions
	if config['intTime'] is None:
		nSubInts = 1
	else:
		nSubInts = int(round(1.0 * config['intTime'] / tInt))
	config['intTime'] = nSubInts * tInt
	
	print nFrames, tInt, 1.0*nFrames*tInt
	
	# Report
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % beginDate
	print "Beam: %i" % beam
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Data Products: %s" % (', '.join(products),)
	print "Frames: %i (%.3f s)" % (nFrames, nFrames*tInt)
	print "---"
	print "Transform Length: %i" % LFFT
	print "Integration: %.3f s" % tInt
	if config['intTime'] != tInt:
		print "-> Resampling to %.3f s integrations" % config['intTime']
	
	# Setup the output file
	outname = filename.replace('.dat', '-waterfall.hdf5')
		
	f = h5py.File(outname, 'w')
	f.attrs['Beam'] = beam
	f.attrs['tInt'] = config['intTime']
	f.attrs['tInt_Units'] = 's'
	f.attrs['sampleRate'] = srate
	f.attrs['sampleRate_Units'] = 'samples/s'
	freq = numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) )
	freq = freq[1:].astype(numpy.float64)
	f.attrs['RBW'] = freq[1]-freq[0]
	f.attrs['RBW_Units'] = 'Hz'
	masterTimes = f.create_dataset('time', (nChunks/nSubInts,), 'f8')
	
	tuning1 = f.create_group('/Tuning1')
	tuning1['freq'] = freq + centralFreq1
	tuning1['freq'].attrs['Units'] = 'Hz'
	
	tuning2 = f.create_group('/Tuning2')
	tuning2['freq'] = freq + centralFreq2
	tuning2['freq'].attrs['Units'] = 'Hz'
	
	spec = {}
	for p in products:
		k = "%s0" % p
		spec[k] = tuning1.create_dataset(p, (nChunks/nSubInts, LFFT-1), 'f4')
		tuning1[p].attrs['axis0'] = 'time'
		tuning1[p].attrs['axis1'] = 'frequency'
		
		k = "%s1" % p
		spec[k] = tuning2.create_dataset(p, (nChunks/nSubInts, LFFT-1), 'f4')
		tuning2[p].attrs['axis0'] = 'time'
		tuning2[p].attrs['axis1'] = 'frequency'
		
	# Loop over DR spectrometer frames to fill in the HDF5 file
	tempTimes = numpy.zeros(nChunks, dtype=numpy.float64)
	tempArray = numpy.zeros((2*nProducts, nChunks, LFFT-1), dtype=numpy.float32)
	
	pbar = progress.ProgressBar(max=nChunks)
	for i in xrange(nChunks):
		frame = drspec.readFrame(fh)
		
		tempTimes[i] = frame.getTime()
		for j,p in enumerate(products):
			tempArray[0*nProducts+j,i,:] = getattr(frame.data, "%s0" % p, None)[1:]
			tempArray[1*nProducts+j,i,:] = getattr(frame.data, "%s1" % p, None)[1:]
			
		pbar.inc()
		if i % 10 == 0:
			sys.stdout.write(pbar.show()+'\r')
			sys.stdout.flush()
	
	if nSubInts is 1:
		sys.stdout.write(pbar.show()+'\n')
	else:
		sys.stdout.write(pbar.show()+'\r')
	sys.stdout.flush()
	
	# Always normalize by the FFT length
	tempArray /= float(LFFT)
	
	# Save, integrating if necessary
	if nSubInts is 1:
		masterTimes[:] = tempTimes
		
		# Save the results to the HDF5 file
		for j,p in enumerate(products):
			spec['%s0' % p][:,:] = tempArray[0*nProducts+j,:,:]
			spec['%s1' % p][:,:] = tempArray[1*nProducts+j,:,:]
		
	else:
		# Calculate the spans to add
		spans = range(0, (tempArray.shape[1]/nSubInts+1)*nSubInts, nSubInts)
		
		# Sample the time at the start of each span
		masterTimes[:] = tempTimes[spans][:-1]
		
		# Create a temporary array to hold the intermediate results for the integrations and Go!
		tempArray2 = numpy.zeros((tempArray.shape[0], len(spans), LFFT-1), dtype=numpy.float32)
		
		pbar = progress.ProgressBar(max=tempArray.shape[0]*tempArray.shape[2])
		for j in xrange(tempArray.shape[0]):
			for i in xrange(tempArray.shape[2]):
				tempArray2[j,:,i] = numpy.add.reduceat(tempArray[j,:,i], spans, dtype=numpy.float64)
				
				pbar.inc()
				if (j*tempArray.shape[2] + i) % 10 == 0:
					sys.stdout.write(pbar.show()+'\r')
					sys.stdout.flush()
					
		sys.stdout.write(pbar.show()+'\n')
		sys.stdout.flush()
		
		# Save the results to the HDF5 file
		for j,p in enumerate(products):
			spec['%s0' % p][:,:] = tempArray[0*nProducts+j,:-1,:]
			spec['%s1' % p][:,:] = tempArray[1*nProducts+j,:-1,:]
			
	# Done
	fh.close()

	# Save the output to a HDF5 file
	f.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	
