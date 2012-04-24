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


def usage(exitCode=None):
	print """drspec2hdf.py - Convert a DR spectrometer file to a HDF5 similar to what
hdfWaterfall.py generates.

Usage: drspec2hdf.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-n, --normalize             Normalize the spectrometer data by the number of fills
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['normalize'] = False
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hn", ["help", "normalize"])
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
		elif opt in ('-n', '--normalize'):
			config['normalize'] = True
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
	
	# Report
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % beginDate
	print "Beam: %i" % beam
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Frames: %i (%.3f s)" % (nFrames, nFrames*tInt)
	print "---"
	print "Transform Length: %i" % LFFT
	print "Integration: %.3f s" % tInt
	
	# Setup the output file
	outname = filename.replace('.dat', '-waterfall.hdf5')
		
	f = h5py.File(outname, 'w')
	f.attrs['Beam'] = beam
	f.attrs['tInt'] = tInt
	f.attrs['tInt_Units'] = 's'
	f.attrs['sampleRate'] = srate
	f.attrs['sampleRate_Units'] = 'samples/s'
	freq = numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) )
	freq = freq[1:].astype(numpy.float64)
	f.attrs['RBW'] = freq[1]-freq[0]
	f.attrs['RBW_Units'] = 'Hz'
	masterTimes = f.create_dataset('time', (nChunks,), 'f8')
		
	tuning1 = f.create_group('/Tuning1')
	tuning1['freq'] = freq + centralFreq1
	tuning1['freq'].attrs['Units'] = 'Hz'
	spec1X = tuning1.create_dataset('X', (nChunks, LFFT-1), 'f4', chunks=True)
	tuning1['X'].attrs['axis0'] = 'time'
	tuning1['X'].attrs['axis1'] = 'frequency'
	spec1Y = tuning1.create_dataset('Y', (nChunks, LFFT-1), 'f4', chunks=True)
	tuning1['Y'].attrs['axis0'] = 'time'
	tuning1['Y'].attrs['axis1'] = 'frequency'

	tuning2 = f.create_group('/Tuning2')
	tuning2['freq'] = freq + centralFreq2
	tuning2['freq'].attrs['Units'] = 'Hz'
	spec2X = tuning2.create_dataset('X', (nChunks, LFFT-1), 'f4', chunks=True)
	tuning2['X'].attrs['axis0'] = 'time'
	tuning2['X'].attrs['axis1'] = 'frequency'
	spec2Y = tuning2.create_dataset('Y', (nChunks, LFFT-1), 'f4', chunks=True)
	tuning2['Y'].attrs['axis0'] = 'time'
	tuning2['Y'].attrs['axis1'] = 'frequency'

	# Loop over DR spectrometer frames to fill in the HDF5 file
	for i in xrange(nChunks):
		frame = drspec.readFrame(fh)
		
		masterTimes[i] = frame.getTime()
		
		if config['normalize']:
			spec1X[i,:] = frame.data.X0[1:] / LFFT / frame.data.fills[0]
			spec1Y[i,:] = frame.data.Y0[1:] / LFFT / frame.data.fills[1]
			spec2X[i,:] = frame.data.X1[1:] / LFFT / frame.data.fills[2]
			spec2Y[i,:] = frame.data.Y1[1:] / LFFT / frame.data.fills[3]
		else:
			spec1X[i,:] = frame.data.X0[1:] / LFFT
			spec1Y[i,:] = frame.data.Y0[1:] / LFFT
			spec2X[i,:] = frame.data.X1[1:] / LFFT
			spec2Y[i,:] = frame.data.Y1[1:] / LFFT

	# Done
	fh.close()

	# Save the output to a HDF5 file
	f.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	
