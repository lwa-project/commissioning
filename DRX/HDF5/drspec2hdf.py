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

from lsl.reader import drx, drspec, errors
from lsl.common import progress

import data as hdfData


def usage(exitCode=None):
	print """drspec2hdf.py - Convert a DR spectrometer file to a HDF5 similar to what
hdfWaterfall.py generates.

Usage: drspec2hdf.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-s, --skip                  Skip foward the specified number of seconds into the file
-m, --metadata              Metadata tarball for additional information
-d, --sdf                   SDF for additional information
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0.0
	config['metadata'] = None
	config['sdf'] = None
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hs:m:d:", ["help", "skip=", "metadata=", "sdf="])
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
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-m', '--metadata'):
			config['metadata'] = value
		elif opt in ('-d', '--sdf'):
			config['sdf'] = value
		else:
			assert False
			
	# Make sure we aren't offsetting when we have metadata
	if config['metadata'] is not None or config['sdf'] is not None:
		config['offset'] = 0.0
		
	# Add in arguments
	config['args'] = args
	
	# Return configuration
	return config


def main(args):
	config = parseOptions(args)

	# Open the file and file good data (not raw DRX data)
	filename = config['args'][0]
	fh = open(filename, 'rb')

	try:
		for i in xrange(5):
			junkFrame = drx.readFrame(fh)
		raise RuntimeError("ERROR: '%s' appears to be a raw DRX file, not a DR spectrometer file" % filename)
	except errors.syncError:
		fh.seek(0)
		
	# Interrogate the file to figure out what frames sizes to expect, now many 
	# frames there are, and what the transform length is
	FrameSize = drspec.getFrameSize(fh)
	nFrames = os.path.getsize(filename) / FrameSize
	nChunks = nFrames
	LFFT = drspec.getTransformSize(fh)

	# Read in the first frame to figure out the DP information
	junkFrame = drspec.readFrame(fh)
	fh.seek(-FrameSize, 1)
	srate = junkFrame.getSampleRate()
	t0 = junkFrame.getTime()
	tInt = junkFrame.header.nInts*LFFT/srate
	
	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(round(config['offset'] / tInt))
	fh.seek(offset*FrameSize, 1)
	
	# Iterate on the offsets until we reach the right point in the file.  This
	# is needed to deal with files that start with only one tuning and/or a 
	# different sample rate.  
	while True:
		## Figure out where in the file we are and what the current tuning/sample 
		## rate is
		junkFrame = drspec.readFrame(fh)
		srate = junkFrame.getSampleRate()
		t1 = junkFrame.getTime()
		tInt = junkFrame.header.nInts*LFFT/srate
		fh.seek(-FrameSize, 1)
		
		## See how far off the current frame is from the target
		tDiff = t1 - (t0 + config['offset'])
		
		## Half that to come up with a new seek parameter
		tCorr = -tDiff / 2.0
		cOffset = int(round(tCorr / tInt))
		offset += cOffset
		
		## If the offset is zero, we are done.  Otherwise, apply the offset
		## and check the location in the file again/
		if cOffset is 0:
			break
		fh.seek(cOffset*FrameSize, 1)
		
	# Update the offset actually used
	config['offset'] = t1 - t0
	nChunks = (os.path.getsize(filename) - fh.tell()) / FrameSize
	
	# Update the file contents
	beam = junkFrame.parseID()
	centralFreq1 = junkFrame.getCentralFreq(1)
	centralFreq2 = junkFrame.getCentralFreq(2)
	srate = junkFrame.getSampleRate()
	dataProducts = junkFrame.getDataProducts()
	t0 = junkFrame.getTime()
	tInt = junkFrame.header.nInts*LFFT/srate
	beginDate = datetime.utcfromtimestamp(junkFrame.getTime())
        
	# Report
	print "Filename: %s" % filename
	if config['metadata'] is not None:
		print "Metadata: %s" % config['metadata']
	elif config['sdf'] is not None:
		print "SDF: %s" % config['sdf']
	print "Date of First Frame: %s" % beginDate
	print "Beam: %i" % beam
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Data Products: %s" % ','.join(dataProducts)
	print "Frames: %i (%.3f s)" % (nFrames, nFrames*tInt)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Transform Length: %i" % LFFT
	print "Integration: %.3f s" % tInt
	
	# Setup the output file
	outname = os.path.split(filename)[1]
	outname = os.path.splitext(outname)[0]
	outname = '%s-waterfall.hdf5' % outname
	
	if os.path.exists(outname):
		yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
		if yn not in ('n', 'N'):
			os.unlink(outname)
		else:
			raise RuntimeError("Output file '%s' already exists" % outname)
			
	f = hdfData.createNewFile(outname)
	obsList = {}
	if config['metadata'] is not None:
		from lsl.common import mcs, metabundle
		sdf = metabundle.getSessionDefinition(config['metadata'])
		
		sdfBeam  = sdf.sessions[0].drxBeam
		spcSetup = sdf.sessions[0].spcSetup
		if sdfBeam != beam:
			raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
			
		for i,obs in enumerate(sdf.sessions[0].observations):
			sdfStart = mcs.mjdmpm2datetime(obs.mjd, obs.mpm)
			sdfStop  = mcs.mjdmpm2datetime(obs.mjd, obs.mpm + obs.dur)
			obsChunks = int(numpy.ceil(obs.dur/1000.0 * drx.filterCodes[obs.filter] / (spcSetup[0]*spcSetup[1])))
			
			obsList[i+1] = (sdfStart, sdfStop, obsChunks)
			
		hdfData.fillFromMetabundle(f, config['metadata'])
		
	elif config['sdf'] is not None:
		from lsl.common import mcs
		from lsl.common.sdf import parseSDF
		sdf = parseSDF(config['sdf'])
		
		sdfBeam  = sdf.sessions[0].drxBeam
		spcSetup = sdf.sessions[0].spcSetup
		if sdfBeam != beam:
			raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
			
		for i,obs in enumerate(sdf.sessions[0].observations):
			sdfStart = mcs.mjdmpm2datetime(obs.mjd, obs.mpm)
			sdfStop  = mcs.mjdmpm2datetime(obs.mjd, obs.mpm + obs.dur)
			obsChunks = int(numpy.ceil(obs.dur/1000.0 * drx.filterCodes[obs.filter] / (spcSetup[0]*spcSetup[1])))
			
			obsList[i+1] = (sdfStart, sdfStop, obsChunks)
			
		hdfData.fillFromSDF(f, config['sdf'])
		
	else:
		obsList[1] = (beginDate, datetime(2222,12,31,23,59,59), nChunks)
		
		hdfData.fillMinimum(f, 1, beam, srate)
		
	dataProducts = junkFrame.getDataProducts()
	for o in sorted(obsList.keys()):
		for t in (1,2):
			hdfData.createDataSets(f, o, t, numpy.arange(LFFT, dtype=numpy.float32), obsList[o][2], dataProducts)
			
	f.attrs['FileGenerator'] = 'drspec2hdf.py'
	f.attrs['InputData'] = os.path.basename(filename)
	
	# Create the various HDF group holders
	ds = {}
	for o in sorted(obsList.keys()):
		obs = hdfData.getObservationSet(f, o)
		
		ds['obs%i' % o] = obs
		ds['obs%i-time' % o] = obs.create_dataset('time', (obsList[o][2],), 'f8')
		
		for t in (1,2):
			ds['obs%i-freq%i' % (o, t)] = hdfData.getDataSet(f, o, t, 'freq')
			for p in dataProducts:
				ds["obs%i-%s%i" % (o, p, t)] = hdfData.getDataSet(f, o, t, p)
			ds['obs%i-Saturation%i' % (o, t)] = hdfData.getDataSet(f, o, t, 'Saturation')
			
	# Loop over DR spectrometer frames to fill in the HDF5 file
	pbar = progress.ProgressBar(max=nChunks)
	o = 1
	j = 0
	
	firstPass = True
	for i in xrange(nChunks):
		frame = drspec.readFrame(fh)
		
		cTime = datetime.utcfromtimestamp(frame.getTime())
		if cTime < obsList[o][0]:
			# Skip over data that occurs before the start of the observation
			continue
		elif cTime > obsList[o][1]:
			# Increment to the next observation
			o += 1
			
			# If we have reached the end, exit...
			try:
				obsList[o]
				
				firstPass = True
			except KeyError:
				sys.stdout.write('%s\r' % (' '*pbar.span))
				sys.stdout.flush()
				print "End of observing block according to SDF, exiting"
				break
				
		else:
			pass
			
		try:
			if frame.getTime() > oTime + 1.001*tInt:
				print 'Warning: Time tag error at frame %i; %.3f > %.3f + %.3f' % (i, frame.getTime(), oTime, tInt)
		except NameError:
			pass
		oTime = frame.getTime()
		
		if firstPass:
			# Otherwise, continue on...
			centralFreq1 = frame.getCentralFreq(1)
			centralFreq2 = frame.getCentralFreq(2)
			srate = frame.getSampleRate()
			tInt  = frame.header.nInts*LFFT/srate
			
			freq = numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) )
			freq = freq.astype(numpy.float64)
			
			sys.stdout.write('%s\r' % (' '*pbar.span))
			sys.stdout.flush()
			print "Switching to Obs. #%i" % o
			print "-> Tunings: %.1f Hz, %.1f Hz" % (centralFreq1, centralFreq2)
			print "-> Sample Rate: %.1f Hz" % srate
			print "-> Integration Time: %.3f s" % tInt
			sys.stdout.write(pbar.show()+'\r')
			sys.stdout.flush()
			
			j = 0
			ds['obs%i-freq1' % o][:] = freq + centralFreq1
			ds['obs%i-freq2' % o][:] = freq + centralFreq2
			
			obs = ds['obs%i' % o]
			obs.attrs['tInt'] = tInt
			obs.attrs['tInt_Units'] = 's'
			obs.attrs['LFFT'] = LFFT
			obs.attrs['nChan'] = LFFT
			obs.attrs['RBW'] = freq[1]-freq[0]
			obs.attrs['RBW_Units'] = 'Hz'
			
			firstPass = False
			
		# Load the data from the spectrometer frame into the HDF5 group
		ds['obs%i-time' % o][j] = frame.getTime()
		
		ds['obs%i-Saturation1' % o][j,:] = frame.data.saturations[0:2]
		ds['obs%i-Saturation2' % o][j,:] = frame.data.saturations[2:4]
		
		for t in (1,2):
			for p in dataProducts:
				ds['obs%i-%s%i' % (o, p, t)][j,:] = getattr(frame.data, "%s%i" % (p, t-1), None)
		j += 1
		
		# Update the progress bar
		pbar.inc()
		if i % 10 == 0:
			sys.stdout.write(pbar.show()+'\r')
			sys.stdout.flush()
			
	sys.stdout.write(pbar.show()+'\n')
	sys.stdout.flush()
	
	# Done
	fh.close()

	# Save the output to a HDF5 file
	f.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	
