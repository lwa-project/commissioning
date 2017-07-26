#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a DRX file, plot the time averaged spectra for each beam output over some 
period.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import math
import numpy
import ephem
import getopt
from datetime import datetime

from lsl.reader import drx, drspec, errors
from lsl.reader.ldp import LWA1DataFile
import lsl.correlator.fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common import progress, stations
from lsl.common import mcs, sdf, metabundle
try:
	from lsl.common import sdfADP, metabundleADP
	adpReady = True
except ImportError:
	adpReady = False

import matplotlib.pyplot as plt

import data as hdfData


def usage(exitCode=None):
	print """hdfWaterfall.py - Read in DRX files and create a collection of 
time-averaged spectra.  These spectra are saved to a HDF5 file called <filename>-waterfall.hdf5.

Usage: hdfWaterfall.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-t, --bartlett              Apply a Bartlett window to the data
-b, --blackman              Apply a Blackman window to the data
-n, --hanning               Apply a Hanning window to the data
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-a, --average               Number of seconds of data to average for spectra 
                            (default = 1)
-d, --duration              Number of seconds to calculate the waterfall for 
                            (default = 0; run the entire file)
-q, --quiet                 Run drxSpectra in silent mode and do not show the plots
-l, --fft-length            Set FFT length (default = 4096)
-c, --clip-level            FFT blanking clipping level in counts (default = 0, 
                            0 disables)
-e, --estimate-clip         Use robust statistics to estimate an appropriate clip 
                            level (overrides the `-c` option)
-m, --metadata              Metadata tarball for additional information
-i, --sdf                   SDF for additional information
-v, --lwasv                 Data is from LWA-SV instead of LWA-1
-f, --force                 Force overwritting of existing HDF5 files
-k, --stokes                Generate Stokes parameters instead of XX and YY
-w, --without-sats          Do not generate saturation counts

Note:  Both the -m/--metadata and -i/--sdf options provide the same additional
       observation information to hdfWaterfall.py so only one needs to be provided.

Note:  Specifying the -m/--metadata or -i/--sdf optiosn overrides the 
       -d/--duration setting and the entire file is reduced.
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0.0
	config['average'] = 1.0
	config['LFFT'] = 4096
	config['freq1'] = 0
	config['freq2'] = 0
	config['maxFrames'] = 28000
	config['window'] = fxc.noWindow
	config['duration'] = 0.0
	config['verbose'] = True
	config['clip'] = 0
	config['estimate'] = False
	config['metadata'] = None
	config['sdf'] = None
	config['site'] = 'lwa1'
	config['force'] = False
	config['linear'] = True
	config['countSats'] = True
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqtbnl:s:a:d:c:em:i:vfkw", ["help", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "skip=", "average=", "duration=", "freq1=", "freq2=", "clip-level=", "estimate-clip", "metadata=", "sdf=", "lwasv", "force", "stokes", "without-sats"])
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
		elif opt in ('-t', '--bartlett'):
			config['window'] = numpy.bartlett
		elif opt in ('-b', '--blackman'):
			config['window'] = numpy.blackman
		elif opt in ('-n', '--hanning'):
			config['window'] = numpy.hanning
		elif opt in ('-l', '--fft-length'):
			config['LFFT'] = int(value)
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-a', '--average'):
			config['average'] = float(value)
		elif opt in ('-d', '--duration'):
			config['duration'] = float(value)
		elif opt in ('-c', '--clip-level'):
			config['clip'] = int(value)
		elif opt in ('-e', '--estimate-clip'):
			config['estimate'] = True
		elif opt in ('-m', '--metadata'):
			config['metadata'] = value
		elif opt in ('-i', '--sdf'):
			config['sdf'] = value
		elif opt in ('-v', '--lwasv'):
			config['site'] = 'lwasv'
		elif opt in ('-f', '--force'):
			config['force'] = True
		elif opt in ('-k', '--stokes'):
			config['linear'] = False
		elif opt in ('-w', '--without-sats'):
			config['countSats'] = False
		else:
			assert False
			
	# Add in arguments
	config['args'] = args
	
	# Return configuration
	return config


def bestFreqUnits(freq):
	"""Given a numpy array of frequencies in Hz, return a new array with the
	frequencies in the best units possible (kHz, MHz, etc.)."""
	
	# Figure out how large the data are
	scale = int(math.log10(freq.max()))
	if scale >= 9:
		divis = 1e9
		units = 'GHz'
	elif scale >= 6:
		divis = 1e6
		units = 'MHz'
	elif scale >= 3:
		divis = 1e3
		units = 'kHz'
	else:
		divis = 1
		units = 'Hz'
		
	# Convert the frequency
	newFreq = freq / divis
	
	# Return units and freq
	return (newFreq, units)


def processDataBatchLinear(idf, antennas, tStart, duration, sampleRate, config, dataSets, obsID=1, clip1=0, clip2=0):
	"""
	Process a chunk of data in a raw DRX file into linear polarization 
	products and add the contents to an HDF5 file.
	"""
	
	# Length of the FFT
	LFFT = config['LFFT']
	
	# Find the start of the observation
	print 'Looking for #%i at %s with sample rate %.1f Hz...' % (obsID, tStart, sampleRate)
	idf.reset()
	
	t0 = idf.getInfo('tStart')
	tDiff = tStart - datetime.utcfromtimestamp(t0)
	offset = idf.offset( tDiff.total_seconds() + 1 )
	t0 = idf.getInfo('tStart')
	srate = idf.getInfo('sampleRate')
	
	print '... Found #%i at %s with sample rate %.1f Hz' % (obsID, datetime.utcfromtimestamp(t0), srate)
	tDiff = datetime.utcfromtimestamp(t0) - tStart
	duration = duration - max([0, tDiff.total_seconds()])
	
	# Number of remaining chunks (and the correction to the number of
	# frames to read in).
	nChunks = int(round(duration / config['average']))
	if nChunks == 0:
		nChunks = 1
		
	# Date & Central Frequency
	beginDate = ephem.Date(unix_to_utcjd(t0) - DJD_OFFSET)
	centralFreq1 = idf.getInfo('freq1')
	centralFreq2 = idf.getInfo('freq2')
	freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d=1/srate))
	if float(fxc.__version__) < 0.8:
		freq = freq[1:]
		
	dataSets['obs%i-freq1' % obsID][:] = freq + centralFreq1
	dataSets['obs%i-freq2' % obsID][:] = freq + centralFreq2
	
	obs = dataSets['obs%i' % obsID]
	obs.attrs['tInt'] = config['average']
	obs.attrs['tInt_Unit'] = 's'
	obs.attrs['LFFT'] = LFFT
	obs.attrs['nChan'] = LFFT-1 if float(fxc.__version__) < 0.8 else LFFT
	obs.attrs['RBW'] = freq[1]-freq[0]
	obs.attrs['RBW_Units'] = 'Hz'
	
	dataProducts = ['XX', 'YY']
	done = False
	for i in xrange(nChunks):
		# Inner loop that actually reads the frames into the data array
		print "Working on chunk %i, %i chunks remaning" % (i+1, nChunks-i-1)
		print "Working on %.1f ms of data" % (config['average']*1000.0,)
		
		tInt, cTime, data = idf.read(config['average'])
		if i == 0:
			print "Actual integration time is %.1f ms" % (tInt*1000.0,)
			
		# Save out some easy stuff
		dataSets['obs%i-time' % obsID][i] = cTime
		
		if config['countSats']:
			sats = ((data.real**2 + data.imag**2) >= 49).sum(axis=1)
			dataSets['obs%i-Saturation1' % obsID][i,:] = sats[0:2]
			dataSets['obs%i-Saturation2' % obsID][i,:] = sats[2:4]
		else:
			dataSets['obs%i-Saturation1' % obsID][i,:] = -1
			dataSets['obs%i-Saturation2' % obsID][i,:] = -1
			
		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		if clip1 == clip2:
			freq, tempSpec1 = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=clip1)
			
			l = 0
			for t in (1,2):
				for p in dataProducts:
					dataSets['obs%i-%s%i' % (obsID, p, t)][i,:] = tempSpec1[l,:]
					l += 1
					
		else:
			freq, tempSpec1 = fxc.SpecMaster(data[:2,:], LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=clip1)
			freq, tempSpec2 = fxc.SpecMaster(data[2:,:], LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=clip2)
			
			for l,p in enumerate(dataProducts):
				dataSets['obs%i-%s%i' % (obsID, p, 1)][i,:] = tempSpec1[l,:]
				dataSets['obs%i-%s%i' % (obsID, p, 2)][i,:] = tempSpec2[l,:]
				
		# We don't really need the data array anymore, so delete it
		del(data)
		
		# Are we done yet?
		if done:
			break
			
	return True


def processDataBatchStokes(idf, antennas, tStart, duration, sampleRate, config, dataSets, obsID=1, clip1=0, clip2=0):
	"""
	Process a chunk of data in a raw DRX file into Stokes parameters and 
	add the contents to an HDF5 file.
	"""
	
	# Length of the FFT
	LFFT = config['LFFT']
	
	# Find the start of the observation
	t0 = idf.getInfo('tStart')
	
	print 'Looking for #%i at %s with sample rate %.1f Hz...' % (obsID, tStart, sampleRate)
	idf.reset()
	
	t0 = idf.getInfo('tStart')
	tDiff = tStart - datetime.utcfromtimestamp(t0)
	offset = idf.offset( tDiff.total_seconds() + 1 )
	t0 = idf.getInfo('tStart')
	srate = idf.getInfo('sampleRate')
	
	print '... Found #%i at %s with sample rate %.1f Hz' % (obsID, datetime.utcfromtimestamp(t0), srate)
	tDiff = datetime.utcfromtimestamp(t0) - tStart
	duration = duration - max([0, tDiff.total_seconds()])
		
	# Number of remaining chunks (and the correction to the number of
	# frames to read in).
	nChunks = int(round(duration / config['average']))
	if nChunks == 0:
		nChunks = 1
		
	# Number of remaining chunks (and the correction to the number of
	# frames to read in).
	nChunks = int(round(duration / config['average']))
	if nChunks == 0:
		nChunks = 1
		
	# Date & Central Frequency
	beginDate = ephem.Date(unix_to_utcjd(t0) - DJD_OFFSET)
	centralFreq1 = idf.getInfo('freq1')
	centralFreq2 = idf.getInfo('freq2')
	freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d=1/srate))
	if float(fxc.__version__) < 0.8:
		freq = freq[1:]
		
	dataSets['obs%i-freq1' % obsID][:] = freq + centralFreq1
	dataSets['obs%i-freq2' % obsID][:] = freq + centralFreq2
	
	obs = dataSets['obs%i' % obsID]
	obs.attrs['tInt'] = config['average']
	obs.attrs['tInt_Unit'] = 's'
	obs.attrs['LFFT'] = LFFT
	obs.attrs['nChan'] = LFFT-1 if float(fxc.__version__) < 0.8 else LFFT
	obs.attrs['RBW'] = freq[1]-freq[0]
	obs.attrs['RBW_Units'] = 'Hz'
	
	dataProducts = ['I', 'Q', 'U', 'V']
	done = False
	for i in xrange(nChunks):
		# Inner loop that actually reads the frames into the data array
		print "Working on chunk %i, %i chunks remaning" % (i+1, nChunks-i-1)
		print "Working on %.1f ms of data" % (config['average']*1000.0,)
		
		tInt, cTime, data = idf.read(config['average'])
		if i == 0:
			print "Actual integration time is %.1f ms" % (tInt*1000.0,)
			
		# Save out some easy stuff
		dataSets['obs%i-time' % obsID][i] = cTime
		
		if config['countSats']:
			sats = ((data.real**2 + data.imag**2) >= 49).sum(axis=1)
			dataSets['obs%i-Saturation1' % obsID][i,:] = sats[0:2]
			dataSets['obs%i-Saturation2' % obsID][i,:] = sats[2:4]
		else:
			dataSets['obs%i-Saturation1' % obsID][i,:] = -1
			dataSets['obs%i-Saturation2' % obsID][i,:] = -1
			
		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		if clip1 == clip2:
			freq, tempSpec1 = fxc.StokesMaster(data, antennas, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=clip1)
			
			for t in (1,2):
				for l,p in enumerate(dataProducts):
					dataSets['obs%i-%s%i' % (obsID, p, t)][i,:] = tempSpec1[l,t-1,:]
					
		else:
			freq, tempSpec1 = fxc.StokesMaster(data[:2,:], antennas[:2], LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=clip1)
			freq, tempSpec2 = fxc.StokesMaster(data[2:,:], antennas[2:], LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=clip2)
			
			for l,p in enumerate(dataProducts):
				dataSets['obs%i-%s%i' % (obsID, p, 1)][i,:] = tempSpec1[l,0,:]
				dataSets['obs%i-%s%i' % (obsID, p, 2)][i,:] = tempSpec2[l,0,:]
				
		# We don't really need the data array anymore, so delete it
		del(data)
		
		# Are we done yet?
		if done:
			break
			
	return True


def main(args):
	# Parse command line options
	config = parseOptions(args)
	
	# Length of the FFT
	LFFT = config['LFFT']
	
	# Open the file and find good data (not spectrometer data)
	filename = config['args'][0]
	fh = open(filename, "rb")
	
	try:
		for i in xrange(5):
			junkFrame = drspec.readFrame(fh)
		raise RuntimeError("ERROR: '%s' appears to be a DR spectrometer file, not a raw DRX file" % filename)
	except errors.syncError:
		fh.seek(0)
		
	# Good, we seem to have a real DRX file, switch over to the LDP interface
	fh.close()
	idf = LWA1DataFile(filename, ignoreTimeTagErrors=False)

	# Offset into the file
	offset = idf.offset(config['offset'])
	
	# Metadata
	nFramesFile = idf.getInfo('nFrames')
	beam = idf.getInfo('beam')
	srate = idf.getInfo('sampleRate')
	beampols = idf.getInfo('beampols')
	beams = max([1, beampols / 4])
	
	# Number of frames to integrate over
	nFramesAvg = int(config['average'] * srate / 4096 * beampols)
	nFramesAvg = int(1.0 * nFramesAvg / beampols*4096/float(LFFT))*LFFT/4096*beampols
	config['average'] = 1.0 * nFramesAvg / beampols * 4096 / srate
	maxFrames = nFramesAvg
	
	# Number of remaining chunks (and the correction to the number of
	# frames to read in).
	if config['metadata'] is not None:
		config['duration'] = 0
	if config['duration'] == 0:
		config['duration'] = 1.0 * nFramesFile / beampols * 4096 / srate
		config['duration'] -= config['offset']
	else:
		config['duration'] = int(round(config['duration'] * srate * beampols / 4096) / beampols * 4096 / srate)
	nChunks = int(round(config['duration'] / config['average']))
	if nChunks == 0:
		nChunks = 1
	nFrames = nFramesAvg*nChunks
	
	# Date & Central Frequency
	t1  = idf.getInfo('tStart')
	beginDate = ephem.Date(unix_to_utcjd(t1) - DJD_OFFSET)
	centralFreq1 = idf.getInfo('freq1')
	centralFreq2 = idf.getInfo('freq2')
	config['freq1'] = centralFreq1
	config['freq2'] = centralFreq2

	# File summary
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % str(beginDate)
	print "Beams: %i" % beams
	print "Tune/Pols: %i" % beampols
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Integration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFramesAvg, nFramesAvg / beampols)
	print "Duration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average']*nChunks, nFrames, nFrames / beampols)
	print "Chunks: %i" % nChunks
	print " "
	
	# Estimate clip level (if needed)
	if config['estimate']:
		estimate = idf.estimateLevels(fh, Sigma=5.0)
		clip1 = (estiamte[0] + estimate[1]) / 2.0
		clip2 = (estiamte[2] + estimate[3]) / 2.0
	else:
		clip1 = config['clip']
		clip2 = config['clip']
		
	# Make the pseudo-antennas for Stokes calculation
	antennas = []
	for i in xrange(4):
		if i / 2 == 0:
			newAnt = stations.Antenna(1)
		else:
			newAnt = stations.Antenna(2)
			
		if i % 2 == 0:
			newAnt.pol = 0
		else:
			newAnt.pol = 1
			
		antennas.append(newAnt)
		
	# Setup the output file
	outname = os.path.split(filename)[1]
	outname = os.path.splitext(outname)[0]
	outname = '%s-waterfall.hdf5' % outname
	
	if os.path.exists(outname):
		if not config['force']:
			yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
		else:
			yn = 'y'
			
		if yn not in ('n', 'N'):
			os.unlink(outname)
		else:
			raise RuntimeError("Output file '%s' already exists" % outname)
			
	f = hdfData.createNewFile(outname)
	
	# Look at the metadata and come up with a list of observations.  If 
	# there are no metadata, create a single "observation" that covers the
	# whole file.
	obsList = {}
	if config['metadata'] is not None:
		try:
			project = metabundle.getSessionDefinition(config['metadata'])
		except Exception as e:
			if adpReady:
				project = metabundleADP.getSessionDefinition(config['metadata'])
			else:
				raise e
				
		sdfBeam  = project.sessions[0].drxBeam
		spcSetup = project.sessions[0].spcSetup
		if sdfBeam != beam:
			raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
			
		for i,obs in enumerate(project.sessions[0].observations):
			sdfStart = mcs.mjdmpm2datetime(obs.mjd, obs.mpm)
			sdfStop  = mcs.mjdmpm2datetime(obs.mjd, obs.mpm + obs.dur)
			obsDur   = obs.dur/1000.0
			obsSR    = drx.filterCodes[obs.filter]
			
			obsList[i+1] = (sdfStart, sdfStop, obsDur, obsSR)
			
		print "Observations:"
		for i in sorted(obsList.keys()):
			obs = obsList[i]
			print " #%i: %s to %s (%.3f s) at %.3f MHz" % (i, obs[0], obs[1], obs[2], obs[3]/1e6)
		print " "
			
		hdfData.fillFromMetabundle(f, config['metadata'])
		
	elif config['sdf'] is not None:
		try:
			project = sdf.parseSDF(config['sdf'])
		except Exception as e:
			if adpReady:
				project = sdfADP.parseSDF(config['sdf'])
			else:
				raise e
				
		sdfBeam  = project.sessions[0].drxBeam
		spcSetup = project.sessions[0].spcSetup
		if sdfBeam != beam:
			raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
			
		for i,obs in enumerate(project.sessions[0].observations):
			sdfStart = mcs.mjdmpm2datetime(obs.mjd, obs.mpm)
			sdfStop  = mcs.mjdmpm2datetime(obs.mjd, obs.mpm + obs.dur)
			obsChunks = int(numpy.ceil(obs.dur/1000.0 * drx.filterCodes[obs.filter] / (spcSetup[0]*spcSetup[1])))
			
			obsList[i+1] = (sdfStart, sdfStop, obsChunks)
			
		hdfData.fillFromSDF(f, config['sdf'], station=config['site'])
		
	else:
		obsList[1] = (datetime.utcfromtimestamp(t1), datetime(2222,12,31,23,59,59), config['duration'], srate)
		
		hdfData.fillMinimum(f, 1, beam, srate, station=config['site'])
		
	if config['linear']:
		dataProducts = ['XX', 'YY']
	else:
		dataProducts = ['I', 'Q', 'U', 'V']
		
	for o in sorted(obsList.keys()):
		for t in (1,2):
			hdfData.createDataSets(f, o, t, numpy.arange(LFFT-1 if float(fxc.__version__) < 0.8 else LFFT, dtype=numpy.float32), int(round(obsList[o][2]/config['average'])), dataProducts)
			
	f.attrs['FileGenerator'] = 'hdfWaterfall.py'
	f.attrs['InputData'] = os.path.basename(filename)
	
	# Create the various HDF group holders
	ds = {}
	for o in sorted(obsList.keys()):
		obs = hdfData.getObservationSet(f, o)
		
		ds['obs%i' % o] = obs
		ds['obs%i-time' % o] = obs.create_dataset('time', (int(round(obsList[o][2]/config['average'])),), 'f8')
		
		for t in (1,2):
			ds['obs%i-freq%i' % (o, t)] = hdfData.getDataSet(f, o, t, 'freq')
			for p in dataProducts:
				ds["obs%i-%s%i" % (o, p, t)] = hdfData.getDataSet(f, o, t, p)
			ds['obs%i-Saturation%i' % (o, t)] = hdfData.getDataSet(f, o, t, 'Saturation')
			
	# Load in the correct analysis function
	if config['linear']:
		processDataBatch = processDataBatchLinear
	else:
		processDataBatch = processDataBatchStokes
		
	# Go!
	for o in sorted(obsList.keys()):
		try:
			processDataBatch(idf, antennas, obsList[o][0], obsList[o][2], obsList[o][3], config, ds, obsID=o, clip1=clip1, clip2=clip2)
		except RuntimeError, e:
			print "Observation #%i: %s, abandoning this observation" % (o, str(e))

	# Save the output to a HDF5 file
	f.close()
	
	# Close out the data file
	idf.close()


if __name__ == "__main__":
	main(sys.argv[1:])
