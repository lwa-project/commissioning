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

import lsl.reader.drx as drx
import lsl.reader.drspec as drspec
import lsl.reader.errors as errors
import lsl.statistics.robust as robust
import lsl.correlator.fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common import progress, stations
from lsl.common import mcs, metabundle

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
-e, --estimate-clip         Use robust statistics to estimate an approprite clip 
                            level (overrides the `-c` option)
-m, --metadata              Metadata tarball for additional information
-k, --stokes                Generate Stokes parameters instead of XX and YY
-w, --without-sats          Do not generate saturation counts

Note:  Specifying the -m/--metadata option overrides the -d/--duration setting 
       and the entire file is reduced.
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
	config['linear'] = True
	config['countSats'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqtbnl:s:a:d:c:em:kw", ["help", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "skip=", "average=", "duration=", "freq1=", "freq2=", "clip-level=", "estimate-clip", "metadata=", "stokes", "without-sats"])
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


def estimateClipLevel(fh, beampols):
	"""
	Read in a set of 100 frames and come up with the 4-sigma clip levels 
	for each tuning.  These clip levels are returned as a two-element 
	tuple.
	"""
	
	filePos = fh.tell()
		
	# Read in the first 100 frames for each tuning/polarization
	count = {0:0, 1:0, 2:0, 3:0}
	data = numpy.zeros((4, 4096*100), dtype=numpy.csingle)
	for i in xrange(beampols*100):
		try:
			cFrame = drx.readFrame(fh, Verbose=False)
		except errors.eofError:
			break
		except errors.syncError:
			continue
		
		beam,tune,pol = cFrame.parseID()
		aStand = 2*(tune-1) + pol
		
		data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
		count[aStand] +=  1
	
	# Go back to where we started
	fh.seek(filePos)
	
	# Compute the robust mean and standard deviation for I and Q for each
	# tuning/polarization
	meanI = []
	meanQ = []
	stdsI = []
	stdsQ = []
	for i in xrange(4):
		meanI.append( robust.mean(data[i,:].real) )
		meanQ.append( robust.mean(data[i,:].imag) )
		
		stdsI.append( robust.std(data[i,:].real) )
		stdsQ.append( robust.std(data[i,:].imag) )
	
	# Report
	print "Statistics:"
	for i in xrange(4):
		print " Mean %i: %.3f + %.3f j" % (i+1, meanI[i], meanQ[i])
		print " Std  %i: %.3f + %.3f j" % (i+1, stdsI[i], stdsQ[i])
	
	# Come up with the clip levels based on 4 sigma
	clip1 = (meanI[0] + meanI[1] + meanQ[0] + meanQ[1]) / 4.0
	clip2 = (meanI[2] + meanI[3] + meanQ[2] + meanQ[3]) / 4.0
	
	clip1 += 5*(stdsI[0] + stdsI[1] + stdsQ[0] + stdsQ[1]) / 4.0
	clip2 += 5*(stdsI[2] + stdsI[3] + stdsQ[2] + stdsQ[3]) / 4.0
	
	clip1 = int(round(clip1))
	clip2 = int(round(clip2))
	
	# Report again
	print "Clip Levels:"
	print " Tuning 1: %i" % clip1
	print " Tuning 2: %i" % clip2
	
	return clip1, clip2


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


def processDataBatchLinear(fh, antennas, tStart, duration, sampleRate, config, dataSets, obsID=1, clip1=0, clip2=0):
	"""
	Process a chunk of data in a raw DRX file into linear polarization 
	products and add the contents to an HDF5 file.
	"""
	
	# Length of the FFT
	LFFT = config['LFFT']
	
	# Find the start of the observation
	junkFrame = drx.readFrame(fh)
	srate = junkFrame.getSampleRate()
	t0 = junkFrame.getTime()
	fh.seek(-drx.FrameSize, 1)
	
	print 'Looking for #%i at %s with sample rate %.1f Hz...' % (obsID, tStart, sampleRate)
	while datetime.utcfromtimestamp(t0) < tStart or srate != sampleRate:
		junkFrame = drx.readFrame(fh)
		srate = junkFrame.getSampleRate()
		t0 = junkFrame.getTime()
	print '... Found #%i at %s with sample rate %.1f Hz' % (obsID, datetime.utcfromtimestamp(t0), srate)
	tDiff = datetime.utcfromtimestamp(t0) - tStart
	try:
		duration = duration - tDiff.total_seconds()
	except:
		duration = duration - (tDiff.seconds + tDiff.microseconds/1e6)
	
	beam,tune,pol = junkFrame.parseID()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol
	
	# Make sure that the file chunk size contains is an integer multiple
	# of the FFT length so that no data gets dropped.  This needs to
	# take into account the number of beampols in the data, the FFT length,
	# and the number of samples per frame.
	maxFrames = int(1.0*config['maxFrames']/beampols*4096/float(LFFT))*LFFT/4096*beampols
	
	# Number of frames to integrate over
	nFramesAvg = int(config['average'] * srate / 4096 * beampols)
	nFramesAvg = int(1.0 * nFramesAvg / beampols*4096/float(LFFT))*LFFT/4096*beampols
	config['average'] = 1.0 * nFramesAvg / beampols * 4096 / srate
	maxFrames = nFramesAvg

	# Number of remaining chunks (and the correction to the number of
	# frames to read in).
	nChunks = int(round(duration / config['average']))
	if nChunks == 0:
		nChunks = 1
	nFrames = nFramesAvg*nChunks
	
	# Date & Central Frequnecy
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	centralFreq1 = 0.0
	centralFreq2 = 0.0
	for i in xrange(4):
		junkFrame = drx.readFrame(fh)
		b,t,p = junkFrame.parseID()
		if p == 0 and t == 1:
			try:
				centralFreq1 = junkFrame.getCentralFreq()
			except AttributeError:
				from lsl.common.dp import fS
				centralFreq1 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
		elif p == 0 and t == 2:
			try:
				centralFreq2 = junkFrame.getCentralFreq()
			except AttributeError:
				from lsl.common.dp import fS
				centralFreq2 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
		else:
			pass
	fh.seek(-4*drx.FrameSize, 1)
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
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = nFrames - i*maxFrames
		if framesRemaining > maxFrames:
			framesWork = maxFrames
		else:
			framesWork = framesRemaining
		print "Working on chunk %i, %i frames remaining" % (i+1, framesRemaining)
		
		count = {0:0, 1:0, 2:0, 3:0}
		data = numpy.zeros((4,framesWork*4096/beampols), dtype=numpy.csingle)
		# If there are fewer frames than we need to fill an FFT, skip this chunk
		if data.shape[1] < LFFT:
			break
			
		# Inner loop that actually reads the frames into the data array
		print "Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0)
		
		for j in xrange(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = drx.readFrame(fh, Verbose=False)
			except errors.eofError:
				done = True
				break
			except errors.syncError:
				continue

			beam,tune,pol = cFrame.parseID()
			aStand = 2*(tune-1) + pol
			if j is 0:
				cTime = cFrame.getTime()
			
			try:
				data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
				count[aStand] +=  1
			except ValueError:
				raise RuntimeError("Invalid Shape")
				
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


def processDataBatchStokes(fh, antennas, tStart, duration, sampleRate, config, dataSets, obsID=1, clip1=0, clip2=0):
	"""
	Process a chunk of data in a raw DRX file into Stokes parameters and 
	add the contents to an HDF5 file.
	"""
	
	# Length of the FFT
	LFFT = config['LFFT']
	
	# Find the start of the observation
	junkFrame = drx.readFrame(fh)
	srate = junkFrame.getSampleRate()
	t0 = junkFrame.getTime()
	fh.seek(-drx.FrameSize, 1)
	
	print 'Looking for #%i at %s with sample rate %.1f Hz...' % (obsID, tStart, sampleRate)
	while datetime.utcfromtimestamp(t0) < tStart or srate != sampleRate:
		junkFrame = drx.readFrame(fh)
		srate = junkFrame.getSampleRate()
		t0 = junkFrame.getTime()
	print '... Found #%i at %s with sample rate %.1f Hz' % (obsID, datetime.utcfromtimestamp(t0), srate)
	tDiff = datetime.utcfromtimestamp(t0) - tStart
	try:
		duration = duration - tDiff.total_seconds()
	except:
		duration = duration - (tDiff.seconds + tDiff.microseconds/1e6)
	
	beam,tune,pol = junkFrame.parseID()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol
	
	# Make sure that the file chunk size contains is an integer multiple
	# of the FFT length so that no data gets dropped.  This needs to
	# take into account the number of beampols in the data, the FFT length,
	# and the number of samples per frame.
	maxFrames = int(1.0*config['maxFrames']/beampols*4096/float(LFFT))*LFFT/4096*beampols
	
	# Number of frames to integrate over
	nFramesAvg = int(config['average'] * srate / 4096 * beampols)
	nFramesAvg = int(1.0 * nFramesAvg / beampols*4096/float(LFFT))*LFFT/4096*beampols
	config['average'] = 1.0 * nFramesAvg / beampols * 4096 / srate
	maxFrames = nFramesAvg

	# Number of remaining chunks (and the correction to the number of
	# frames to read in).
	nChunks = int(round(duration / config['average']))
	if nChunks == 0:
		nChunks = 1
	nFrames = nFramesAvg*nChunks
	
	# Date & Central Frequnecy
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	centralFreq1 = 0.0
	centralFreq2 = 0.0
	for i in xrange(4):
		junkFrame = drx.readFrame(fh)
		b,t,p = junkFrame.parseID()
		if p == 0 and t == 1:
			try:
				centralFreq1 = junkFrame.getCentralFreq()
			except AttributeError:
				from lsl.common.dp import fS
				centralFreq1 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
		elif p == 0 and t == 2:
			try:
				centralFreq2 = junkFrame.getCentralFreq()
			except AttributeError:
				from lsl.common.dp import fS
				centralFreq2 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
		else:
			pass
	fh.seek(-4*drx.FrameSize, 1)
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
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = nFrames - i*maxFrames
		if framesRemaining > maxFrames:
			framesWork = maxFrames
		else:
			framesWork = framesRemaining
		print "Working on chunk %i, %i frames remaining" % (i+1, framesRemaining)
		
		count = {0:0, 1:0, 2:0, 3:0}
		data = numpy.zeros((4,framesWork*4096/beampols), dtype=numpy.csingle)
		# If there are fewer frames than we need to fill an FFT, skip this chunk
		if data.shape[1] < LFFT:
			break
			
		# Inner loop that actually reads the frames into the data array
		print "Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0)
		
		for j in xrange(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = drx.readFrame(fh, Verbose=False)
			except errors.eofError:
				done = True
				break
			except errors.syncError:
				continue

			beam,tune,pol = cFrame.parseID()
			aStand = 2*(tune-1) + pol
			if j is 0:
				cTime = cFrame.getTime()
			
			try:
				data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
				count[aStand] +=  1
			except ValueError:
				raise RuntimeError("Invalid Shape")

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

	# Open the file and find good data (not spectreomter data)
	filename = config['args'][0]
	fh = open(filename, "rb")
	nFramesFile = os.path.getsize(filename) / drx.FrameSize
	
	try:
		for i in xrange(5):
			junkFrame = drspec.readFrame(fh)
		raise RuntimeError("ERROR: '%s' appears to be a DR spectrometer file, not a raw DRX file" % filename)
	except errors.syncError:
		fh.seek(0)
		
	while True:
		try:
			junkFrame = drx.readFrame(fh)
			try:
				srate = junkFrame.getSampleRate()
				t0 = junkFrame.getTime()
				break
			except ZeroDivisionError:
				pass
		except errors.syncError:
			fh.seek(-drx.FrameSize+1, 1)
			
	fh.seek(-drx.FrameSize, 1)
	
	beam,tune,pol = junkFrame.parseID()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol

	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(config['offset'] * srate / 4096 * beampols)
	offset = int(1.0 * offset / beampols) * beampols
	fh.seek(offset*drx.FrameSize, 1)
	
	# Iterate on the offsets until we reach the right point in the file.  This
	# is needed to deal with files that start with only one tuning and/or a 
	# different sample rate.  
	while True:
		## Figure out where in the file we are and what the current tuning/sample 
		## rate is
		junkFrame = drx.readFrame(fh)
		srate = junkFrame.getSampleRate()
		t1 = junkFrame.getTime()
		tunepols = drx.getFramesPerObs(fh)
		tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
		beampols = tunepol
		fh.seek(-drx.FrameSize, 1)
		
		## See how far off the current frame is from the target
		tDiff = t1 - (t0 + config['offset'])
		
		## Half that to come up with a new seek parameter
		tCorr = -tDiff / 2.0
		cOffset = int(tCorr * srate / 4096 * beampols)
		cOffset = int(1.0 * cOffset / beampols) * beampols
		offset += cOffset
		
		## If the offset is zero, we are done.  Otherwise, apply the offset
		## and check the location in the file again/
		if cOffset is 0:
			break
		fh.seek(cOffset*drx.FrameSize, 1)
	
	# Update the offset actually used
	config['offset'] = t1 - t0
	offset = int(round(config['offset'] * srate / 4096 * beampols))
	offset = int(1.0 * offset / beampols) * beampols

	# Make sure that the file chunk size contains is an integer multiple
	# of the FFT length so that no data gets dropped.  This needs to
	# take into account the number of beampols in the data, the FFT length,
	# and the number of samples per frame.
	maxFrames = int(1.0*config['maxFrames']/beampols*4096/float(LFFT))*LFFT/4096*beampols

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
	else:
		config['duration'] = int(round(config['duration'] * srate * beampols / 4096) / beampols * 4096 / srate)
	nChunks = int(round(config['duration'] / config['average']))
	if nChunks == 0:
		nChunks = 1
	nFrames = nFramesAvg*nChunks
	
	# Date & Central Frequnecy
	t1  = junkFrame.getTime()
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	centralFreq1 = 0.0
	centralFreq2 = 0.0
	for i in xrange(4):
		junkFrame = drx.readFrame(fh)
		b,t,p = junkFrame.parseID()
		if p == 0 and t == 1:
			try:
				centralFreq1 = junkFrame.getCentralFreq()
			except AttributeError:
				from lsl.common.dp import fS
				centralFreq1 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
		elif p == 0 and t == 2:
			try:
				centralFreq2 = junkFrame.getCentralFreq()
			except AttributeError:
				from lsl.common.dp import fS
				centralFreq2 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
		else:
			pass
	fh.seek(-4*drx.FrameSize, 1)
	
	config['freq1'] = centralFreq1
	config['freq2'] = centralFreq2

	# File summary
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % str(beginDate)
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
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
		clip1, clip2 = estimateClipLevel(fh, beampols)
	else:
		clip1 = config['clip']
		clip2 = config['clip']
		
	# Make the pseudo-antennas for Stokes calculation
	antennas = []
	for i in xrange(beampols):
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
	f = hdfData.createNewFile(outname)
	
	# Look at the metadata and come up with a list of observations.  If 
	# there are no metadata, create a single "observation" that covers the
	# whole file.
	obsList = {}
	if config['metadata'] is not None:
		sdf = metabundle.getSessionDefinition(config['metadata'])
		
		sdfBeam  = sdf.sessions[0].drxBeam
		spcSetup = sdf.sessions[0].spcSetup
		if sdfBeam != beam:
			raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
			
		for i,obs in enumerate(sdf.sessions[0].observations):
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
	else:
		obsList[1] = (datetime.utcfromtimestamp(t1), datetime(2222,12,31,23,59,59), config['duration'], srate)
		
		hdfData.fillMinimum(f, 1, beam, srate)
		
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
			processDataBatch(fh, antennas, obsList[o][0], obsList[o][2], obsList[o][3], config, ds, obsID=o, clip1=clip1, clip2=clip2)
		except RuntimeError, e:
			print "Observation #%i: %s, abondonding this observation" % (o, str(e))

	# Save the output to a HDF5 file
	f.close()


if __name__ == "__main__":
	main(sys.argv[1:])
