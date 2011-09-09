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
import math
import numpy
import ephem
import getopt

import lsl.reader.drx as drx
import lsl.reader.errors as errors
import lsl.statistics.robust as robust
import lsl.correlator.fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """drxWaterfall.py - Read in DRX files and create a collection of 
time-averaged spectra.  These spectra are saved to a NPZ file called <filename>-waterfall.npz.

Usage: drxWaterfall.py [OPTIONS] file

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
                            (default = 10)
-q, --quiet                 Run drxSpectra in silent mode and do not show the plots
-l, --fft-length            Set FFT length (default = 4096)
-c, --clip-level            FFT blanking clipping level in counts (default = 0, 
                            0 disables)
-e, --estimate-clip         Use robust statistics to estimate an approprite clip 
                            level (overrides the `-c` option)
-o, --output                Output file name for waterfall image
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
	config['LFFT'] = 1024
	config['maxFrames'] = 28000
	config['window'] = fxc.noWindow
	config['output'] = None
	config['duration'] = 10.0
	config['verbose'] = True
	config['clip'] = 0
	config['estimate'] = False
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqtbnl:o:s:a:d:c:e", ["help", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "output=", "skip=", "average=", "duration=", "clip-level=", "estimate-clip"])
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
		elif opt in ('-o', '--output'):
			config['output'] = value
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


def main(args):
	# Parse command line options
	config = parseOptions(args)
	print config['clip']

	# Length of the FFT
	LFFT = config['LFFT']

	fh = open(config['args'][0], "rb")
	nFramesFile = os.path.getsize(config['args'][0]) / drx.FrameSize
	junkFrame = drx.readFrame(fh)
	fh.seek(0)
	
	srate = junkFrame.getSampleRate()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol

	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(config['offset'] * srate / 4096 * beampols)
	offset = int(1.0 * offset / beampols) * beampols
	config['offset'] = 1.0 * offset / beampols * 4096 / srate
	fh.seek(offset*drx.FrameSize)

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
	nChunks = int(round(config['duration'] / config['average']))
	if nChunks == 0:
		nChunks = 1
	nFrames = nFramesAvg*nChunks
	
	# Date
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Integration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFramesAvg, nFramesAvg / beampols)
	print "Duration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average']*nChunks, nFrames, nFrames / beampols)
	print "Chunks: %i" % nChunks

	# Sanity check
	if offset > nFramesFile:
		raise RuntimeError("Requested offset is greater than file length")
	if nFrames > (nFramesFile - offset):
		raise RuntimeError("Requested integration time+offset is greater than file length")
	
	# Estimate clip level (if needed)
	if config['estimate']:
		filePos = fh.tell()
		
		# Read in the first 100 frames for each tuning/polarization
		count = {0:0, 1:0, 2:0, 3:0}
		data = numpy.zeros((4, 4096*100), dtype=numpy.csingle)
		for i in xrange(4*100):
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
		
	else:
		clipLevel1 = config['clip']
		clipLevel2 = config['clip']

	# Master loop over all of the file chunks
	masterCount = {0:0, 1:0, 2:0, 3:0}
	masterWeight = numpy.zeros((nChunks, beampols, LFFT-1))
	masterSpectra = numpy.zeros((nChunks, beampols, LFFT-1))
	masterTimes = numpy.zeros(nChunks)
	for i in xrange(nChunks):
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = nFrames - i*maxFrames
		if framesRemaining > maxFrames:
			framesWork = maxFrames
		else:
			framesWork = framesRemaining
		print "Working on chunk %i, %i frames remaining" % (i, framesRemaining)
		
		count = {0:0, 1:0, 2:0, 3:0}
		data = numpy.zeros((beampols,framesWork*4096/beampols), dtype=numpy.csingle)
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
				break
			except errors.syncError:
				continue

			beam,tune,pol = cFrame.parseID()
			aStand = 2*(tune-1) + pol
			if j is 0:
				cTime = cFrame.getTime()
			
			data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
			# Update the counters so that we can average properly later on
			count[aStand] +=  1
			masterCount[aStand] += 1

		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		freq, tempSpec1 = fxc.SpecMaster(data[:2,:], LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=clip1)
		
		freq, tempSpec2 = fxc.SpecMaster(data[2:,:], LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=clip2)
		
		# Save the results to the various master arrays
		masterTimes[i] = cTime
		
		masterSpectra[i,0,:] = tempSpec1[0,:]
		masterSpectra[i,1,:] = tempSpec1[1,:]
		masterSpectra[i,2,:] = tempSpec2[0,:]
		masterSpectra[i,3,:] = tempSpec2[1,:]
		
		masterWeight[i,0,:] = count[0]
		masterWeight[i,1,:] = count[1]
		masterWeight[i,2,:] = count[2]
		masterWeight[i,3,:] = count[3]

		# We don't really need the data array anymore, so delete it
		del(data)

	# Now that we have read through all of the chunks, perform the final averaging by
	# dividing by all of the chunks
	outname = config['args'][0].replace('.dat', '-waterfall.npz')
	numpy.savez(outname, freq=freq, times=masterTimes, spec=masterSpectra, tInt=(maxFrames*4096/beampols/srate), standMapper=[4*(beam-1) + i for i in xrange(masterSpectra.shape[1])])
	spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )

	# The plots:  This is setup for the current configuration of 20 beampols
	fig = plt.figure()
	figsX = int(round(math.sqrt(beampols)))
	figsY = beampols / figsX
	
	# Put the frequencies in the best units possible
	freq, units = bestFreqUnits(freq)

	for i in xrange(masterSpectra.shape[1]):
		ax = fig.add_subplot(figsX,figsY,i+1)
		currSpectra = numpy.squeeze( numpy.log10(masterSpectra[:,i,:])*10.0 )
		currSpectra = numpy.where( numpy.isfinite(currSpectra), currSpectra, -10)
		
		#ax.plot(freq, currSpectra, label='%i (avg)' % (i+1))
		ax.imshow(currSpectra, extent=(freq.min(), freq.max(), 0, config['average']*nChunks))
		print currSpectra.min(), currSpectra.max()

		ax.axis('auto')
		ax.set_title('Beam %i, Tune. %i, Pol. %i' % (beam, i+1, i%2))
		ax.set_xlabel('Frequency Offset [%s]' % units)
		ax.set_ylabel('Time [s]')
		ax.set_xlim([freq.min(), freq.max()])

	print "RBW: %.4f %s" % ((freq[1]-freq[0]), units)
	plt.subplots_adjust(hspace=0.35, wspace=0.30)
	if config['verbose']:
		plt.show()

	# Save spectra image if requested
	if config['output'] is not None:
		fig.savefig(config['output'])


if __name__ == "__main__":
	main(sys.argv[1:])
