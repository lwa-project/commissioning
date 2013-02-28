#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a DRX file, plot the time averaged spectra for each beam output.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import math
import ephem
import numpy
import getopt

import lsl.reader.drx as drx
import lsl.reader.errors as errors
import lsl.correlator.fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """drxSpectra.py - Read in DRX files and create a collection of 
time-averaged spectra.

Usage: drxSpectra.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-t, --bartlett              Apply a Bartlett window to the data
-b, --blackman              Apply a Blackman window to the data
-n, --hanning               Apply a Hanning window to the data
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-a, --average               Number of seconds of data to average for spectra 
                            (default = 10)
-q, --quiet                 Run drxSpectra in silent mode
-l, --fft-length            Set FFT length (default = 4096)
-c, --clip-level            FFT blanking clipping level in counts (default = 0, 
                            0 disables)
-o, --output                Output file name for spectra image
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0.0
	config['average'] = 10.0
	config['LFFT'] = 4096
	config['freq1'] = 0
	config['freq2'] = 0
	config['maxFrames'] = 19144*4
	config['window'] = fxc.noWindow
	config['output'] = None
	config['displayChunks'] = True
	config['verbose'] = True
	config['clip'] = 0
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqtbnl:o:s:a:c:", ["help", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "output=", "skip=", "average=", "clip-level="])
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
		elif opt in ('-c', '--clip-level'):
			config['clip'] = int(value)
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
	scale = int(math.log10((freq - freq.mean()).max()))
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

	# Length of the FFT
	LFFT = config['LFFT']

	fh = open(config['args'][0], "rb")
	nFramesFile = os.path.getsize(config['args'][0]) / drx.FrameSize
	
	while True:
		junkFrame = drx.readFrame(fh)
		try:
			srate = junkFrame.getSampleRate()
			t0 = junkFrame.getTime()
			break
		except ZeroDivisionError:
			pass
	fh.seek(-drx.FrameSize, 1)
	
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol

	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(config['offset'] * srate / 4096 * beampols)
	offset = int(1.0 * offset / beampols) * beampols
	fh.seek(offset*drx.FrameSize)
	
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

	# Make sure that the file chunk size contains is an intger multiple
	# of the FFT length so that no data gets dropped.  This needs to
	# take into account the number of beampols in the data, the FFT length,
	# and the number of samples per frame.
	maxFrames = int(1.0*config['maxFrames']/beampols*4096/float(LFFT))*LFFT/4096*beampols

	# Number of frames to integrate over
	nFrames = int(config['average'] * srate / 4096 * beampols)
	nFrames = int(1.0 * nFrames / beampols*4096/float(LFFT))*LFFT/4096*beampols
	config['average'] = 1.0 * nFrames / beampols * 4096 / srate

	# Number of remaining chunks
	nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))
	
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
	
	config['freq1'] = centralFreq1
	config['freq2'] = centralFreq2

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Integration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFrames, nFrames / beampols)
	print "Chunk Time: %.3f s (%i frames; %i frames per beam/tune/pol)" % (4096.0*maxFrames/beampols/srate, maxFrames, maxFrames/beampols)
	print "Chunks: %i" % nChunks

	# Sanity check
	if nFrames > (nFramesFile - offset):
		raise RuntimeError("Requestion integration time+offset is greater than file length")

	# Master loop over all of the file chuncks
	masterWeight = numpy.zeros((nChunks, 4, LFFT-1))
	masterSpectra = numpy.zeros((nChunks, 4, LFFT-1))
	for i in range(nChunks):
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
		data = numpy.zeros((4,framesWork*4096/beampols), dtype=numpy.csingle)
		# If there are fewer frames than we need to fill an FFT, skip this chunk
		if data.shape[1] < LFFT:
			break

		# Inner loop that actually reads the frames into the data array
		print "Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0)

		for j in range(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = drx.readFrame(fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				continue
			except errors.numpyError:
				break

			beam,tune,pol = cFrame.parseID()
			aStand = 2*(tune-1) + pol

			data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
			count[aStand] += 1

		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		freq1, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, CentralFreq=config['freq1'], ClipLevel=config['clip'])
		for stand in count.keys():
			masterSpectra[i,stand,:] = tempSpec[stand,:]
			masterWeight[i,stand,:] = int(count[stand] * 4096 / LFFT)

		# We don't really need the data array anymore, so delete it
		del(data)

	# Create the frequency array for the second tuning
	freq2 = freq1 - config['freq1'] + config['freq2']

	# Now that we have read through all of the chunks, peform the final averaging by
	# dividing by all of the chunks
	spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )

	# The plots:  This is setup for the current configuration of 20 beampols
	fig = plt.figure()
	figsX = int(round(math.sqrt(4)))
	figsY = 4 / figsX
	# Put the freqencies in the best units possible
	freq1, units1 = bestFreqUnits(freq1)
	freq2, units2 = bestFreqUnits(freq2)

	for i in xrange(masterSpectra.shape[1]):
		if i/2+1 == 1:
			freq = freq1
			units = units1
		else:
			freq = freq2
			units = units2

		ax = fig.add_subplot(figsX,figsY,i+1)
		currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
		ax.plot(freq, currSpectra, label='%i (avg)' % (i+1))

		ax.set_title('Beam %i, Tune. %i, Pol. %i' % (beam, i/2+1, i%2))
		if freq.min() < 0:
			ax.set_xlabel('Frequency Offset [%s]' % units)
		else:
			ax.set_xlabel('Frequency [%s]' % units)
		ax.set_ylabel('P.S.D. [dB/RBW]')
		ax.set_xlim([freq.min(), freq.max()])
		ax.legend(loc=0)

	print "RBW 1: %.4f %s" % ((freq1[1]-freq1[0]), units1)
	print "RBW 2: %.4f %s" % ((freq2[1]-freq2[0]), units2)
	plt.subplots_adjust(hspace=0.35, wspace=0.30)
	plt.show()
	
	outfile = os.path.split(config['args'][0])[1]
	outfile = os.path.splitext(outfile)[0]
	outfile = '%s.npz' % outfile
	numpy.savez(outfile, freq1=freq1, freq2=freq2, units1=units1, units2=units2, spec=spec, standMapper=[4*(beam-1) + i for i in xrange(masterSpectra.shape[1])])

	# Save spectra image if requested
	if config['output'] is not None:
		fig.savefig(config['output'])

if __name__ == "__main__":
	main(sys.argv[1:])
