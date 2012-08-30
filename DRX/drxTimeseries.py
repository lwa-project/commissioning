#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a DRX file, plot the time series I and Q data as a function of time.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import math
import time
import numpy
import getopt

import lsl.reader.drx as drx
import lsl.reader.errors as errors

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """drxTimeseries.py - Read in DRX files and create a collection of 
timeseries (I/Q) plots.

Usage: drxTimeseries.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-p, --plot-range            Number of seconds of data to show in the I/Q plots
                            (default = 0.0001)
-i, --instanenous-power     Plot I*I + Q*Q instead of the raw samples
-m, --mark-frames           Mark the frame bounaries in time
-q, --quiet                 Run drxSpectra in silent mode
-o, --output                Output file name for time series image
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0.0
	config['average'] = 0.0001
	config['maxFrames'] = 19144
	config['output'] = None
	config['verbose'] = True
	config['doPower'] = False
	config['markFrames'] = False
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqo:s:p:im", ["help", "quiet", "output=", "skip=", "plot-range=", "instanenous-power", "mark-frames"])
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
		elif opt in ('-o', '--output'):
			config['output'] = value
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-p', '--plot-range'):
			config['average'] = float(value)
		elif opt in ('-i', '--instanenous-power'):
			config['doPower'] = True
		elif opt in ('-m', '--mark-frames'):
			config['markFrames'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	# Parse command line options
	config = parseOptions(args)
	
	fh = open(config['args'][0], "rb")
	nFramesFile = os.path.getsize(config['args'][0]) / drx.FrameSize
	
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
	
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol

	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(round(config['offset'] * srate / 4096 * beampols))
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
	offset = int(config['offset'] * srate / 4096 * beampols)

	# Make sure that the file chunk size contains is an intger multiple
	# of the beampols.
	maxFrames = int(config['maxFrames']/beampols)*beampols

	# Number of frames to integrate over
	toClip = False
	oldAverage = config['average']
	if config['average'] < 4096/srate:		
		toClip = True
		config['average'] = 4096/srate
	nFrames = int(config['average'] * srate / 4096 * beampols)
	nFrames = int(1.0 * nFrames / beampols) * beampols
	config['average'] = 1.0 * nFrames / beampols * 4096 / srate

	# Number of remaining chunks
	nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Plot time: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFrames, nFrames / beampols)
	print "Chunks: %i" % nChunks

	# Sanity check
	if nFrames > (nFramesFile - offset):
		raise RuntimeError("Requested integration time+offset is greater than file length")

	# Align the file handle so that the first frame read in the
	# main analysis loop is from tuning 1, polarization 0
	junkFrame = drx.readFrame(fh)
	b,t,p = junkFrame.parseID()
	while 2*(t-1)+p != 0:
		junkFrame = drx.readFrame(fh)
		b,t,p = junkFrame.parseID()
	fh.seek(-drx.FrameSize, 1)

	# Master loop over all of the file chuncks
	standMapper = []
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
		tt = numpy.zeros((beampols,framesWork/beampols), dtype=numpy.int64) - 1
		data = numpy.zeros((beampols,framesWork*4096/beampols), dtype=numpy.csingle)
		
		# Inner loop that actually reads the frames into the data array
		print "Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0)
		t0 = time.time()
		
		for j in xrange(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = drx.readFrame(fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/drx.FrameSize-1)
				continue
			except errors.numpyError:
				break
			
			beam,tune,pol = cFrame.parseID()
			aStand = 2*(tune-1) + pol
			
			tt[aStand, count[aStand]] = cFrame.data.timeTag
			if config['doPower']:
				data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = numpy.abs(cFrame.data.iq)**2
			else:
				data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
			
			# Update the counters so that we can average properly later on
			count[aStand] += 1
			
		# The plots:  This is setup for the current configuration of 20 beampols
		fig = plt.figure()
		figsX = int(round(math.sqrt(beampols)))
		figsY = beampols / figsX

		samples = int(oldAverage * srate)
		if toClip:
			print "Plotting only the first %i samples (%.3f ms) of data" % (samples, oldAverage*1000.0)
			
		sortedMapper = sorted(standMapper)
		for i in xrange(data.shape[0]):
			ax = fig.add_subplot(figsX,figsY,i+1)
			if config['doPower']:
				limits = (-10, 100)
				if toClip:
					ax.plot(config['offset'] + numpy.arange(0,samples)/srate, data[i,0:samples])
				else:
					ax.plot(config['offset'] + numpy.arange(0,data.shape[1])/srate, data[i,:])
			else:
				limits = (-8, 8)
				if toClip:
					ax.plot(config['offset'] + numpy.arange(0,samples)/srate, data[i,0:samples].real, label='I')
					ax.plot(config['offset'] + numpy.arange(0,samples)/srate, data[i,0:samples].imag, label='Q')
				else:
					ax.plot(config['offset'] + numpy.arange(0,data.shape[1])/srate, data[i,:].real, label='I')
					ax.plot(config['offset'] + numpy.arange(0,data.shape[1])/srate, data[i,:].imag, label='Q')
				ax.legend(loc=0)

			if config['markFrames']:
				for j in xrange(0, samples-4096, 4096):
					ax.vlines(float(j)/srate, limits[0], limits[1], color='k', label='%i' % tt[i,j/4096])

			ax.set_ylim(limits)
			ax.set_title('Beam %i, Tune. %i, Pol. %i' % (beam, i/2+1,i%2))
			ax.set_xlabel('Time [seconds]')
			if config['doPower']:
				ax.set_ylabel('I$^2$ + Q$^2$')
			else:
				ax.set_ylabel('Output Level')
		plt.show()

		# Save image if requested
		if config['output'] is not None:
			fig.savefig(config['output'])


if __name__ == "__main__":
	main(sys.argv[1:])
