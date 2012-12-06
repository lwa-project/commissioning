#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a DRX file, plot the instantaneous power as a function of time.

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
import ephem

from lsl import astro
import lsl.reader.drx as drx
import lsl.reader.errors as errors

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """drxPower.py - Read in DRX files and create a collection of 
timeseries power (I*I + Q*Q) plots.  These power measurements are saved to a NPZ 
file called <filename>-power.npz.  Also, create a series of polarization and 
tuning diagnostic plots.

Usage: drxPower.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-t, --trim-level            Trim level for power analysis with clipping
                            (default = 49)
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-a, --average               Number of seconds of data to average together for power
                            (default = 0.0002 = 0.2 ms)
-d, --duration              Number of seconds to calculate the waterfall for 
                            (default = 10)
-q, --quiet                 Run drxSpectra in silent mode
-o, --output                Output file name for average power image
-w, --write-npz             Write a NPZ file of the averaged power
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0.0
	config['average'] = 0.0002
	config['duration'] = 10.0
	config['maxFrames'] = 19144*3
	config['output'] = None
	config['verbose'] = True
	config['trimLevel'] = 49
	config['writeNPZ'] = False
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqt:o:s:a:d:w", ["help", "quiet", "trim-level=", "output=", "skip=", "average=", "duration=", "write-npz"])
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
		elif opt in ('-t', '--trim-level'):
			config['trim'] = int(value)
		elif opt in ('-o', '--output'):
			config['output'] = value
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-a', '--average'):
			config['average'] = float(value)
		elif opt in ('-d', '--duration'):
			config['duration'] = float(value)
		elif opt in ('-w', '--write-npz'):
			config['writeNPZ'] = True
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
	offset = int(round(config['offset'] * srate / 4096 * beampols))
	offset = int(1.0 * offset / beampols) * beampols

	# Make sure that the file chunk size contains is an intger multiple
	# of the beampols.
	maxFrames = int(round(config['average']*srate/4096))*beampols
	if maxFrames < beampols:
		maxFrames = beampols
	config['average'] = 1.0*maxFrames/beampols*4096/srate

	# Number of remaining chunks
	nChunks = int(round(config['duration'] / config['average']))
	nFrames = maxFrames * nChunks

	# Store the information about the first frame and convert the timetag to 
	# an ephem.Date object.
	prevTime = junkFrame.data.timeTag
	prevDate = ephem.Date(astro.unix_to_utcjd(junkFrame.getTime()) - astro.DJD_OFFSET)

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	print "Sample Rate: %i Hz" % srate
	print "Date of first frame: %i -> %s" % (prevTime, str(prevDate))
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Integration: %.4f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], maxFrames, maxFrames / beampols)
	print "Duration: %.4f s (%i frames; %i frames per beam/tune/pol)" % (config['average']*nChunks, nFrames, nFrames / beampols)
	print " "

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
	masterTimes = numpy.zeros((nChunks, 4), dtype=numpy.float64)
	masterData  = numpy.zeros((nChunks, 4), dtype=numpy.float32)
	masterData2 = numpy.zeros((nChunks, 4), dtype=numpy.float32)
	masterProcessed = [0, 0, 0, 0]
	masterRemoved = [0, 0, 0, 0]
	for i in range(nChunks):
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = nFrames - i*maxFrames
		if framesRemaining > maxFrames:
			framesWork = maxFrames
		else:
			framesWork = framesRemaining
		#print "Working on chunk %i, %i frames remaining" % (i, framesRemaining)
		
		count = {0:0, 1:0, 2:0, 3:0}
		data  = numpy.zeros((4,framesWork*4096/beampols), dtype=numpy.float32)
		data2 = numpy.zeros((4,framesWork*4096/beampols), dtype=numpy.float32)
		
		## Inner loop that actually reads the frames into the data array
		#print "Working on %.2f ms of data" % ((framesWork*4096/beampols/srate)*1000.0)
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
			
			if j < 4:
				masterTimes[i,aStand] = cFrame.getTime()

			try:
				framePower = numpy.abs(cFrame.data.iq)**2
				
				# Calculate the clipping
				mask = numpy.where( framePower <= config['trimLevel'], 1, 0 )
				
				data[ aStand,  count[aStand]*4096:(count[aStand]+1)*4096] = framePower
				data2[aStand,  count[aStand]*4096:(count[aStand]+1)*4096] = framePower*mask
				
				masterProcessed[aStand] += mask.size
				masterRemoved[aStand] += (mask.size - mask.sum())
			
				# Update the counters so that we can average properly later on
				count[aStand] += 1
			except ValueError:
				pass

		# Save the data
		masterData[i,:]  = data.sum(axis=1)
		masterData2[i,:] = data2.sum(axis=1)
		
	# Really save the data to a NPZ file
	if config['writeNPZ']:
		outfile = os.path.splitext(config['args'][0])[0]
		outfile = '%s-power.npz' % outfile
		numpy.savez(outfile, beam=beam, avgPowerFull=masterData, avgPowerTrim=masterData2, times=masterTimes, trimLevel=config['trimLevel'])

	# Report on the clipping
	print "Summary:"
	for i in xrange(4):
		print "  Tuning %i, Pol. %s:" % (i/2+1, 'X' if i%2 else 'Y')
		print "    Processed: %i samples" % masterProcessed[i]
		print "    Clipped:   %i samples" % masterRemoved[i]
		print "      -> %.1f%% blanked" % (100.0*masterRemoved[i]/masterProcessed[i],)

	# The plots:  This is setup for the current configuration of 20 beampols
	fig = plt.figure()
	figsX = int(round(math.sqrt(4)))
	figsY = 4 / figsX

	for i in xrange(masterData.shape[1]):
		ax = fig.add_subplot(figsX,figsY,i+1)
		ax.plot(numpy.arange(0, masterData.shape[0] )*config['average'], masterData[:,i],  label='Full')
		ax.plot(numpy.arange(0, masterData2.shape[0])*config['average'], masterData2[:,i], label='Trimmed')
		ax.set_ylim([0, masterData.max()])
		
		ax.set_title('Beam %i, Tune. %i, Pol. %i' % (beam, i/2+1, i%2))
		ax.set_xlabel('Time [seconds]')
		ax.set_ylabel('Output Power Level')

		ax.legend(loc=0)

	# Part 2, polarization stuff
	fig2 = plt.figure()
	ax = fig2.add_subplot(3, 2, 1)
	ax.plot(numpy.arange(0, masterData.shape[0])*config['average'], numpy.sqrt(masterData[:,0]**2 + masterData[:,1]**2))
	ax.set_title('$\\sqrt{X1^2 + Y1^2}$')
	ax.set_xlabel('Time [seconds]')

	ax = fig2.add_subplot(3, 2, 2)
	ax.plot(numpy.arange(0, masterData.shape[0])*config['average'], masterData[:,1] / masterData[:,0])
	ax.set_title('$Y1 / X1$')
	ax.set_xlabel('Time [seconds]')

	ax = fig2.add_subplot(3, 2, 3)
	ax.plot(numpy.arange(0, masterData.shape[0])*config['average'], numpy.sqrt(masterData[:,2]**2 + masterData[:,3]**2))
	ax.set_title('$\\sqrt{X2^2 + Y2^2}$')
	ax.set_xlabel('Time [seconds]')

	ax = fig2.add_subplot(3, 2, 4)
	ax.plot(numpy.arange(0, masterData.shape[0])*config['average'], masterData[:,3] / masterData[:,2])
	ax.set_title('$Y2 / X2$')
	ax.set_xlabel('Time [seconds]')

	ax = fig2.add_subplot(3, 2, 5)
	ax.plot(numpy.arange(0, masterData.shape[0])*config['average'], numpy.sqrt(masterData[:,2]**2 + masterData[:,3]**2) / numpy.sqrt(masterData[:,0]**2 + masterData[:,1]**2))
	ax.set_title('$\\sqrt{X2^2 + Y2^2} / \\sqrt{X1^2 + Y1^2}$')
	ax.set_xlabel('Time [seconds]')

	ax = fig2.add_subplot(3, 2, 6)
	ax.plot(numpy.arange(0, masterData.shape[0])*config['average'], (masterData[:,3]/masterData[:,2]) / (masterData[:,1]/masterData[:,0]))
	ax.set_title('$(Y2 / X2) / (Y1 / X1)$')
	ax.set_xlabel('Time [seconds]')

	plt.show()

	# Save image if requested
	if config['output'] is not None:
		fig.savefig(config['output'])


if __name__ == "__main__":
	main(sys.argv[1:])
