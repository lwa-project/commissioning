#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a TBN file, check for missing frames (or frames considered missing by the
ring buffer) and plot what the missing packet rate and what packets might be 
missing.  Rather than do this for a whole file, it is done for some small portion
of the file that is controlled by the -s/--skip and -a/--average flags.

$Rev$
$LastChangedBy$
$LastChagnedDate$
"""

import os
import sys
import math
import numpy
import ephem
import getopt

from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.correlator import fx as fxc
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """checkMissing.py - Read in a TBN file and find out what is missing or
appears to be missing.

Usage: checkMissing.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-v, --lwasv                 Use mapping from LWA-SV instead of LWA1
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-a, --average               Number of seconds of data to examine for frame loss
                            (default = 10)
-q, --quiet                 Run checkMissing in silent mode
-o, --output                Output file name for diagnostic image
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['site'] = 'lwa1'
	config['maxFrames'] = 200*520
	config['offset'] = 0.0
	config['average'] = 10.0
	config['output'] = None
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqvo:s:a:", ["help", "quiet", "lwasv", "output=", "skip=", "average="])
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
		elif opt in ('-v', '--lwasv'):
			config['site'] = 'lwasv'
			config['maxFrames'] = 200*512
		elif opt in ('-o', '--output'):
			config['output'] = value
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-a', '--average'):
			config['average'] = float(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def plotMissing(ax1, ax2, missingPackets, missingList, antpols):
	d = ax1.imshow(missingPackets, vmin=0, vmax=1, cmap=plt.cm.gray)
	ax1.set_xlabel('Missing Packets Set')
	ax1.set_ylabel('Digitizer-1')
	#fig.colorbar(d, ax=ax1)
					
	ax2.plot(numpy.array(missingList)/float(antpols))
	ax2.set_xlabel('Frame Count Relative to File Start')
	ax2.set_ylabel('Fraction Missing Packets')


def main(args):
	# Parse command line options
	config = parseOptions(args)
	
	# Set the station
	if config['site'] == 'lwa1':
		station = stations.lwa1
	elif config['site'] == 'lwasv':
		station = stations.lwasv
	else:
		raise RuntimeError("Unknown site name: %s" % config['site'])
	antennas = station.getAntennas()
	
	fh = open(config['args'][0], "rb")
	nFramesFile = os.path.getsize(config['args'][0]) / tbn.FrameSize
	srate = tbn.getSampleRate(fh)
	#antpols = tbn.getFramesPerObs(fh)
	antpols = len(antennas)

	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(config['offset'] * srate / 512 * antpols)
	offset = int(1.0 * offset / antpols) * antpols
	config['offset'] = 1.0 * offset / antpols * 512 / srate
	fh.seek(offset*tbn.FrameSize)

	# Number of frames to integrate over
	nFrames = int(config['average'] * srate / 512 * antpols)
	config['average'] = 1.0 * nFrames / antpols * 512 / srate

	# Number of remaining chunks
	nChunks = int(math.ceil(1.0*(nFrames)/config['maxFrames']))

	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbn.readFrame(fh)
	fh.seek(0)
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Integration: %.3f s (%i frames; %i frames per stand/pol)" % (config['average'], nFrames, nFrames / antpols)
	print "Chunks: %i" % nChunks

	# Sanity check
	if offset > nFramesFile:
		raise RuntimeError("Requested offset is greater than file length")
	if nFrames > (nFramesFile - offset):
		raise RuntimeError("Requested integration time+offset is greater than file length")

	# Create the FrameBuffer instance
	buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1])

	# Master loop over all of the file chunks
	masterCount = [0 for a in xrange(len(antennas))]
	
	# Missing packet control variables
	missingPackets = numpy.ones((antpols, 2048), dtype=numpy.int8)
	pc = 0
	missing = 0
	missingList = []
	
	# Figure
	fig = plt.figure()
	ax1 = fig.add_subplot(1, 2, 1)
	ax2 = fig.add_subplot(1, 2, 2)
	
	k = 0
	for i in xrange(nChunks):
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (config['maxFrames']),
		# only deal with that chunk
		framesRemaining = nFrames - k
		if framesRemaining > config['maxFrames']:
			framesWork = config['maxFrames']
		else:
			framesWork = framesRemaining
		
		count = [0 for a in xrange(len(antennas))]
		
		j = 0
		fillsWork = framesWork / antpols
		# Inner loop that actually reads the frames into the data array
		while j < fillsWork:
			try:
				cFrame = tbn.readFrame(fh)
				k = k + 1
			except errors.eofError:
				break
			except errors.syncError:
				print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FrameSize-1)
				continue
			
			#print cFrame.header.frameCount, cFrame.data.timeTag, cFrame.parseID()
			
			buffer.append(cFrame)
			cFrames = buffer.get()

			if cFrames is None:
				continue
			
			valid = reduce(lambda x,y: x+int(y.valid), cFrames, 0)
			print "Frame #%5i:  %.4f seconds with %i valid ant/pols%s" % (cFrames[0].header.frameCount, cFrames[0].getTime(), valid, '!' if valid != antpols else '')
			if valid != antpols:
				bad = []
				for cFrame in cFrames:
					if not cFrame.valid:
						bad.append(cFrame.parseID())
						
						missingPackets[2*(bad[-1][0]-1)+bad[-1][1], pc] = 0
				bad.sort()
				
				pc += 1
				if pc == missingPackets.shape[1]:
					plotMissing(ax1, ax2, missingPackets, missingList, antpols)
					plt.show()
					sys.exit(0)

				missing += (antpols-valid)
				total = (buffer.full + buffer.partial)*antpols
				#print j, valid, antpols-valid, cFrames[0].header.frameCount, 1.0*missing / total* 100, bad[0], bad[-1], buffer.dropped
				#print buffer.status()
				
				missingList.append( antpols - valid )
			else:
				total = (buffer.full + buffer.partial)*antpols
				missingList.append(0)
				
			times = numpy.array([f.data.timeTag for f in cFrames], dtype=numpy.int64)
			#print cFrames[0].header.frameCount, times.min(), times.max(), times.max()-times.min(), "%6.3f%%" % (1.0*missing/total*100,)
			for cFrame in cFrames:
				stand,pol = cFrame.header.parseID()
				
				# In the current configuration, stands start at 1 and go up to 260.  So, we
				# can use this little trick to populate the data array
				aStand = 2*(stand-1)+pol
				
				# Update the counters so that we can average properly later on
				count[aStand] = count[aStand] + 1
				masterCount[aStand] = masterCount[aStand] + 1
			
			j += 1
	
	# Empty the remaining portion of the buffer and integrate what's left
	for cFrames in buffer.flush():
		valid = reduce(lambda x,y: x+int(y.valid), cFrames, 0)
		print "Frame #%5i:  %.4f seconds with %i valid ant/pols" % (cFrames[0].header.frameCount, cFrames[0].getTime(), valid)
		if valid != antpols:
			bad = []
			for cFrame in cFrames:
				if not cFrame.valid:
					bad.append(cFrame.parseID())
					
					missingPackets[2*(bad[-1][0]-1)+bad[-1][1], pc] = 0
			bad.sort()
			
			pc += 1
			if pc == missingPackets.shape[1]:
				plotMissing(ax1, ax2, missingPackets, missingList, antpols)
				plt.show()
				sys.exit(0)

			missing += (antpols-valid)
			total = (buffer.full + buffer.partial)*antpols
			#print j, valid, antpols-valid, cFrames[0].header.frameCount, 1.0*missing / total* 100, bad[0], bad[-1], buffer.dropped
			#print buffer.status()
			
			missingList.append( antpols - valid )
		else:
			total = (buffer.full + buffer.partial)*antpols
			missingList.append(0)
		
		# Inner loop that actually reads the frames into the data array
		for cFrame in cFrames:
			stand,pol = cFrame.header.parseID()
			# In the current configuration, stands start at 1 and go up to 10.  So, we
			# can use this little trick to populate the data array
			aStand = 2*(stand-1)+pol
							
			# Update the counters so that we can average properly later on
			count[aStand] = count[aStand] + 1
			masterCount[aStand] = masterCount[aStand] + 1

		j += 1

	plotMissing(ax1, ax2, missingPackets, missingList, antpols)
	plt.show()
	sys.exit(0)


if __name__ == "__main__":
	main(sys.argv[1:])
