#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a TBF file, check the time tags.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import math
import ephem
import numpy

from lsl.common import stations
from lsl.reader import tbf
from lsl.reader import errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt

def main(args):
	filename = args[0]
	
	fh = open(filename, "rb")
	nFrames = os.path.getsize(filename) / tbf.FrameSize
	
	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbf.readFrame(fh)
	fh.seek(0)
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	
	# Figure out how many frames there are per observation and the number of
	# channels that are in the file
	nFramesPerObs = tbf.getFramesPerObs(fh)
	nChannels = tbf.getChannelCount(fh)
	nSamples = 7840
	
	# Figure out how many chunks we need to work with
	nChunks = nFrames / nFramesPerObs
	
	# Pre-load the channel mapper
	mapper = []
	for i in xrange(2*nFramesPerObs):
		cFrame = tbf.readFrame(fh)
		if cFrame.header.firstChan not in mapper:
			mapper.append( cFrame.header.firstChan )
	fh.seek(-2*nFramesPerObs*tbf.FrameSize, 1)
	mapper.sort()
	
	# File summary
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % str(beginDate)
	print "Frames per Observation: %i" % nFramesPerObs
	print "Channel Count: %i" % nChannels
	print "Frames: %i" % nFrames
	print "==="
	print "Chunks: %i" % nChunks
	
	# Master loop over all of the file chunks
	timeTags = numpy.zeros((nFramesPerObs, nChunks), dtype=numpy.int64) - 1
	for i in xrange(nChunks):
		# Inner loop that actually reads the frames into the data array
		for j in xrange(nFramesPerObs):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbf.readFrame(fh)
			except errors.eofError:
				break
			except errors.syncError:
				print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbf.FrameSize-1)
				continue
			if not cFrame.header.isTBF():
				continue
				
			firstChan = cFrame.header.firstChan
			
			# Figure out where to map the channel sequence to
			try:
				aStand = mapper.index(firstChan)
			except ValueError:
				mapper.append(firstChan)
				aStand = mapper.index(firstChan)
			
			if cFrame.header.frameCount % 10000 == 0:
				print "%4i -> %4i  %7i  %i" % (firstChan, aStand, cFrame.header.frameCount, cFrame.data.timeTag)
				
			# Actually load the data.  x pol goes into the even numbers, y pol into the 
			# odd numbers
			count = cFrame.header.frameCount - 1
			timeTags[aStand,   count] = cFrame.data.timeTag
			
	# Check for missing frames
	missing = numpy.where( timeTags < 0 )
	if len(missing) != 0:
		print "Found %i missing frames.  Missing data from:" % len(missing[0])
		for i,f in zip(missing[0], missing[1]):
			print "  channel set %4i @ frame %5i" % (mapper[i], f+1)
			
	# Check time tags to make sure every ant/pol as the same time as each frame
	for f in xrange(timeTags.shape[1]):
		## For each frame count value, get the median time tag and use this for comparison.
		## If things are really bad, we will get a lot of errors.
		frameTime = numpy.median( timeTags[:,f] )

		## Compare all of the antpols at a particular frame count, ignoring the ones that
		## are missing.
		missing = numpy.where( (timeTags[:,f] != frameTime) & (timeTags[:,f]>=0) )[0]

		## Report any errors
		for m in missing:
			print "ERROR: t.t. %i @ frame %i != frame median of %i" % (timeTags[m,f], f+1, frameTime)
			print "       -> difference: %i" % (timeTags[m,f]-frameTime,)

	# Check time tags to make sure the times increment correctly between frames
	for i in xrange(timeTags.shape[0]):
		for f in xrange(1,timeTags.shape[1]):
			## Skip missing frames since they always fail
			if timeTags[i,f] < 0 or timeTags[i,f-1] < 0:
				continue

			## Compare the current time tag with previous and report an error if there
			## is a discrepancy between the two modulo the expected skip.
			if timeTags[i,f] > (timeTags[i,f-1] + nSamples):
				## Too far into the future
				print "ERROR: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (timeTags[i,f], f+1, timeTags[i,f-1], f)
				print "       -> difference: %i" % (timeTags[i,f]-timeTags[i,f-1],)
			elif timeTags[i,f] < (timeTags[i,f-1] + nSamples):
				## Not far enough into the future
				print "ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i + skip" % (timeTags[i,f], f+1, timeTags[i,f-1], f)
				print "       -> difference: %i" % (timeTags[i,f]-timeTags[i,f-1],)
			else:
				## Everything is good if we make it here
				pass


if __name__ == "__main__":
	main(sys.argv[1:])

