#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a TBW file, check the time tags.

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
from lsl.reader import tbw
from lsl.reader import errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt

def main(args):
	filename = args[0]

	# Set the station
	station = stations.lwa1
	antennas = station.getAntennas()

	fh = open(filename, "rb")
	nFrames = os.path.getsize(filename) / tbw.FrameSize
	dataBits = tbw.getDataBits(fh)
	# The number of ant/pols in the file is hard coded because I cannot figure out 
	# a way to get this number in a systematic fashion
	maxFrames = 30000*260
	antpols = len(antennas)
	nChunks = int(math.ceil(1.0*nFrames/maxFrames))
	if dataBits == 12:
		nSamples = 400
	else:
		nSamples = 1200

	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbw.readFrame(fh)
	fh.seek(0)
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)

	# File summary
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Length: %i-bit" % dataBits
	print "Frames: %i" % nFrames
	print "Chunks: %i" % nChunks
	print "==="

	nChunks = 1

	# Skip over any non-TBW frames at the beginning of the file
	i = 0
	junkFrame = tbw.readFrame(fh)
	while not junkFrame.header.isTBW():
		junkFrame = tbw.readFrame(fh)
		i += 1
	fh.seek(-tbw.FrameSize, 1)
	print "Skipped %i non-TBW frames at the beginning of the file" % i

	# Master loop over all of the file chunks
	timeTags = numpy.zeros((antpols, 30000), dtype=numpy.int64) - 1
	for i in range(nChunks):
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = nFrames - i*maxFrames
		if framesRemaining > maxFrames:
			framesWork = maxFrames
		else:
			framesWork = nFrames
		print "Working on chunk %i, %i frames remaining" % ((i+1), framesRemaining)

		# Inner loop that actually reads the frames into the data array
		for j in range(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbw.readFrame(fh)
			except errors.eofError:
				break
			except errors.syncError:
				print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbw.FrameSize-1)
				continue
			if not cFrame.header.isTBW():
				continue
			
			stand = cFrame.header.parseID()
			# In the current configuration, stands start at 1 and go up to 10.  So, we
			# can use this little trick to populate the data array
			aStand = 2*(stand-1)
			if cFrame.header.frameCount % 10000 == 0:
				print "%3i -> %3i  %5i  %i" % (stand, aStand, cFrame.header.frameCount, cFrame.data.timeTag)

			# Actually load the data.  x pol goes into the even numbers, y pol into the 
			# odd numbers
			count = cFrame.header.frameCount - 1
			timeTags[aStand,   count] = cFrame.data.timeTag
			timeTags[aStand+1, count] = cFrame.data.timeTag

	# Check for missing frames
	missing = numpy.where( timeTags < 0 )
	if len(missing) != 0:
		dp1Boards = []
		print "Found %i missing frames (%i missing time tags).  Missing data from:" % (len(missing[0])/2, len(missing[0]))
		for i,f in zip(missing[0], missing[1]):
			if antennas[i].board not in dp1Boards:
				dp1Boards.append(antennas[i].board)

			print "  stand %3i, pol. %1i (dig. %3i) @ frame %5i" % (antennas[i].stand.id, antennas[i].pol, antennas[i].digitizer, f+1)
		print "-> DP1 boards with missing frames: %s" % ','.join([str(b) for b in dp1Boards])

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

