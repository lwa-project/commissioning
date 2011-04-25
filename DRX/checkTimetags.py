#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import ephem

from lsl import astro
from lsl.reader import drx
from lsl.reader import errors
from lsl.common.dp import fS

def main(args):
	fh = open(args[0], "rb", buffering=drx.FrameSize*10000)

	# Get the first frame and find out what the firt time tag is, which the
	# first frame number is, and what the sample rate it.  From the sample 
	# rate, estimate how the time tag should advance between frames.
	junkFrame = drx.readFrame(fh)
	sampleRate = junkFrame.getSampleRate()
	tagSkip = fS / sampleRate * junkFrame.data.iq.shape[0]
	fh.seek(0)

	# Store the information about the first frame and convert the timetag to 
	# an ephem.Date object.
	prevTime = junkFrame.data.timeTag
	prevDate = ephem.Date(astro.unix_to_utcjd(junkFrame.getTime()) - astro.DJD_OFFSET)
	prevFrame = junkFrame.header.frameCount

	# Report on the file
	print "Filename: %s" % os.path.basename(args[0])
	print "Date of first frame: %i -> %s" % (prevTime, str(prevDate))
	print "Sample rate: %i Hz" % sampleRate
	print "Time tag skip per frame: %i" % tagSkip

	k = 0
	while True:
		try:
			currFrame = drx.readFrame(fh)
		except errors.eofError:
			break
		except errors.syncError:
			continue
		
		beam, tune, pol = currFrame.parseID()
		currTime = currFrame.data.timeTag
		currDate = ephem.Date(astro.unix_to_utcjd(currFrame.getTime()) - astro.DJD_OFFSET)
		currFrame = 1 + k / 4
		
		#print k, currFrame, beam, tune, pol, currTime
		#print beam, tune, pol, prevTime, currTime, currTime-prevTime

		if k == 0 or (currFrame % 50000 == 0 and tune == 1 and pol == 0):
			print "At beam %i, tuning %i, pol %i:  frame %i -> %i (%s)" % (beam, tune, pol, currFrame, currTime, currDate)

		if currTime < prevTime:
			print "ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i" % (currTime, currFrame, prevTime, prevFrame)
			print "       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate))
		elif (currTime-prevTime) > tagSkip:
			print "Error: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (currTime, currFrame, prevTime, prevFrame)
			print "       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate))
		elif (currTime-prevTime) == 0 or (currTime-prevTime) == tagSkip:
			pass
		else:
			print "Warning: t.t %i @ frame %i is more or less than expected" % (currTime, currFrame)
			print "         -> difference: %i (%.5f seconds; %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate))
			print "         -> %i %i %i" % (beam, tune, pol)
		
		prevTime = currTime
		prevFrame = currFrame
		k = k + 1

	fh.close()


if __name__ == "__main__":
	main(sys.argv[1:])
