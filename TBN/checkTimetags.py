#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Check the time tags in an old 20-stand prototype system data file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import ephem

from lsl import astro
from lsl.reader import tbn
from lsl.reader import errors
from lsl.common.dp import fS

def main(args):
	fh = open(args[0], "rb", buffering=tbn.FrameSize*10000)

	# Get the first frame and find out what the firt time tag is, which the
	# first frame number is, and what the sample rate it.  From the sample 
	# rate, estimate how the time tag should advance between frames.
	junkFrame = tbn.readFrame(fh)
	sampleRate = tbn.getSampleRate(fh)
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
			currFrame = tbn.readFrame(fh)
		except errors.eofError:
			break
		except errors.syncError:
			continue
		except errors.numpyError:
			break
		
		stand, pol = currFrame.parseID()
		currTime = currFrame.data.timeTag
		currDate = ephem.Date(astro.unix_to_utcjd(currFrame.getTime()) - astro.DJD_OFFSET)
		currFrame = currFrame.header.frameCount

		if k == 0 or (currFrame % 5000 == 0 and stand == 1 and pol == 0):
			print "At stand %i, pol %i:  frame %i -> %s" % (stand, pol, currFrame, currDate)

		if currTime < prevTime:
			print "ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i" % (currTime, currFrame, prevTime, prevFrame)
			print "       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate))
		if (currTime-prevTime) > tagSkip:
			print "Error: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (currTime, currFrame, prevTime, prevFrame)
			print "       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate))
		
		prevTime = currTime
		prevFrame = currFrame
		k = k + 1

	fh.close()


if __name__ == "__main__":
	main(sys.argv[1:])
