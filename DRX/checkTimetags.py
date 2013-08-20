#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Check the time times in a DRX file for flow.  This script should be immune to the 
various DRX cross-tuning time tag issues because it does comparisons on a tuning/
polarization basis.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import ephem
import gc

from lsl import astro
from lsl.reader import drx
from lsl.reader import errors
from lsl.common.dp import fS

def main(args):
	if args[0] == '-s':
		skip = float(args[1])
		fh = open(args[2], "rb")
	else:
		skip = 0
		fh = open(args[0], "rb")

	# Get the first frame and find out what the firt time tag is, which the
	# first frame number is, and what the sample rate it.  From the sample 
	# rate, estimate how the time tag should advance between frames.
	while True:
		junkFrame = drx.readFrame(fh)
		try:
			sampleRate = junkFrame.getSampleRate()
			break
		except ZeroDivisionError:
			pass
	tagSkip = int(fS / sampleRate * junkFrame.data.iq.shape[0])
	fh.seek(-drx.FrameSize, 1)

	# Store the information about the first frame and convert the timetag to 
	# an ephem.Date object.
	prevTime = junkFrame.data.timeTag
	prevDate = ephem.Date(astro.unix_to_utcjd(junkFrame.getTime()) - astro.DJD_OFFSET)
	prevFrame = junkFrame.header.frameCount

	# Skip ahead
	fh.seek(int(skip*sampleRate/4096)*4*drx.FrameSize)

	# Report on the file
	print "Filename: %s" % os.path.basename(args[0])
	print "Date of first frame: %i -> %s" % (prevTime, str(prevDate))
	print "Sample rate: %i Hz" % sampleRate
	print "Time tag skip per frame: %i" % tagSkip
	if skip != 0:
		print "Skipping ahead %i frames (%.6f seconds)" % (int(skip*sampleRate/4096)*4, int(skip*sampleRate/4096)*4096/sampleRate)

	k = 0
	#k = 1
	prevTime = [0, 0, 0, 0]
	prevDate = ['', '', '', '']
	prevNumb = [0, 0, 0, 0]
	for i in xrange(4):
		currFrame = drx.readFrame(fh)
		beam, tune, pol = currFrame.parseID()
		rID = 2*(tune-1) + pol

		prevTime[rID] = currFrame.data.timeTag
		prevDate[rID] = ephem.Date(astro.unix_to_utcjd(currFrame.getTime()) - astro.DJD_OFFSET)
		prevNumb[rID] = 1 + k / 4
		#prevNumb[rID] = k
		
		k += 1
	
	while True:
		try:
			currFrame = drx.readFrame(fh)
		except errors.eofError:
			break
		except errors.syncError:
			currNumb = 1 + k / 4
			
			print "ERROR: invalid frame (sync. word error) @ frame %8i" % currNumb
			continue
		
		beam, tune, pol = currFrame.parseID()
		rID = 2*(tune-1) + pol
		currTime = currFrame.data.timeTag
		currDate = ephem.Date(astro.unix_to_utcjd(currFrame.getTime()) - astro.DJD_OFFSET)
		currNumb = 1 + k / 4
		#currNumb = k

		if tune == 1 and pol == 0 and currNumb % 50000 == 0:
			print "Beam %i, tune %i, pol %i: frame %8i -> %i (%s)" % (beam, tune, pol, currNumb, currTime, currDate)

		if currTime < prevTime[rID]:
			print "ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i" % (currTime, currNumb, prevTime[rID], prevNumb[rID])
			print "       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime[rID], float(currTime-prevTime[rID])/fS, str(currDate))
			print "       -> beam %i, tuning %i, pol %i" % (beam, tune, pol)
		elif currTime > (prevTime[rID] + tagSkip):
			print "ERROR: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (currTime, currNumb, prevTime[rID], prevNumb[rID])
			print "       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime[rID], float(currTime-prevTime[rID])/fS, str(currDate))
			print "       -> beam %i, tuning %i, pol %i" % (beam, tune, pol)
		elif currTime < (prevTime[rID] + tagSkip):
			print "ERROR: t.t %i @ frame %i < t.t. %i @ frame %i + skip" % (currTime, currNumb, prevTime[rID], prevNumb[rID])
			print "       -> difference: %i (%.5f seconds; %s" % (currTime-prevTime[rID], float(currTime-prevTime[rID])/fS, str(currDate))
			print "       -> beam %i, tuning %i, pol %i" % (beam, tune, pol)
		else:
			pass
		
		prevTime[rID] = currTime
		prevDate[rID] = currDate
		prevNumb[rID] = currNumb
		k += 1
		
		del currFrame
		if k % 100 == 0:
			gc.collect()

	fh.close()


if __name__ == "__main__":
	main(sys.argv[1:])
