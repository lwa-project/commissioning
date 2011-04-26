#!/usr/bin/env python

import os
import sys
import numpy

from lsl.common.dp import fS
from lsl.reader import drx
from lsl.reader import errors

def main(args):
	filename = args[0]

	fh = open(filename, "rb")
	nFramesFile = os.path.getsize(filename) / drx.FrameSize
	junkFrame = drx.readFrame(fh)
	beam,tune,pol = junkFrame.parseID()
	while 2*(tune-1)+pol != 0:
		junkFrame = drx.readFrame(fh)
		beam,tune,pol = junkFrame.parseID()
	fh.seek(fh.tell() - drx.FrameSize)

	srate = junkFrame.getSampleRate()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol

	# File summary
	print "Filename: %s" % filename
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "==="

	beamIDs = [0,0,0,0]
	timeTags = numpy.zeros(4, dtype=numpy.int64) - 1
	timeOffsets = numpy.zeros(4, dtype=numpy.int64) - 1
	timeValues = numpy.zeros(4, dtype=numpy.float64)
	for i in xrange(16):
		# Read in the next frame and anticipate any problems that could occur
		try:
			cFrame = drx.readFrame(fh, Verbose=False)
		except errors.eofError:
			break
		except errors.syncError:
			#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/drx.FrameSize-1)
			continue
		
		## Save the time time, time offset, and computed time values
		beam,tune,pol = cFrame.parseID()
		aStand = 2*(tune-1) + pol
		if timeTags[aStand] == -1:
			beamIDs[aStand] = (beam,tune,pol)
			timeTags[aStand] = cFrame.data.timeTag
			timeOffsets[aStand] = cFrame.header.timeOffset
			timeValues[aStand] = cFrame.getTime()

	for id,tt,to,tv in zip(beamIDs, timeTags, timeOffsets, timeValues):
		b,t,p = id
		print "B%i, T%i, P%i: t.t. is %i with offset %i -> %.9f" % (b,t,p,tt,to,tv)

	print "==="
	t1t = timeTags[0] - timeOffsets[0]
	t2t = timeTags[3] - timeOffsets[3]
	if (t2t-t1t) != 0:
		print "T2 time tags appear to be offset from T1 by %i (%.9f s)" % (t2t-t1t, (t2t-t1t)/fS)
	else:
		print "T2 and T1 do not appear to have a time tag offset"


if __name__ == "__main__":
	main(sys.argv[1:])

