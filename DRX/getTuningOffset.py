#!/usr/bin/env python

import os
import sys
import numpy
import curses
import string

from lsl.common.dp import fS
from lsl.reader import drx
from lsl.reader import errors

display = string.Template("""$preamble

B$b0, T$t0, P$p0: t.t. is $tt0 with offset $os0 -> $tv0
B$b1, T$t1, P$p1: t.t. is $tt1 with offset $os1 -> $tv1
B$b2, T$t2, P$p2: t.t. is $tt2 with offset $os2 -> $tv2
B$b3, T$t3, P$p3: t.t. is $tt3 with offset $os3 -> $tv3
-> T2 time tags appear to be offset from T1 by $ttd ($tvd s)
""")

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
	out = "Filename: %s" % filename
	out += "\nBeams: %i" % beams
	out += "\nTune/Pols: %i %i %i %i" % tunepols
	out += "\nSample Rate: %i Hz" % srate
	out += "\nFrames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	out += "\n==="
	print out

	tuningOffset = numpy.zeros(nFramesFile/8, dtype=numpy.int64)
	try:
		screen = curses.initscr()
		curses.noecho()
		curses.cbreak()
		screen.nodelay(1)

		strdict = {'preamble': out}
		for i in xrange(tuningOffset.size):
			screen.clear()

			beamIDs = [0,0,0,0]
			timeTags = numpy.zeros(4, dtype=numpy.int64) - 1
			timeOffsets = numpy.zeros(4, dtype=numpy.int64) - 1
			timeValues = numpy.zeros(4, dtype=numpy.float64)
			for j in xrange(4):
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

			k = 0
			for id,tt,to,tv in zip(beamIDs, timeTags, timeOffsets, timeValues):
				strdict['b%i' % k] = id[0]
				strdict['t%i' % k] = id[1]
				strdict['p%i' % k] = id[2]

				strdict['tt%i' % k] = tt
				strdict['os%i' % k] = to
				strdict['tv%i' % k] = tv

				k += 1

			t1t = timeTags[0] - timeOffsets[0]
			t2t = timeTags[3] - timeOffsets[3]
			tuningOffset[i] = t2t - t1t

			strdict['ttd'] = t2t - t1t
			strdict['tvd'] = (t2t-t1t)/fS

			screen.addstr(0,0,display.safe_substitute(strdict))
			screen.refresh()
				
			# Check for keypress and exit if Q
			c = screen.getch()
			if (c > 0):
				if chr(c) == 'q': 
					break
				if chr(c) == 'Q': 
					break

		curses.nocbreak()
    		curses.echo()
    		curses.endwin()

	except KeyboardInterrupt:
		curses.nocbreak()
		curses.echo()
		curses.endwin()

		tuningOffset = tuningOffset[0:i]

	print display.safe_substitute(strdict)

	print "T2-T1 time tag offset range: %i to %i (based on %i sets of frames)" % (tuningOffset.min(), tuningOffset.max(), len(tuningOffset))

if __name__ == "__main__":
	main(sys.argv[1:])

