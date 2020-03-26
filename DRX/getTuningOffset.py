#!/usr/bin/env python

"""
Read in a DRX file and look at the time tag difference between tuning 1 and 
tuning 2 to and if that difference changes throughout a file.
"""

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
    nFramesFile = os.path.getsize(filename) / drx.FRAME_SIZE
    junkFrame = drx.read_frame(fh)
    beam,tune,pol = junkFrame.id
    while 2*(tune-1)+pol != 0:
        junkFrame = drx.read_frame(fh)
        beam,tune,pol = junkFrame.id
    fh.seek(fh.tell() - drx.FRAME_SIZE)

    srate = junkFrame.sample_rate
    beams = drx.get_beam_count(fh)
    tunepols = drx.get_frames_per_obs(fh)
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
            timetags = numpy.zeros(4, dtype=numpy.int64) - 1
            time_offsets = numpy.zeros(4, dtype=numpy.int64) - 1
            timeValues = numpy.zeros(4, dtype=numpy.float64)
            for j in xrange(4):
                # Read in the next frame and anticipate any problems that could occur
                try:
                    cFrame = drx.read_frame(fh, verbose=False)
                except errors.EOFError:
                    break
                except errors.SyncError:
                    #print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/drx.FRAME_SIZE-1)
                    continue
        
                ## Save the time time, time offset, and computed time values
                beam,tune,pol = cFrame.id
                aStand = 2*(tune-1) + pol
                if timetags[aStand] == -1:
                    beamIDs[aStand] = (beam,tune,pol)
                    timetags[aStand] = cFrame.data.timetag
                    time_offsets[aStand] = cFrame.header.time_offset
                    timeValues[aStand] = cFrame.get_time()

            k = 0
            for id,tt,to,tv in zip(beamIDs, timetags, time_offsets, timeValues):
                strdict['b%i' % k] = id[0]
                strdict['t%i' % k] = id[1]
                strdict['p%i' % k] = id[2]

                strdict['tt%i' % k] = tt
                strdict['os%i' % k] = to
                strdict['tv%i' % k] = tv

                k += 1

            t1t = timetags[0] - time_offsets[0]
            t2t = timetags[3] - time_offsets[3]
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

