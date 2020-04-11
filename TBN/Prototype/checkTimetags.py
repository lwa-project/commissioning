#!/usr/bin/env python

"""
Check the time tags in a TBN file from the prototype system at the
north arm.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import ephem

from lsl import astro
from lsl.reader import tbn
from lsl.reader import errors
from lsl.common.dp import fS

def main(args):
    fh = open(args[0], "rb", buffering=tbn.FRAME_SIZE*10000)

    # Get the first frame and find out what the firt time tag is, which the
    # first frame number is, and what the sample rate it.  From the sample 
    # rate, estimate how the time tag should advance between frames.
    junkFrame = tbn.read_frame(fh)
    sample_rate = tbn.get_sample_rate(fh)
    tagSkip = fS // sample_rate * junkFrame.payload.data.shape[0]
    fh.seek(0)

    # Store the information about the first frame and convert the timetag to 
    # an ephem.Date object.
    prevTime = junkFrame.payload.timetag
    prevDate = ephem.Date(astro.unix_to_utcjd(sum(junkFrame.time, 0.0)) - astro.DJD_OFFSET)
    prevFrame = junkFrame.header.frame_count

    # Report on the file
    print("Filename: %s" % os.path.basename(args[0]))
    print("Date of first frame: %i -> %s" % (prevTime, str(prevDate)))
    print("Sample rate: %i Hz" % sample_rate)
    print("Time tag skip per frame: %i" % tagSkip)

    k = 0
    while True:
        try:
            currFrame = tbn.read_frame(fh)
        except errors.EOFError:
            break
        except errors.SyncError:
            continue
        
        stand, pol = currFrame.id
        currTime = currFrame.payload.timetag
        currDate = ephem.Date(astro.unix_to_utcjd(sum(currFrame.time, 0.0)) - astro.DJD_OFFSET)
        currFrame = currFrame.header.frame_count

        if k == 0 or (currFrame % 5000 == 0 and stand == 1 and pol == 0):
            print("At stand %i, pol %i:  frame %i -> %s" % (stand, pol, currFrame, currDate))

        if currTime < prevTime:
            print("ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i" % (currTime, currFrame, prevTime, prevFrame))
            print("       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate)))
        if (currTime-prevTime) > tagSkip:
            print("ERROR: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (currTime, currFrame, prevTime, prevFrame))
            print("       -> difference: %i (%.5f seconds); %s" % (currTime-prevTime, float(currTime-prevTime)/fS, str(currDate)))
        
        prevTime = currTime
        prevFrame = currFrame
        k = k + 1

    fh.close()


if __name__ == "__main__":
    main(sys.argv[1:])
