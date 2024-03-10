#!/usr/bin/env python3

"""
Check the time times in a DRX file for flow.  This script should be immune to the 
various DRX cross-tuning time tag issues because it does comparisons on a tuning/
polarization basis.
"""

import os
import sys
import argparse

from lsl import astro
from lsl.reader import drx
from lsl.reader import errors
from lsl.common.dp import fS
from lsl.misc import parser as aph


def main(args):
    skip = args.skip
    fh = open(args.filename, "rb")
    
    # Get the first frame and find out what the firt time tag is, which the
    # first frame number is, and what the sample rate it.  From the sample 
    # rate, estimate how the time tag should advance between frames.
    while True:
        junkFrame = drx.read_frame(fh)
        try:
            sample_rate = junkFrame.sample_rate
            break
        except ZeroDivisionError:
            pass
    tagSkip = int(fS / sample_rate * junkFrame.payload.data.shape[0])
    fh.seek(-drx.FRAME_SIZE, 1)

    # Store the information about the first frame.
    prevTime = junkFrame.payload.timetag
    prevDate = junkFrame.time.datetime
    prevFrame = junkFrame.header.frame_count
    
    # Skip ahead
    fh.seek(int(skip*sample_rate/4096)*4*drx.FRAME_SIZE)

    # Report on the file
    print("Filename: %s" % os.path.basename(args.filename))
    print("Date of first frame: %i -> %s" % (prevTime, str(prevDate)))
    print("Sample rate: %i Hz" % sample_rate)
    print("Time tag skip per frame: %i" % tagSkip)
    if skip != 0:
        print("Skipping ahead %i frames (%.6f seconds)" % (int(skip*sample_rate/4096)*4, int(skip*sample_rate/4096)*4096/sample_rate))
        
    k = 0
    prevTime = [-1, -1, -1, -1]
    prevDate = ['', '', '', '']
    prevNumb = [0, 0, 0, 0]
    
    while True:
        try:
            currFrame = drx.read_frame(fh)
        except errors.EOFError:
            break
        except errors.SyncError:
            currNumb = 1 + k // 4
            
            print("ERROR: invalid frame (sync. word error) @ frame %8i" % currNumb)
            continue
        
        beam, tune, pol = currFrame.id
        rID = 2*(tune-1) + pol
        currTime = currFrame.payload.timetag
        currDate = currFrame.time.datetime
        currNumb = prevNumb[rID] + 1

        if tune == 1 and pol == 0 and currNumb % 50000 == 0:
            print("Beam %i, tune %i, pol %i: frame %8i -> %i (%s)" % (beam, tune, pol, currNumb, currTime, currDate))
            
        if currTime < prevTime[rID] and prevTime[rID] >= 0:
            print("ERROR: t.t. %i @ frame %i < t.t. %i @ frame %i" % (currTime, currNumb, prevTime[rID], prevNumb[rID]))
            print("       -> difference: %i (%.3f frames; %.5f seconds); %s" % (currTime-prevTime[rID], float(currTime-prevTime[rID])/tagSkip, float(currTime-prevTime[rID])/fS, str(currDate)))
            print("       -> beam %i, tuning %i, pol %i" % (beam, tune, pol))
        elif currTime > (prevTime[rID] + tagSkip) and prevTime[rID] >= 0:
            print("ERROR: t.t. %i @ frame %i > t.t. %i @ frame %i + skip" % (currTime, currNumb, prevTime[rID], prevNumb[rID]))
            print("       -> difference: %i (%.3f frames; %.5f seconds); %s" % (currTime-prevTime[rID], float(currTime-prevTime[rID])/tagSkip, float(currTime-prevTime[rID])/fS, str(currDate)))
            print("       -> beam %i, tuning %i, pol %i" % (beam, tune, pol))
        elif currTime < (prevTime[rID] + tagSkip) and prevTime[rID] >= 0:
            print("ERROR: t.t %i @ frame %i < t.t. %i @ frame %i + skip" % (currTime, currNumb, prevTime[rID], prevNumb[rID]))
            print("       -> difference: %i (%.3f frames; %.5f seconds; %s" % (currTime-prevTime[rID], float(currTime-prevTime[rID])/tagSkip, float(currTime-prevTime[rID])/fS, str(currDate)))
            print("       -> beam %i, tuning %i, pol %i" % (beam, tune, pol))
        else:
            pass
            
        k += 1
        prevTime[rID] = currTime
        prevDate[rID] = currDate
        prevNumb[rID] = currNumb
        
    fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a DRX file and check the flow of time', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to check')
    parser.add_argument('-s', '--skip', type=aph.positive_float, default=0.0, 
                        help='skip period in seconds between chunks')
    args = parser.parse_args()
    main(args)
    
