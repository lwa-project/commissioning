#!/usr/bin/env python

"""
Given a COR filles created by ADP, combine the files together into a single 
file that can be used like a standard DR-recorded COR file
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import copy
import struct
import argparse
from collections import deque

from lsl.reader.ldp import CORFile
from lsl.reader import cor, errors, buffer


class RawCORFrame(object):
    """
    Class to help hold and work with a raw (packed) COR frame.
    """
    
    def __init__(self, contents):
        self.contents = bytearray(contents)
        if len(self.contents) != cor.FRAME_SIZE:
            raise errors.EOFError
        if self.contents[0] != 0xDE or self.contents[1] != 0xC0 or self.contents[2] != 0xDE or self.contents[3] != 0x5c:
            raise errors.SyncError
            
    def __getitem__(self, key):
        return self.contents[key]
        
    def __setitem__(self, key, value):
        self.contents[key] = value
        
    @property
    def timetag(self):
        timetag = 0
        timetag |= self.contents[16] << 56
        timetag |= self.contents[17] << 48
        timetag |= self.contents[18] << 40
        timetag |= self.contents[19] << 32
        timetag |= self.contents[20] << 24
        timetag |= self.contents[21] << 16
        timetag |= self.contents[22] <<  8
        timetag |= self.contents[23]
        return timetag
        
    @property
    def first_chan(self):
        chan0 = (self.contents[12] << 8) | self.contents[13]
        return chan0
        
    @property
    def stand0(self):
        stand = (self.contents[28] << 8) | self.contents[29]
        return stand
        
    @property
    def stand1(self):
        stand = (self.contents[30] << 8) | self.contents[31]
        return stand
        
    def parse_id(self):
        return (self.stand0, self.stand1)


class RawCORFrameBuffer(buffer.FrameBufferBase):
    """
    A sub-type of FrameBufferBase specifically for dealing with raw (packed) COR
    frames.  See :class:`lsl.reader.buffer.FrameBufferBase` for a description of 
    how the buffering is implemented.
    
    Keywords:
      stands
        list of stands to expect baseline data for
    
    chans
        list of start channel numbers to expect data for
    
      nsegments
        number of ring segments to use for the buffer (default is 5)
    
      reorder
        whether or not to reorder frames returned by get() or flush() by 
        start channel (default is False)
    
    The number of segements in the ring can be converted to a buffer time in 
    seconds:
    
    +----------+--------+
    | Segments |  Time  |
    +----------+--------+
    |     1    |   10   |
    +----------+--------+
    |     2    |   20   |
    +----------+--------+
    |     5    |   50   |
    +----------+--------+
    
    """
    
    def __init__(self, chans, nsegments=5, reorder=False):
        super(RawCORFrameBuffer, self).__init__(mode='COR', stands=list(range(1,256+1)), chans=chans, nsegments=nsegments, reorder=reorder)
        
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.stands)*(len(self.stands)+1)/2 * len(self.chans)
        for stand0 in self.stands:
            for stand1 in self.stands:
                if stand1 < stand0:
                    continue
                for chan in self.chans:
                    frameList.append((stand0,stand1,chan))
                    
        return (nFrames, frameList)
        
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for sorting frames.  For COR this is:
        frame.payload.timetag
        """
        
        return frame.timetag
        
    def frameID(self, frame):
        """
        ID value or tuple for a given frame.
        """
        
        return frame.id+(frame.first_chan,)
        
    def createFill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = RawCORFrame( copy.deepcopy(self.buffer[key][0].contents) )
        
        # Get out the frame parameters and fix-up the header
        stand0, stand1, chan = frameParameters
        fillFrame[12] = (chan & 0xFF00) >> 8
        fillFrame[13] = (chan & 0x00FF)
        fillFrame[28] = (stand0 & 0xFF00) >> 8
        fillFrame[29] = (stand0 & 0x00FF)
        fillFrame[30] = (stand1 & 0xFF00) >> 8
        fillFrame[31] = (stand1 & 0x00FF)
        
        # Zero the data for the fill packet
        fillFrame[32:] = '\x00'*(72*4*4)
        
        return fillFrame


def main(args):
    # Parse the command line
    filenames = args.filename
    filenames.sort()
    if len(filenames) < 2:
        raise RuntimeError("Need at least two files to combine")
        
    # Open them up and make sure we have a continuous range of frequencies
    idf = [CORFile(filename) for filename in filenames]
    chans = []
    for i in idf:
        chans.extend( i.buffer.chans )
    chans.sort()
    for i in xrange(1, len(chans)):
        if chans[i] != chans[i-1] + cor.FRAME_CHANNEL_COUNT:
            raise RuntimeError("Unexpected channel increment: %i != %s" % (chans[i]-chans[i-1], cor.FRAME_CHANNEL_COUNT))
            
    # Setup the buffer
    buffer = RawCORFrameBuffer(chans=chans, reorder=False)
    
    # Setup the output filename
    if args.output is None:
        names = [os.path.basename(filename) for filename in filenames]
        common = names[0][-1]
        
        valid = True
        while valid and len(common) < len(names[0]):
            for name in names:
                if name[-len(common):] != common:
                    valid = False
                    break
            if valid:
                common = name[-len(common)-1:]
        common = common[1:]
        if common[0] == '_':
            common = common[1:]
        args.output = common
        
    print("Writing combined file to '%s'" % os.path.basename(args.output))
    oh = open(args.output, 'wb')
    
    # Go!
    fh = [i.fh for i in idf]
    eofFound = [False for i in idf]
    while not all(eofFound):
        ## Read in a frame from all input files
        rFrames = deque()
        for i,f in enumerate(fh):
            try:
                rFrames.append( RawCORFrame(f.read(cor.FRAME_SIZE)) )
            except errors.EOFError:
                eofFound[i] = True
                continue
            except errors.SyncError:
                continue
                
        ## Add the frames to the buffer
        buffer.append(rFrames)
        rFrames = buffer.get()
        
        ## Continue adding frames if nothing comes out.
        if rFrames is None:
            continue
            
        ## Write the re-ordered frames to the output file
        for rFrame in rFrames:
            oh.write(rFrame.contents)
    # Empty the buffer
    for rFrames in buffer.flush():
        for rFrame in rFrames:
            oh.write(rFrame.contents)
            
    # Done
    oh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='given a COR files created by ADP, combine the files together into a single file that can be used like a standard DR-recorded COR file', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+', 
                        help='filename to combine')
    parser.add_argument('-o', '--output', type=str, 
                        help='write the combined file to the provided filename, auto-determine if not provided')
    args = parser.parse_args()
    main(args)
    