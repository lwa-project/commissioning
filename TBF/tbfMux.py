#!/usr/bin/env python

"""
Given a TBF filles created the the on-line triggering system on ADP, combine 
the files together into a single file that can be used like a standard 
DR-recorded TBF file

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import struct
from collections import deque

from lsl.reader.ldp import TBFFile
from lsl.reader import tbf, errors, buffer


class RawTBFFrame(object):
    """
    Class to help hold and work with a raw (packed) TBF frame.
    """
    
    def __init__(self, contents):
        self.contents = bytearray(contents)
        if len(self.contents) != tbf.FrameSize:
            raise errors.eofError
        if self.contents[0] != 0xDE or self.contents[1] != 0xC0 or self.contents[2] != 0xDE or self.contents[3] != 0x5c:
            raise errors.syncError
            
    def __getitem__(self, key):
        return self.contents[key]
        
    def __setitem__(self, key, value):
        self.contents[key] = value
        
    @property
    def timeTag(self):
        timeTag = 0L
        timeTag |= self.contents[16] << 56
        timeTag |= self.contents[17] << 48
        timeTag |= self.contents[18] << 40
        timeTag |= self.contents[19] << 32
        timeTag |= self.contents[20] << 24
        timeTag |= self.contents[21] << 16
        timeTag |= self.contents[22] <<  8
        timeTag |= self.contents[23]
        return timeTag
        
    @property
    def firstChan(self):
        chan0 = (self.contents[12] << 8) | self.contents[13]
        return chan0


class RawTBFFrameBuffer(buffer.FrameBuffer):
    """
    A sub-type of FrameBuffer specifically for dealing with raw (packed) TBF
    frames.  See :class:`lsl.reader.buffer.FrameBuffer` for a description of 
    how the buffering is implemented.
    
    Keywords:
      chans
        list of start channel numbers to expect data for
    
      nSegments
        number of ring segments to use for the buffer (default is 25)
    
      ReorderFrames
        whether or not to reorder frames returned by get() or flush() by 
        start channel (default is False)
    
    The number of segements in the ring can be converted to a buffer time in 
    seconds:
    
    +----------+--------+
    | Segments |  Time  |
    +----------+--------+
    |    10    | 0.0004 |
    +----------+--------+
    |    25    | 0.001  |
    +----------+--------+
    |    50    | 0.002  |
    +----------+--------+
    |   100    | 0.004  |
    +----------+--------+
    |   200    | 0.008  |
    +----------+--------+
    
    """
    
    def __init__(self, chans, nSegments=25, ReorderFrames=False):
        super(RawTBFFrameBuffer, self).__init__(mode='TBF', chans=chans, nSegments=nSegments, ReorderFrames=ReorderFrames)
        
    def calcFrames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.chans)
        for chans in self.chans:
            frameList.append(chans)
            
        return (nFrames, frameList)
        
    def figureOfMerit(self, frame):
        """
        Figure of merit for sorting frames.  For TBF this is:
        frame.data.timeTag
        """
        
        return frame.timeTag
        
    def createFill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = copy.deepcopy(self.buffer[key][0])
        
        # Get out the frame parameters and fix-up the header
        chan = frameParameters
        fillFrame[12] = (chan & 0xFF00) >> 8
        filLFrame[14] = (chan & 0x00FF)
        
        # Zero the data for the fill packet
        fillFrame[32:] = '\x00'*(12*256*2)
        
        return fillFrame


def main(args):
    # Get the files
    filenames = args
    
    # Open them up and make sure we have a continuous range of frequencies
    idf = [TBFFile(filename) for filename in filenames]
    chans = []
    for i in idf:
        chans.extend( i.buffer.chans )
    chans.sort()
    for i in xrange(1, len(chans)):
        if chans[i] != chans[i-1] + 12:
            raise RuntimeError("Unexpected channel increment: %i != 12" % (chans[i]-chans[i-1],))
            
    # Setup the buffer
    buffer = RawTBFFrameBuffer(chans=chans, ReorderFrames=False)
    
    # Setup the output filename
    outname = os.path.basename(filenames[0])
    outname = outname.rsplit('_', 1)[1]
    oh = open(outname, 'wb')
    
    # Go!
    fh = [i.fh for i in idf]
    eofFound = [False for i in idf]
    while not all(eofFound):
        ## Read in a frame from all input files
        rFrames = deque()
        for i,f in enumerate(fh):
            try:
                rFrames.append( RawTBFFrame(f.read(tbf.FrameSize)) )
            except errors.eofError:
                eofFound[i] = True
                continue
            except errors.syncError:
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
    main(sys.argv[1:])
    