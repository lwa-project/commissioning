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
import getopt
import struct
from collections import deque

from lsl.reader.ldp import TBFFile
from lsl.reader import tbf, errors, buffer


def usage(exitCode=None):
    print """tbfMux.py - Given a TBF filles created the the on-line triggering system on ADP, 
combine the files together into a single file that can be used like a standard 
DR-recorded TBF file

Usage: tbfMux.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-o, --output                Write the combined file to the provided filename
                            (Default = auto-deterine the filename)
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['output'] = None
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "ho:", ["help", "output="])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-o', '--output'):
            config['output'] = value
        else:
            assert False
            
    # Add in arguments
    config['args'] = args
    
    # Return configuration
    return config


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
    # Parse the command line
    config = parseOptions(args)
    filenames = config['args']
    
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
    if config['output'] is None:
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
        config['output'] = common
        
    print "Writing combined file to '%s'" % os.path.basename(config['output'])
    oh = open(config['output'], 'wb')
    
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
    